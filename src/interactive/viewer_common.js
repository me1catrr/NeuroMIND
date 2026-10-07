/* NeuroMIND — shared canvas helpers for longitudinal / transversal viewers */
(function (global) {
  'use strict';

  const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"];
  const effectValue = row => Number(
    row && row.effect_dz != null ? row.effect_dz :
      (row && row.effect_d_pooled != null ? row.effect_d_pooled : 0)
  ) || 0;
  const isFdrResult = row => row && row.is_fdr != null
    ? row.is_fdr === true || String(row.is_fdr) === 'true'
    : false;
  const isNominalResult = row => row && row.is_nominal != null
    ? row.is_nominal === true || String(row.is_nominal) === 'true'
    : false;

  function sortBandsPhysio(rows, key) {
    key = key || 'band';
    const order = Object.fromEntries(BAND_ORDER.map((b, i) => [b, i]));
    return [...(rows || [])].sort((a, b) => (order[a[key]] ?? 999) - (order[b[key]] ?? 999));
  }

  function resizeCanvas(cv, opts) {
    opts = opts || {};
    const dpr = window.devicePixelRatio || 1;
    const parent = cv.parentElement;
    let cssW = parent ? parent.clientWidth - 8 : 0;
    if (cssW < 40) cssW = opts.fallbackW || 640;
    let cssH = opts.height || Number(cv.getAttribute('data-h')) || 280;
    if (opts.square) cssH = Math.min(cssW, opts.maxSquare || 420);
    cv.style.width = cssW + 'px';
    cv.style.height = cssH + 'px';
    cv.width = Math.round(cssW * dpr);
    cv.height = Math.round(cssH * dpr);
    const ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    return { ctx, w: cssW, h: cssH, dpr };
  }

  function viridis(t) {
    t = Math.max(0, Math.min(1, t));
    const stops = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
    const x = t * (stops.length - 1), i = Math.floor(x), f = x - i;
    const a = stops[i], b = stops[Math.min(i + 1, stops.length - 1)];
    return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
  }

  /* t in [-1,1]: neg=blue, pos=red (matches Makie :RdBu_r) */
  function rdBu(t) {
    const u = (t + 1) / 2;
    if (u < 0.5) {
      const f = u / 0.5;
      return [33 + f * (247 - 33), 102 + f * (247 - 102), 172 + f * (247 - 172)];
    }
    const f = (u - 0.5) / 0.5;
    return [247 + f * (178 - 247), 247 + f * (24 - 247), 247 + f * (43 - 247)];
  }

  function drawColorbar(ctx, x, y, w, h, lim, diverging) {
    for (let i = 0; i < h; i++) {
      const t = 1 - i / h;
      let col;
      if (diverging) {
        const [r, g, b] = rdBu((t * 2 - 1));
        col = `rgb(${r | 0},${g | 0},${b | 0})`;
      } else {
        const [r, g, b] = viridis(t);
        col = `rgb(${r | 0},${g | 0},${b | 0})`;
      }
      ctx.fillStyle = col;
      ctx.fillRect(x, y + i, w, 1);
    }
    ctx.strokeStyle = '#94a3b8';
    ctx.strokeRect(x, y, w, h);
    ctx.fillStyle = '#64748b';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'left';
    if (diverging) {
      ctx.fillText(('+' + lim.toFixed(3)), x + w + 4, y + 8);
      ctx.fillText('0', x + w + 4, y + h / 2 + 3);
      ctx.fillText(('-' + lim.toFixed(3)), x + w + 4, y + h - 2);
    } else {
      ctx.fillText(lim.toFixed(3), x + w + 4, y + 8);
      ctx.fillText('0', x + w + 4, y + h - 2);
    }
  }

  function drawHeatmap(cv, flat, opts) {
    opts = opts || {};
    const { ctx, w, h } = resizeCanvas(cv, { square: true, maxSquare: opts.maxSquare || 380, fallbackW: 320 });
    if (!flat || !opts.n) return;
    const n = opts.n, ch = opts.channels || [];
    const cbarW = 14, padR = 56, padL = 36, padB = 36, padT = 10;
    const cell = Math.min((w - padL - padR) / n, (h - padT - padB) / n);
    const ox = padL + Math.max(0, (w - padL - padR - cell * n) / 2);
    const oy = padT;
    let lim = opts.lim;
    if (!(lim > 0)) {
      lim = 1e-9;
      for (let k = 0; k < flat.length; k++) lim = Math.max(lim, Math.abs(flat[k]));
    }
    const diverging = !!opts.diverging;
    const fdrSet = opts.fdrSet || null;
    // 'lower' = i>j (below diagonal); 'upper' = i<j; null/full = both
    const tri = opts.triangle || null;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const skip = (tri === 'lower' && i < j) || (tri === 'upper' && i > j);
        if (skip) {
          ctx.fillStyle = '#f1f5f9';
          ctx.fillRect(ox + j * cell, oy + i * cell, cell + 0.5, cell + 0.5);
          continue;
        }
        const v = flat[j * n + i];
        let col;
        if (i === j) col = [226, 232, 240];
        else if (diverging) {
          const [r, g, b] = rdBu(v / lim);
          col = [r, g, b];
        } else {
          const [r, g, b] = viridis(Math.max(0, v) / lim);
          col = [r, g, b];
        }
        ctx.fillStyle = `rgb(${col[0] | 0},${col[1] | 0},${col[2] | 0})`;
        ctx.fillRect(ox + j * cell, oy + i * cell, cell + 0.5, cell + 0.5);
        if (fdrSet && i !== j) {
          const a = ch[i], b = ch[j];
          const key = a < b ? a + '|' + b : b + '|' + a;
          if (fdrSet.has(key)) {
            ctx.strokeStyle = '#0f172a';
            ctx.lineWidth = 1.2;
            ctx.strokeRect(ox + j * cell + 0.5, oy + i * cell + 0.5, cell - 1, cell - 1);
          }
        }
      }
    }
    ctx.fillStyle = '#64748b';
    ctx.font = Math.max(7, Math.min(9, cell * 0.55)) + 'px sans-serif';
    for (let i = 0; i < n; i++) {
      ctx.save();
      ctx.translate(ox + i * cell + cell / 2, oy + n * cell + 10);
      ctx.rotate(0.7);
      ctx.fillText(ch[i] || '', 0, 0);
      ctx.restore();
      ctx.textAlign = 'right';
      ctx.fillText(ch[i] || '', ox - 3, oy + i * cell + cell * 0.7);
      ctx.textAlign = 'left';
    }
    drawColorbar(ctx, ox + n * cell + 8, oy, cbarW, Math.max(40, n * cell), lim, diverging);
    cv._hm = { ox, oy, cell, n, ch, flat, diverging, lim, triangle: tri };
  }

  function bindHeatmap(cv, edgeLookupFn) {
    cv.addEventListener('mousemove', ev => {
      const info = cv._hm;
      if (!info) return;
      const rect = cv.getBoundingClientRect();
      const x = ev.clientX - rect.left;
      const y = ev.clientY - rect.top;
      const j = Math.floor((x - info.ox) / info.cell);
      const i = Math.floor((y - info.oy) / info.cell);
      if (i < 0 || j < 0 || i >= info.n || j >= info.n) {
        hideTip();
        return;
      }
      if ((info.triangle === 'lower' && i < j) || (info.triangle === 'upper' && i > j)) {
        hideTip();
        return;
      }
      const a = info.ch[i], b = info.ch[j], v = info.flat[j * info.n + i];
      if (i === j) {
        showTip(ev.clientX, ev.clientY, a + ' (diag)');
        return;
      }
      const ed = edgeLookupFn ? edgeLookupFn(a, b) : null;
      let t = `${a}–${b} = ${Number(v).toFixed(4)}`;
      if (ed) {
        const isLong = ed.t1_mean != null || ed.t2_mean != null;
        if (isLong) {
          t += `\nT1=${Number(ed.t1_mean).toFixed(4)}  T2=${Number(ed.t2_mean).toFixed(4)}` +
               `  Δ=${Number(ed.diff).toFixed(4)}`;
          t += `\np=${Number(ed.p_value).toFixed(4)}  q=${Number(ed.q_value).toFixed(4)}  dz=${effectValue(ed).toFixed(3)}`;
          if (ed.n != null) t += `\nn=${ed.n}`;
        } else {
          t += `\nCtrl=${Number(ed.ctrl_mean).toFixed(4)}  MS=${Number(ed.ms_mean).toFixed(4)}`;
          t += `\np=${Number(ed.p_value).toFixed(4)}  q=${Number(ed.q_value).toFixed(4)}  d=${effectValue(ed).toFixed(3)}`;
          if (ed.effect_rrb != null && ed.effect_rrb !== '' && !isNaN(Number(ed.effect_rrb))) {
            t += `  r_rb=${Number(ed.effect_rrb).toFixed(3)}`;
          }
          if (ed.n_ms != null) t += `\nn_MS=${ed.n_ms}  n_Ctrl=${ed.n_ctrl != null ? ed.n_ctrl : '—'}`;
        }
      }
      showTip(ev.clientX, ev.clientY, t);
    });
    cv.addEventListener('mouseleave', hideTip);
  }

  function showTip(x, y, txt) {
    const t = document.getElementById('tip');
    if (!t) return;
    t.style.display = 'block';
    t.style.left = (x + 12) + 'px';
    t.style.top = (y + 12) + 'px';
    t.textContent = txt;
  }
  function hideTip() {
    const t = document.getElementById('tip');
    if (t) t.style.display = 'none';
  }

  function signedCandidate(values, cis, opts) {
    opts = opts || {};
    if (opts.signed) return true;
    if (values.some(v => v < 0)) return true;
    return !!(cis && cis.some(function (c) { return c && (c[0] < 0 || c[1] < 0); }));
  }
  function niceCeil(x) {
    if (!(x > 0)) return 1;
    const exp = Math.pow(10, Math.floor(Math.log10(x)));
    const f = x / exp;
    const nf = f <= 1 ? 1 : f <= 1.5 ? 1.5 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 3 ? 3 : f <= 4 ? 4 : f <= 5 ? 5 : 10;
    return nf * exp;
  }

  function barChart(cv, labels, values, colorFn, opts) {
    opts = opts || {};
    const { ctx, w, h } = resizeCanvas(cv, {
      height: opts.height || Number(cv.getAttribute('data-h')) || 220,
      fallbackW: 640
    });
    if (!labels.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin datos', 20, 40);
      return;
    }
    const pad = { l: opts.padL || 48, r: 14, t: opts.padT || (opts.yLabel ? 28 : 22), b: 58 };
    const cis = opts.cis || null; // optional [[lo,hi], ...]
    let vmax = Math.max(...values.map(Math.abs), 1e-9);
    if (cis) {
      cis.forEach(function (ci) {
        if (!ci) return;
        vmax = Math.max(vmax, Math.abs(ci[0] || 0), Math.abs(ci[1] || 0));
      });
    }
    // Nice symmetric limit for signed charts
    if (opts.niceLimit && (signedCandidate(values, cis, opts))) {
      vmax = niceCeil(vmax);
    }
    const bw = (w - pad.l - pad.r) / labels.length;
    const signed = signedCandidate(values, cis, opts);
    const y0 = pad.t + (h - pad.t - pad.b) / 2;
    const base = signed ? y0 : (h - pad.b);
    const scale = signed ? (h - pad.t - pad.b) / 2 / vmax : (h - pad.t - pad.b) / vmax;
    const showVals = opts.showValues !== false;
    const showAxes = opts.showAxes !== false;
    const nTicks = opts.nTicks || (signed ? 5 : 3);

    if (opts.yLabel) {
      ctx.fillStyle = '#334155';
      ctx.font = '11px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText(opts.yLabel, w / 2, 14);
    }

    if (showAxes) {
      ctx.strokeStyle = '#e2e8f0';
      ctx.fillStyle = '#94a3b8';
      ctx.font = '9px sans-serif';
      ctx.textAlign = 'right';
      const ticks = [];
      if (signed) {
        for (let i = 0; i < nTicks; i++) {
          ticks.push(-vmax + (2 * vmax) * i / (nTicks - 1));
        }
      } else {
        for (let i = 0; i < nTicks; i++) {
          ticks.push(vmax * i / (nTicks - 1));
        }
      }
      ticks.forEach(tv => {
        const yy = signed
          ? y0 - tv * scale
          : (h - pad.b) - tv * scale;
        ctx.beginPath();
        ctx.moveTo(pad.l, yy);
        ctx.lineTo(w - pad.r, yy);
        ctx.stroke();
        const lab = Math.abs(tv) >= 10 ? tv.toFixed(0) : (Math.abs(tv) >= 1 ? tv.toFixed(2) : tv.toFixed(2));
        ctx.fillText(lab, pad.l - 4, yy + 3);
      });
    }

    ctx.strokeStyle = '#cbd5e1';
    ctx.beginPath();
    if (signed) {
      ctx.moveTo(pad.l, y0);
      ctx.lineTo(w - pad.r, y0);
    } else {
      ctx.moveTo(pad.l, h - pad.b);
      ctx.lineTo(w - pad.r, h - pad.b);
    }
    ctx.stroke();

    labels.forEach((lab, i) => {
      const v = values[i];
      const bh = v * scale;
      const x = pad.l + i * bw + bw * 0.15;
      const xc = x + bw * 0.35;
      ctx.fillStyle = colorFn(v, i);
      if (signed) ctx.fillRect(x, v >= 0 ? y0 - bh : y0, bw * 0.7, Math.abs(bh));
      else ctx.fillRect(x, base - Math.abs(bh), bw * 0.7, Math.abs(bh));

      if (cis && cis[i]) {
        const lo = cis[i][0], hi = cis[i][1];
        const yLo = signed ? y0 - lo * scale : (h - pad.b) - lo * scale;
        const yHi = signed ? y0 - hi * scale : (h - pad.b) - hi * scale;
        ctx.strokeStyle = '#0f172a';
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.moveTo(xc, yLo);
        ctx.lineTo(xc, yHi);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(xc - 4, yLo);
        ctx.lineTo(xc + 4, yLo);
        ctx.moveTo(xc - 4, yHi);
        ctx.lineTo(xc + 4, yHi);
        ctx.stroke();
      }

      if (showVals) {
        ctx.fillStyle = '#334155';
        ctx.font = '9px sans-serif';
        ctx.textAlign = 'center';
        const txt = Math.abs(v) >= 10 ? v.toFixed(0) : (Math.abs(v) >= 1 ? v.toFixed(2) : v.toFixed(3));
        let ty = signed
          ? (v >= 0 ? y0 - Math.abs(bh) - 4 : y0 + Math.abs(bh) + 10)
          : (base - Math.abs(bh) - 4);
        if (cis && cis[i]) {
          const yHi = signed ? y0 - Math.max(cis[i][0], cis[i][1]) * scale : (h - pad.b) - Math.max(cis[i][0], cis[i][1]) * scale;
          if (signed && v >= 0) ty = Math.min(ty, yHi - 4);
        }
        ctx.fillText(txt, xc, ty);
      }

      ctx.save();
      ctx.translate(pad.l + i * bw + bw / 2, h - 10);
      ctx.rotate(-0.45);
      ctx.fillStyle = '#475569';
      ctx.font = '10px sans-serif';
      ctx.textAlign = 'right';
      const short = String(lab).replace('BETA_', 'β').replace('DELTA', 'δ')
        .replace('THETA', 'θ').replace('ALPHA', 'α').replace('GAMMA', 'γ')
        .replace('LOW', 'L').replace('MID', 'M').replace('HIGH', 'H');
      ctx.fillText(short, 0, 0);
      ctx.restore();
    });
  }

  /**
   * Topographic graph of filtered edges.
   * opts.mode: 'fdr' | 'top_d' | 'p05'
   * FDR edges (q<0.05): vivid red/blue; exploratory non-FDR: amber/slate.
   */
  function drawGraph(cv, S, edges, opts) {
    opts = opts || {};
    const { ctx, w, h } = resizeCanvas(cv, { square: true, maxSquare: 440, fallbackW: 420 });
    const cx = w / 2, cy = h / 2 - (opts.legend !== false ? 8 : 0);
    const R = Math.min(w, h) * 0.36;
    ctx.strokeStyle = '#cbd5e1';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.arc(cx, cy, R * 1.08, 0, Math.PI * 2);
    ctx.stroke();
    const pos = S.positions || {};
    const xy = {};
    (S.channels || []).forEach(c => {
      const p = pos[c.toUpperCase()] || pos[c] || [0, 0];
      xy[c] = { x: cx + p[0] * R, y: cy - p[1] * R };
    });
    const list = edges || [];
    const thickKey = opts.thicknessKey || 'diff'; // 'diff' | 'effect_dz' | 'effect_d_pooled'
    const mag = e => Math.abs(
      thickKey === 'effect_dz' || thickKey === 'effect_d_pooled'
        ? effectValue(e) : Number(e.diff) || 0
    );
    const dmax = Math.max(1e-9, ...list.map(mag));

    // Draw non-FDR first so FDR sits on top
    const sorted = [...list].sort((a, b) => {
      const af = isFdrResult(a) ? 1 : 0;
      const bf = isFdrResult(b) ? 1 : 0;
      return af - bf;
    });

    sorted.forEach(e => {
      const A = xy[e.ch_a], B = xy[e.ch_b];
      if (!A || !B) return;
      const d = Number(e.diff) || 0;
      const isFdr = isFdrResult(e);
      const tmag = mag(e) / dmax;
      if (isFdr) {
        ctx.strokeStyle = d >= 0 ? 'rgba(185,28,28,0.88)' : 'rgba(29,78,216,0.88)';
        ctx.lineWidth = 1.8 + 3.2 * tmag;
        ctx.setLineDash([]);
      } else {
        ctx.strokeStyle = d >= 0 ? 'rgba(194,65,12,0.55)' : 'rgba(71,85,105,0.5)';
        ctx.lineWidth = 1.0 + 2.0 * tmag;
        ctx.setLineDash([4, 3]);
      }
      ctx.beginPath();
      ctx.moveTo(A.x, A.y);
      ctx.lineTo(B.x, B.y);
      ctx.stroke();
      ctx.setLineDash([]);
    });

    (S.channels || []).forEach(c => {
      const p = xy[c];
      if (!p) return;
      ctx.beginPath();
      ctx.arc(p.x, p.y, 5, 0, Math.PI * 2);
      ctx.fillStyle = '#fff';
      ctx.fill();
      ctx.strokeStyle = '#0f172a';
      ctx.lineWidth = 1.2;
      ctx.stroke();
      ctx.fillStyle = '#0f172a';
      ctx.font = '9px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText(c, p.x, p.y - 8);
    });

    if (!list.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.font = '13px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('Sin edges con el filtro actual', cx, cy - 8);
      ctx.font = '11px sans-serif';
      ctx.fillText('Prueba Top-N |d| para exploración', cx, cy + 12);
    } else if (opts.legend !== false) {
      const ly = h - 18;
      ctx.font = '10px sans-serif';
      ctx.textAlign = 'left';
      ctx.setLineDash([]);
      ctx.strokeStyle = 'rgba(185,28,28,0.9)';
      ctx.lineWidth = 3;
      ctx.beginPath(); ctx.moveTo(16, ly); ctx.lineTo(36, ly); ctx.stroke();
      ctx.fillStyle = '#334155';
      ctx.fillText(opts.legendFdr || 'FDR q<0.05', 40, ly + 3);
      ctx.strokeStyle = 'rgba(194,65,12,0.7)';
      ctx.lineWidth = 2;
      ctx.setLineDash([4, 3]);
      ctx.beginPath(); ctx.moveTo(120, ly); ctx.lineTo(140, ly); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillText(opts.legendExplore || 'Exploratorio', 144, ly + 3);
      ctx.fillStyle = '#64748b';
      const thickLab = thickKey === 'effect_dz' || thickKey === 'effect_d_pooled'
        ? 'grosor ∝ |dz|' : 'grosor ∝ |Δ|';
      ctx.fillText(opts.legendThick || thickLab, w - 100, ly + 3);
      if (opts.legendDir) {
        ctx.fillText(opts.legendDir, 16, ly - 12);
      }
    }
  }

  function drawTopo(cv, rows, positions, opts) {
    opts = opts || {};
    const showN = opts.showN !== false;
    const { ctx, w, h } = resizeCanvas(cv, { square: true, maxSquare: 430, fallbackW: 430 });
    if (!rows || !rows.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin datos espectrales', 20, 40);
      return;
    }
    const plotW = Math.max(220, w - 72);
    const cx = plotW / 2 + 6, cy = h / 2, R = Math.min(plotW, h) * 0.39;
    ctx.strokeStyle = '#94a3b8';
    ctx.beginPath();
    ctx.arc(cx, cy, R * 1.05, 0, Math.PI * 2);
    ctx.stroke();
    const pos = positions || {};
    const lim = Math.max(1e-9, ...rows.map(r => Math.abs(Number(r.diff) || 0)));
    const pts = [];
    rows.forEach(r => {
      const c = String(r.channel);
      const p = pos[c.toUpperCase()] || pos[c] || [0, 0];
      const x = cx + p[0] * R, y = cy - p[1] * R;
      const [rr, gg, bb] = rdBu((Number(r.diff) || 0) / lim);
      const sig = isFdrResult(r);
      const n = r.n != null ? Number(r.n) : NaN;
      const coverage = r.coverage_pct != null ? Number(r.coverage_pct) : 100;
      const lowN = r.coverage_low === true;
      const alpha = Math.max(0.35, Math.min(1, coverage / 100));
      const rad = sig ? 10 : (lowN ? 8 : 7);
      if (lowN) {
        ctx.setLineDash([3, 2]);
        ctx.strokeStyle = '#a16207';
        ctx.lineWidth = 1.6;
        ctx.beginPath();
        ctx.arc(x, y, rad + 3, 0, Math.PI * 2);
        ctx.stroke();
        ctx.setLineDash([]);
      }
      ctx.beginPath();
      ctx.arc(x, y, rad, 0, Math.PI * 2);
      ctx.fillStyle = `rgba(${rr | 0},${gg | 0},${bb | 0},${alpha.toFixed(3)})`;
      ctx.fill();
      ctx.strokeStyle = sig ? '#0f172a' : (lowN ? '#a16207' : '#94a3b8');
      ctx.lineWidth = sig ? 2.4 : 1;
      ctx.stroke();
      ctx.fillStyle = '#0f172a';
      ctx.font = '10px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText(c, x, y - 11);
      if (showN && !isNaN(n)) {
        ctx.fillStyle = lowN ? '#a16207' : '#64748b';
        ctx.font = '9px sans-serif';
        ctx.fillText('n=' + n, x, y + rad + 10);
      }
      pts.push({ x, y, c, n, coverage, diff: Number(r.diff),
        p: Number(r.p_value), q: Number(r.q_value), lowN });
    });
    drawColorbar(ctx, w - 62, 20, 12, h - 56, lim, true);
    ctx.fillStyle = '#64748b';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText('Opacidad = cobertura · ámbar = cobertura <70% · negro = q<0,05', 12, h - 8);

    cv._topo = { pts };
    if (!cv._topoBound) {
      cv._topoBound = true;
      cv.addEventListener('mousemove', ev => {
        const info = cv._topo;
        if (!info) return;
        const rect = cv.getBoundingClientRect();
        const mx = ev.clientX - rect.left, my = ev.clientY - rect.top;
        let best = null, bd = 1e9;
        info.pts.forEach(p => {
          const dx = p.x - mx, dy = p.y - my;
          const d2 = dx * dx + dy * dy;
          if (d2 < bd) { bd = d2; best = p; }
        });
        if (best && bd < 140) {
          let t = best.c + '\nΔ=' + best.diff.toFixed(4);
          if (!isNaN(best.n)) t += '  n válido=' + best.n;
          if (!isNaN(best.coverage)) t += '  cobertura=' + best.coverage.toFixed(1) + '%';
          if (!isNaN(best.p)) t += '\np=' + best.p.toFixed(4) + '  q=' + best.q.toFixed(4);
          if (best.lowN) t += '\n⚠ cobertura <70%';
          showTip(ev.clientX, ev.clientY, t);
        } else hideTip();
      });
      cv.addEventListener('mouseleave', hideTip);
    }
  }

  function drawStrength(cv, ptsA, ptsB, labelA, labelB) {
    const { ctx, w, h } = resizeCanvas(cv, { height: 360, fallbackW: 400 });
    if (!ptsA.length || !ptsB.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin network_metrics', 20, 40);
      return;
    }
    const mapB = {};
    ptsB.forEach(r => { mapB[r.channel] = Number(r.strength); });
    const pts = ptsA.map(r => ({
      c: r.channel,
      x: Number(r.strength),
      y: mapB[r.channel],
      d: mapB[r.channel] - Number(r.strength)
    })).filter(p => p.y != null && !isNaN(p.y));
    if (!pts.length) return;
    const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
    const lo = Math.min(...xs, ...ys), hi = Math.max(...xs, ...ys);
    const pad = 48;
    const X = v => pad + (v - lo) / Math.max(hi - lo, 1e-9) * (w - 2 * pad);
    const Y = v => h - pad - (v - lo) / Math.max(hi - lo, 1e-9) * (h - 2 * pad);
    ctx.strokeStyle = '#e2e8f0';
    ctx.beginPath();
    ctx.moveTo(pad, h - pad);
    ctx.lineTo(w - pad, h - pad);
    ctx.lineTo(w - pad, pad);
    ctx.stroke();
    ctx.strokeStyle = '#94a3b8';
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(X(lo), Y(lo));
    ctx.lineTo(X(hi), Y(hi));
    ctx.stroke();
    ctx.setLineDash([]);
    const dmax = Math.max(1e-9, ...pts.map(p => Math.abs(p.d)));
    pts.forEach(p => {
      const t = p.d / dmax;
      const [r, g, b] = rdBu(t);
      const rad = 3.5 + 3.5 * Math.abs(t);
      ctx.beginPath();
      ctx.arc(X(p.x), Y(p.y), rad, 0, Math.PI * 2);
      ctx.fillStyle = `rgb(${r | 0},${g | 0},${b | 0})`;
      ctx.fill();
      ctx.strokeStyle = 'rgba(15,23,42,0.35)';
      ctx.lineWidth = 0.8;
      ctx.stroke();
    });
    ctx.fillStyle = '#64748b';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText((labelA || 'A') + ' strength →', w / 2, h - 12);
    ctx.save();
    ctx.translate(14, h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText((labelB || 'B') + ' strength', 0, 0);
    ctx.restore();
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText('Descriptivo (canales, no inferencia por nodo) · color = Δ strength', pad, 18);

    cv._str = { pts, X, Y, pad, w, h };
    if (!cv._strBound) {
      cv._strBound = true;
      cv.addEventListener('mousemove', ev => {
        const info = cv._str;
        if (!info) return;
        const rect = cv.getBoundingClientRect();
        const mx = ev.clientX - rect.left, my = ev.clientY - rect.top;
        let best = null, bd = 1e9;
        info.pts.forEach(p => {
          const dx = info.X(p.x) - mx, dy = info.Y(p.y) - my;
          const d2 = dx * dx + dy * dy;
          if (d2 < bd) { bd = d2; best = p; }
        });
        if (best && bd < 100) {
          showTip(ev.clientX, ev.clientY,
            `${best.c}\n${labelA}=${best.x.toFixed(3)}  ${labelB}=${best.y.toFixed(3)}\nΔ=${best.d.toFixed(3)}`);
        } else hideTip();
      });
      cv.addEventListener('mouseleave', hideTip);
    }
  }

  /**
   * Volcano plot. opts.yKey: 'p' (default) or 'q'
   * opts.xKey: 'd' (Cohen d, default) or 'rrb' (rank-biserial)
   */
  function drawVolcano(cv, edges, opts) {
    opts = opts || {};
    const yKey = opts.yKey || 'p';
    const xKey = opts.xKey || 'd';
    const { ctx, w, h } = resizeCanvas(cv, { height: 300, fallbackW: 480 });
    const rows = edges || [];
    if (!rows.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin estadísticas de edges', 20, 40);
      return;
    }
    const xs = rows.map(r => {
      if (xKey === 'rrb') {
        const v = Number(r.effect_rrb);
        return isNaN(v) ? effectValue(r) : v;
      }
      return effectValue(r);
    });
    const rawY = rows.map(r => {
      const v = yKey === 'q' ? Number(r.q_value) : Number(r.p_value);
      return -Math.log10(Math.max(v || 1e-12, 1e-12));
    });
    const pad = { l: 52, r: 20, t: 36, b: 42 };
    const absMax = Math.max(...xs.map(Math.abs), 0.2);
    const xmin = -absMax, xmax = absMax; // symmetric for negative effects
    const ymin = 0, ymax = Math.max(...rawY, 1.5) * 1.12;
    const X = v => pad.l + (v - xmin) / Math.max(xmax - xmin, 1e-9) * (w - pad.l - pad.r);
    const Y = v => pad.t + (1 - (v - ymin) / Math.max(ymax - ymin, 1e-9)) * (h - pad.t - pad.b);

    // grid
    ctx.strokeStyle = '#f1f5f9';
    ctx.fillStyle = '#94a3b8';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'right';
    for (let i = 0; i <= 4; i++) {
      const yv = ymin + (ymax - ymin) * i / 4;
      const yy = Y(yv);
      ctx.beginPath();
      ctx.moveTo(pad.l, yy);
      ctx.lineTo(w - pad.r, yy);
      ctx.stroke();
      ctx.fillText(yv.toFixed(1), pad.l - 4, yy + 3);
    }
    ctx.textAlign = 'center';
    [-absMax, 0, absMax].forEach(xv => {
      const xx = X(xv);
      ctx.beginPath();
      ctx.moveTo(xx, pad.t);
      ctx.lineTo(xx, h - pad.b);
      ctx.stroke();
      ctx.fillText(xv.toFixed(2), xx, h - pad.b + 14);
    });

    // axes
    ctx.strokeStyle = '#cbd5e1';
    ctx.beginPath();
    ctx.moveTo(pad.l, h - pad.b);
    ctx.lineTo(w - pad.r, h - pad.b);
    ctx.lineTo(w - pad.r, pad.t);
    ctx.stroke();

    // zero line
    ctx.strokeStyle = '#94a3b8';
    ctx.setLineDash([2, 2]);
    ctx.beginPath();
    ctx.moveTo(X(0), pad.t);
    ctx.lineTo(X(0), h - pad.b);
    ctx.stroke();

    // threshold line at 0.05
    const y05 = Y(-Math.log10(0.05));
    ctx.strokeStyle = '#b45309';
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(pad.l, y05);
    ctx.lineTo(w - pad.r, y05);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#b45309';
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(yKey === 'q' ? 'q=0.05' : 'p=0.05', pad.l + 4, y05 - 4);

    rows.forEach((r, i) => {
      const sig = isFdrResult(r);
      const nom = isNominalResult(r);
      ctx.beginPath();
      ctx.arc(X(xs[i]), Y(rawY[i]), sig ? 4.5 : 2.5, 0, Math.PI * 2);
      if (sig) ctx.fillStyle = '#b91c1c';
      else if (nom) ctx.fillStyle = 'rgba(194,65,12,0.65)';
      else ctx.fillStyle = 'rgba(100,116,139,0.4)';
      ctx.fill();
    });

    ctx.fillStyle = '#64748b';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText((opts.xLabel || (xKey === 'rrb' ? 'r_rb (rank-biserial)' : 'Cohen d')) + ' →', w / 2, h - 8);
    ctx.save();
    ctx.translate(14, h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(yKey === 'q' ? '−log10(q)' : '−log10(p)', 0, 0);
    ctx.restore();
    // legend
    const lx = w - pad.r - 118, ly0 = pad.t + 4;
    ctx.font = '9px sans-serif';
    ctx.textAlign = 'left';
    [['#b91c1c', 'q<0.05'], ['rgba(194,65,12,0.85)', 'p<0.05 · q≥0.05'],
     ['rgba(100,116,139,0.7)', 'p≥0.05']].forEach(function (row, i) {
      ctx.fillStyle = row[0];
      ctx.beginPath();
      ctx.arc(lx, ly0 + i * 14, 3.5, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillStyle = '#475569';
      ctx.fillText(row[1], lx + 8, ly0 + i * 14 + 3);
    });
  }

  function buildFdrSet(edgeStats) {
    const s = new Set();
    (edgeStats || []).forEach(e => {
      if (isFdrResult(e)) {
        const a = e.ch_a, b = e.ch_b;
        s.add(a < b ? a + '|' + b : b + '|' + a);
      }
    });
    return s;
  }

  function edgeLookupFactory(edgeStats) {
    const map = new Map();
    (edgeStats || []).forEach(e => {
      map.set(e.ch_a + '|' + e.ch_b, e);
      map.set(e.ch_b + '|' + e.ch_a, e);
    });
    return (a, b) => map.get(a + '|' + b) || null;
  }

  function fillTable(tbl, cols, rows, fmt) {
    const thead = tbl.querySelector('thead');
    const tbody = tbl.querySelector('tbody');
    thead.innerHTML = '<tr>' + cols.map(c => `<th>${c}</th>`).join('') + '</tr>';
    tbody.innerHTML = (rows || []).map(r => '<tr>' + cols.map(c => {
      let v = r[c];
      if (fmt && fmt[c]) v = fmt[c](v);
      else if (v == null || v === '' || v === 'missing') v = '—';
      return `<td>${v ?? '—'}</td>`;
    }).join('') + '</tr>').join('') ||
      '<tr><td colspan="' + cols.length + '" class="muted">Sin filas</td></tr>';
  }

  global.NMV = {
    BAND_ORDER, sortBandsPhysio, resizeCanvas, viridis, rdBu,
    drawHeatmap, bindHeatmap, barChart, drawGraph, drawTopo,
    drawStrength, drawVolcano, buildFdrSet, edgeLookupFactory,
    fillTable, showTip, hideTip
  };
})(window);
