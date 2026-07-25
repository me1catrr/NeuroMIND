/* NeuroMIND — shared canvas helpers for longitudinal / transversal viewers */
(function (global) {
  'use strict';

  const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"];

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
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
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
    cv._hm = { ox, oy, cell, n, ch, flat, diverging, lim };
  }

  function bindHeatmap(cv, edgeLookupFn) {
    cv.addEventListener('mousemove', ev => {
      const info = cv._hm;
      if (!info) return;
      const rect = cv.getBoundingClientRect();
      const scaleX = info.n ? (cv.clientWidth / (cv.width / (window.devicePixelRatio || 1))) : 1;
      const x = ev.clientX - rect.left;
      const y = ev.clientY - rect.top;
      const j = Math.floor((x - info.ox) / info.cell);
      const i = Math.floor((y - info.oy) / info.cell);
      if (i < 0 || j < 0 || i >= info.n || j >= info.n) {
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
        t += `\np=${Number(ed.p_value).toFixed(4)}  q=${Number(ed.q_value).toFixed(4)}  d=${Number(ed.effect_d).toFixed(3)}`;
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

  function barChart(cv, labels, values, colorFn, opts) {
    opts = opts || {};
    const { ctx, w, h } = resizeCanvas(cv, { height: opts.height || Number(cv.getAttribute('data-h')) || 220, fallbackW: 640 });
    if (!labels.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin datos', 20, 40);
      return;
    }
    const pad = { l: 44, r: 16, t: 16, b: 52 };
    const vmax = Math.max(...values.map(Math.abs), 1e-9);
    const bw = (w - pad.l - pad.r) / labels.length;
    const signed = values.some(v => v < 0) || opts.signed;
    const y0 = pad.t + (h - pad.t - pad.b) / 2;
    const base = signed ? y0 : (h - pad.b);
    const scale = signed ? (h - pad.t - pad.b) / 2 / vmax : (h - pad.t - pad.b) / vmax;
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
      ctx.fillStyle = colorFn(v, i);
      if (signed) ctx.fillRect(x, v >= 0 ? y0 - bh : y0, bw * 0.7, Math.abs(bh));
      else ctx.fillRect(x, base - Math.abs(bh), bw * 0.7, Math.abs(bh));
      ctx.save();
      ctx.translate(pad.l + i * bw + bw / 2, h - 8);
      ctx.rotate(-0.55);
      ctx.fillStyle = '#64748b';
      ctx.font = '10px sans-serif';
      ctx.textAlign = 'right';
      ctx.fillText(lab, 0, 0);
      ctx.restore();
    });
  }

  function drawGraph(cv, S, edges) {
    const { ctx, w, h } = resizeCanvas(cv, { square: true, maxSquare: 440, fallbackW: 420 });
    const cx = w / 2, cy = h / 2, R = Math.min(w, h) * 0.38;
    ctx.strokeStyle = '#94a3b8';
    ctx.beginPath();
    ctx.arc(cx, cy, R * 1.05, 0, Math.PI * 2);
    ctx.stroke();
    const pos = S.positions || {};
    const xy = {};
    (S.channels || []).forEach(c => {
      const p = pos[c.toUpperCase()] || pos[c] || [0, 0];
      xy[c] = { x: cx + p[0] * R, y: cy - p[1] * R };
    });
    const list = edges || [];
    const dmax = Math.max(1e-9, ...list.map(e => Math.abs(Number(e.diff) || 0)));
    list.forEach(e => {
      const A = xy[e.ch_a], B = xy[e.ch_b];
      if (!A || !B) return;
      const d = Number(e.diff) || 0;
      ctx.strokeStyle = d >= 0 ? 'rgba(185,28,28,0.75)' : 'rgba(29,78,216,0.75)';
      ctx.lineWidth = 1 + 2.5 * Math.abs(d) / dmax;
      ctx.beginPath();
      ctx.moveTo(A.x, A.y);
      ctx.lineTo(B.x, B.y);
      ctx.stroke();
    });
    (S.channels || []).forEach(c => {
      const p = xy[c];
      if (!p) return;
      ctx.beginPath();
      ctx.arc(p.x, p.y, 5, 0, Math.PI * 2);
      ctx.fillStyle = '#fff';
      ctx.fill();
      ctx.strokeStyle = '#0f172a';
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
      ctx.fillText('Cambia a Top-N |d| para exploración', cx, cy + 12);
    }
  }

  function drawTopo(cv, rows, positions) {
    const { ctx, w, h } = resizeCanvas(cv, { square: true, maxSquare: 400, fallbackW: 400 });
    if (!rows || !rows.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin datos espectrales', 20, 40);
      return;
    }
    const cx = w / 2, cy = h / 2, R = Math.min(w, h) * 0.38;
    ctx.strokeStyle = '#94a3b8';
    ctx.beginPath();
    ctx.arc(cx, cy, R * 1.05, 0, Math.PI * 2);
    ctx.stroke();
    const pos = positions || {};
    const lim = Math.max(1e-9, ...rows.map(r => Math.abs(Number(r.diff) || 0)));
    rows.forEach(r => {
      const c = String(r.channel);
      const p = pos[c.toUpperCase()] || pos[c] || [0, 0];
      const x = cx + p[0] * R, y = cy - p[1] * R;
      const [rr, gg, bb] = rdBu((Number(r.diff) || 0) / lim);
      const sig = Number(r.q_value) < 0.05;
      ctx.beginPath();
      ctx.arc(x, y, sig ? 9 : 6, 0, Math.PI * 2);
      ctx.fillStyle = `rgb(${rr | 0},${gg | 0},${bb | 0})`;
      ctx.fill();
      ctx.strokeStyle = sig ? '#0f172a' : '#64748b';
      ctx.lineWidth = sig ? 2 : 1;
      ctx.stroke();
      ctx.fillStyle = '#0f172a';
      ctx.font = '9px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText(c, x, y - 10);
    });
    drawColorbar(ctx, w - 40, 20, 12, h - 50, lim, true);
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
    const pts = ptsA.map(r => ({ c: r.channel, x: Number(r.strength), y: mapB[r.channel] }))
      .filter(p => p.y != null && !isNaN(p.y));
    if (!pts.length) return;
    const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
    const lo = Math.min(...xs, ...ys), hi = Math.max(...xs, ...ys);
    const pad = 48;
    const X = v => pad + (v - lo) / Math.max(hi - lo, 1e-9) * (w - 2 * pad);
    const Y = v => h - pad - (v - lo) / Math.max(hi - lo, 1e-9) * (h - 2 * pad);
    ctx.strokeStyle = '#cbd5e1';
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
    pts.forEach(p => {
      ctx.beginPath();
      ctx.arc(X(p.x), Y(p.y), 4, 0, Math.PI * 2);
      ctx.fillStyle = '#0f766e';
      ctx.fill();
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
  }

  function drawVolcano(cv, edges) {
    const { ctx, w, h } = resizeCanvas(cv, { height: 300, fallbackW: 480 });
    const rows = edges || [];
    if (!rows.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin estadísticas de edges', 20, 40);
      return;
    }
    const xs = rows.map(r => Number(r.effect_d) || 0);
    const ys = rows.map(r => -Math.log10(Math.max(Number(r.p_value) || 1e-12, 1e-12)));
    const pad = { l: 48, r: 20, t: 20, b: 40 };
    const xmin = Math.min(...xs), xmax = Math.max(...xs);
    const ymin = 0, ymax = Math.max(...ys, 1);
    const X = v => pad.l + (v - xmin) / Math.max(xmax - xmin, 1e-9) * (w - pad.l - pad.r);
    const Y = v => pad.t + (1 - (v - ymin) / Math.max(ymax - ymin, 1e-9)) * (h - pad.t - pad.b);
    ctx.strokeStyle = '#e2e8f0';
    ctx.beginPath();
    ctx.moveTo(pad.l, h - pad.b);
    ctx.lineTo(w - pad.r, h - pad.b);
    ctx.lineTo(w - pad.r, pad.t);
    ctx.stroke();
    // p=0.05 line
    const y05 = Y(-Math.log10(0.05));
    ctx.strokeStyle = '#94a3b8';
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(pad.l, y05);
    ctx.lineTo(w - pad.r, y05);
    ctx.stroke();
    ctx.setLineDash([]);
    rows.forEach((r, i) => {
      const sig = Number(r.q_value) < 0.05;
      ctx.beginPath();
      ctx.arc(X(xs[i]), Y(ys[i]), sig ? 4.5 : 2.5, 0, Math.PI * 2);
      ctx.fillStyle = sig ? '#b91c1c' : 'rgba(100,116,139,0.45)';
      ctx.fill();
    });
    ctx.fillStyle = '#64748b';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('Cohen d →', w / 2, h - 10);
    ctx.save();
    ctx.translate(14, h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('−log10(p)', 0, 0);
    ctx.restore();
  }

  function meanSem(arr) {
    if (!arr.length) return { m: NaN, sem: 0 };
    const m = arr.reduce((s, v) => s + v, 0) / arr.length;
    if (arr.length < 2) return { m, sem: 0 };
    const v = arr.reduce((s, x) => s + (x - m) * (x - m), 0) / (arr.length - 1);
    return { m, sem: Math.sqrt(v) / Math.sqrt(arr.length) };
  }

  function buildFdrSet(edgeStats) {
    const s = new Set();
    (edgeStats || []).forEach(e => {
      if (Number(e.q_value) < 0.05) {
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
      return `<td>${v ?? ''}</td>`;
    }).join('') + '</tr>').join('') ||
      '<tr><td colspan="' + cols.length + '" class="muted">Sin filas</td></tr>';
  }

  global.NMV = {
    BAND_ORDER, sortBandsPhysio, resizeCanvas, viridis, rdBu,
    drawHeatmap, bindHeatmap, barChart, drawGraph, drawTopo,
    drawStrength, drawVolcano, meanSem, buildFdrSet, edgeLookupFactory,
    fillTable, showTip, hideTip
  };
})(window);
