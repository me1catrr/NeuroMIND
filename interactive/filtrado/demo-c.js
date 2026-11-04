/* Simulación C — `filt` (causal) frente a `filtfilt` (fase cero) con un Butterworth del tipo y orden del pipeline.
 * Datos 100 % sintéticos generados en el navegador (fórmulas de interactive/filtrado/filtfilt.js, validadas contra DSP.jl en
 * tests/verificar.mjs, sección N). Script clásico, sin red, sin almacenamiento y sin dependencias; no lee ni escribe nada del proyecto. */
(function () {
  'use strict';
  var T = window.NMFiltfilt, raiz = document.getElementById('demo-c');
  if (!T || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };

  var FS = T.FS_SENAL, NREG = 4000, NYQ = FS / 2, TOT = NREG / FS, T0 = T.T_EVENTO, FMIN = 0.01;
  var ZONAS = {
    todo: [0, TOT, 'todo el registro (0–8 s)'], inicio: [0, 1.5, 'inicio (0–1,5 s)'], final: [TOT - 1.5, TOT, 'final (6,5–8 s)'],
    evento: [T0 - 0.5, T0 + 1.5, 'alrededor del evento (3,5–5,5 s)'], zoom: [T0 - 0.1, T0 + 0.3, 'zoom del evento (3,9–4,3 s)'],
    micro: [T0 - 0.03, T0 + 0.07, 'micro-zoom del evento (3,97–4,07 s)']
  };
  var ZONA_DEF = { dos: 'todo', borde: 'inicio', escalon: 'evento', impulso: 'zoom' };
  var ESC_NOMBRE = {
    dos: 'dos componentes (rápida + lenta)', borde: 'offset al inicio + dos componentes',
    escalon: 'escalón de ' + T.A_ESCALON + ' µV en t = 4 s', impulso: 'impulso de ' + T.A_IMPULSO + ' µV en t = 4 s'
  };
  var DEF = { esc: 'dos', tipo: 'hp', fc: 0.5, N: 4, f1: 10, f2: 1, off: 50, zona: 'todo', ver: { 'in': 1, filt: 1, ff: 1 } };
  var DEF_TIPO = { hp: 0.5, lp: 150 };
  var LIM = { fc: [0.05, 0.98 * NYQ], f1: [2, 100], f2: [0.1, 5], off: [0, 100], N: [1, 8] };
  var st;

  function copiar() { return JSON.parse(JSON.stringify(DEF)); }
  function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }
  function lmap(v, lo, hi) { return lo * Math.pow(hi / lo, v / 1000); }
  function linv(x, lo, hi) { return Math.round(1000 * Math.log(x / lo) / Math.log(hi / lo)); }
  function num(x, d) {                        // coma decimal y signo menos tipográfico
    var s = Math.abs(x) < 0.5 * Math.pow(10, -d) ? (0).toFixed(d) : x.toFixed(d);
    return s.replace('-', '−').replace('.', ',');
  }
  function sig(x) { return num(x, x >= 100 ? 0 : x >= 10 ? 1 : x >= 1 ? 2 : x >= 0.1 ? 2 : 3); }
  function con(x, d) { return (x > 0 ? '+' : '') + num(x, d); }

  function clampEstado() {
    st.N = clamp(Math.round(st.N), LIM.N[0], LIM.N[1]);
    ['fc', 'f1', 'f2', 'off'].forEach(function (k) { st[k] = clamp(st[k], LIM[k][0], LIM[k][1]); });
  }
  function cfg() { return { tipo: st.tipo, fs: FS, N: st.N, fc: st.fc }; }

  // ───── controles ─────
  function escribirControles() {
    $('c-esc').value = st.esc; $('c-tipo').value = st.tipo; $('c-zona').value = st.zona;
    $('c-fc').value = linv(st.fc, LIM.fc[0], LIM.fc[1]); $('c-fc-val').textContent = sig(st.fc) + ' Hz';
    $('c-n').value = st.N; $('c-n-val').textContent = String(st.N);
    $('c-f1').value = linv(st.f1, LIM.f1[0], LIM.f1[1]); $('c-f1-val').textContent = sig(st.f1) + ' Hz';
    $('c-f2').value = linv(st.f2, LIM.f2[0], LIM.f2[1]); $('c-f2-val').textContent = sig(st.f2) + ' Hz';
    $('c-off').value = st.off; $('c-off-val').textContent = num(st.off, 0) + ' µV';
    var senoides = st.esc === 'dos' || st.esc === 'borde';
    $('c-f1-etq').hidden = !senoides; $('c-f2-etq').hidden = !senoides; $('c-off-etq').hidden = st.esc !== 'borde';
    $('c-ver-in').checked = !!st.ver['in']; $('c-ver-filt').checked = !!st.ver.filt; $('c-ver-ff').checked = !!st.ver.ff;
    $('c-lim').textContent = 'Límites (Fs fija en ' + FS + ' Hz, Nyquist ' + NYQ + ' Hz): frecuencia de corte de ' + num(LIM.fc[0], 2) + ' a ' + num(LIM.fc[1], 1) +
      ' Hz; orden de diseño entero de 1 a 8; componente rápida de ' + LIM.f1[0] + ' a ' + LIM.f1[1] + ' Hz; componente lenta de ' + num(LIM.f2[0], 1) + ' a ' + LIM.f2[1] +
      ' Hz; offset de 0 a ' + LIM.off[1] + ' µV. Los controles no admiten valores fuera de rango.';
  }

  // ───── cálculo ─────
  function calcular() {
    var d = T.disenar(cfg()), x = T.senal(st.esc, { f1: st.f1, f2: st.f2, offset: st.off }, NREG, FS);
    return { d: d, x: x, yc: T.filt(d, x), yz: T.filtfilt(d, x) };
  }

  // ───── dibujo ─────
  var W = 760, L = 58, R = 18;
  function medir() {
    var w = $('c-tiempo').getBoundingClientRect().width;
    W = Math.max(300, Math.round(w || 760)); L = W < 480 ? 54 : 58; R = W < 480 ? 12 : 18;
  }
  function trazo(xs, ys) { var d = ''; for (var i = 0; i < xs.length; i++) d += (i ? 'L' : 'M') + xs[i].toFixed(1) + ',' + ys[i].toFixed(1) + ' '; return d; }
  function pasoNice(rango, max) {
    var p = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
    for (var i = 0; i < p.length; i++) if (rango / p[i] <= (max || 6)) return p[i];
    return 2000;
  }
  function pasoT(span) {
    var p = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2];
    for (var i = 0; i < p.length; i++) if (span / p[i] <= 8) return p[i];
    return 2;
  }

  function dibujarTiempo(S) {
    var Z = ZONAS[st.zona], t0 = Z[0], t1 = Z[1], H = 300, B = 36, Tp = 12, h = H - B - Tp;
    var i0 = Math.round(t0 * FS), i1 = Math.min(NREG, Math.round(t1 * FS)), idx = [], i, j;
    for (i = i0; i < i1; i++) idx.push(i);
    var m = 0.5;
    [[S.x, st.ver['in']], [S.yc, st.ver.filt], [S.yz, st.ver.ff]].forEach(function (q) {
      if (!q[1]) return;
      for (j = 0; j < idx.length; j++) m = Math.max(m, Math.abs(q[0][idx[j]]));
    });
    var pas = pasoNice(2 * m * 1.08), ymax = Math.ceil(m * 1.08 / pas) * pas;
    function X(t) { return L + (t - t0) / (t1 - t0) * (W - L - R); }
    function Y(v) { return Tp + (ymax - clamp(v, -ymax, ymax)) / (2 * ymax) * h; }
    var s = '', v, k, decY = pas < 1 ? 1 : 0;
    for (k = -Math.round(ymax / pas); k <= Math.round(ymax / pas); k++) {
      v = k * pas;
      s += '<line class="rejilla" x1="' + L + '" x2="' + (W - R) + '" y1="' + Y(v) + '" y2="' + Y(v) + '"/><text x="' + (L - 6) + '" y="' + (Y(v) + 4) + '" text-anchor="end" font-size="12">' + num(v, decY) + '</text>';
    }
    var tp = pasoT(t1 - t0), decT = tp < 0.01 ? 3 : tp < 0.1 ? 2 : tp < 1 ? 1 : 0;
    for (k = Math.ceil(t0 / tp - 1e-9); k * tp <= t1 + 1e-9; k++) {
      v = k * tp;
      s += '<line class="rejilla" x1="' + X(v) + '" x2="' + X(v) + '" y1="' + Tp + '" y2="' + (H - B) + '"/><text x="' + X(v) + '" y="' + (H - B + 16) + '" text-anchor="middle" font-size="12">' + num(v, decT) + '</text>';
    }
    s += '<text x="' + ((L + W - R) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" font-size="12">' + (W < 480 ? 'Tiempo (s)' : 'Tiempo (s) — zona: ' + Z[2]) + '</text>';
    s += '<text transform="translate(14 ' + (Tp + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Amplitud (µV)</text>';
    s += '<rect x="' + L + '" y="' + Tp + '" width="' + (W - L - R) + '" height="' + h + '" fill="none" class="eje"/>';
    if ((st.esc === 'escalon' || st.esc === 'impulso') && T0 >= t0 && T0 <= t1) {
      s += '<line class="guia" x1="' + X(T0) + '" x2="' + X(T0) + '" y1="' + Tp + '" y2="' + (H - B) + '"/><text x="' + (X(T0) + 4) + '" y="' + (Tp + 12) + '" font-size="12">' + (st.esc === 'escalon' ? 'escalón' : 'impulso') + '</text>';
    }
    var xs = idx.map(function (q) { return X(q / FS); });
    if (st.ver['in']) s += '<path class="curva-0" d="' + trazo(xs, idx.map(function (q) { return Y(S.x[q]); })) + '"/>';
    if (st.ver.filt) s += '<path class="curva-1" d="' + trazo(xs, idx.map(function (q) { return Y(S.yc[q]); })) + '"/>';
    if (st.ver.ff) s += '<path class="curva-2" d="' + trazo(xs, idx.map(function (q) { return Y(S.yz[q]); })) + '"/>';
    $('c-tiempo').setAttribute('viewBox', '0 0 ' + W + ' ' + H); $('c-tiempo').innerHTML = s;
  }

  function dibujarRetardo(S) {
    var H = 250, B = 34, Tp = 12, h = H - B - Tp, d = S.d, i, k;
    function xpix(f) { return L + (Math.log10(f) - Math.log10(FMIN)) / (Math.log10(NYQ) - Math.log10(FMIN)) * (W - L - R); }
    var hz = [], g = [];
    for (i = 0; i < 500; i++) { var f = FMIN * Math.pow(0.9995 * NYQ / FMIN, i / 499); hz.push(f); g.push(T.retardoGrupo(d, f) * 1000); }
    var gmax = 1e-9, gmin = 0;
    g.forEach(function (q) { gmax = Math.max(gmax, q); gmin = Math.min(gmin, q); });
    var seg = gmax >= 3000, fu = seg ? 1000 : 1, un = seg ? 's' : 'ms';
    var pas = pasoNice((gmax - gmin) / fu * 1.05, 5), top = Math.ceil(gmax / fu * 1.05 / pas) * pas, bot = gmin < 0 ? Math.floor(gmin / fu / pas) * pas : 0;
    function Y(v) { return Tp + (top - clamp(v / fu, bot, top)) / (top - bot) * h; }
    var s = '', dec = pas < 1 ? 1 : 0, v;
    for (k = Math.round(bot / pas); k * pas <= top + 1e-9; k++) {
      v = k * pas;
      s += '<line class="rejilla" x1="' + L + '" x2="' + (W - R) + '" y1="' + Y(v * fu) + '" y2="' + Y(v * fu) + '"/><text x="' + (L - 6) + '" y="' + (Y(v * fu) + 4) + '" text-anchor="end" font-size="12">' + num(v, dec) + '</text>';
    }
    s += '<text transform="translate(14 ' + (Tp + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Retardo de grupo (' + un + ')</text>';
    for (k = -2; k <= 2; k++) {
      f = Math.pow(10, k);
      if (f < FMIN || f > NYQ) continue;
      s += '<line class="rejilla" x1="' + xpix(f) + '" x2="' + xpix(f) + '" y1="' + Tp + '" y2="' + (H - B) + '"/>';
      if (xpix(NYQ) - xpix(f) > 90) s += '<text x="' + xpix(f) + '" y="' + (H - B + 16) + '" text-anchor="middle" font-size="12">' + num(f, k < 0 ? -k : 0) + '</text>';
    }
    s += '<line class="rejilla" x1="' + xpix(NYQ) + '" x2="' + xpix(NYQ) + '" y1="' + Tp + '" y2="' + (H - B) + '"/>';
    s += '<text x="' + xpix(NYQ) + '" y="' + (H - B + 16) + '" text-anchor="end" font-size="12">Nyquist ' + NYQ + '</text>';
    s += '<text x="' + ((L + W - R) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" font-size="12">Frecuencia (Hz, escala logarítmica)</text>';
    s += '<rect x="' + L + '" y="' + Tp + '" width="' + (W - L - R) + '" height="' + h + '" fill="none" class="eje"/>';
    s += '<path class="curva-1" d="' + trazo(hz.map(xpix), g.map(Y)) + '"/>';
    s += '<line class="curva-2" x1="' + L + '" x2="' + (W - R) + '" y1="' + Y(0) + '" y2="' + Y(0) + '"/>';
    function marca(f, txt, dy, col) {
      if (f <= FMIN || f >= NYQ) return '';
      return '<line x1="' + xpix(f) + '" x2="' + xpix(f) + '" y1="' + Tp + '" y2="' + (H - B) + '" style="stroke:' + col + ';stroke-width:1.5;stroke-dasharray:1 4"/>' +
        '<text x="' + (xpix(f) + 4) + '" y="' + (Tp + dy) + '" font-size="12" style="fill:' + (col === 'var(--nm-plum)' ? 'var(--nm-t-plum)' : 'var(--nm-text)') + '">' + txt + '</text>' +
        '<circle cx="' + xpix(f) + '" cy="' + Y(T.retardoGrupo(d, f) * 1000) + '" r="4.5" style="fill:var(--nm-orange)"/>';
    }
    s += marca(st.fc, 'fc', 12, 'var(--nm-muted)');
    if (st.esc === 'dos' || st.esc === 'borde') s += marca(st.f1, 'f₁', 26, 'var(--nm-plum)') + marca(st.f2, 'f₂', 40, 'var(--nm-plum)');
    $('c-retardo').setAttribute('viewBox', '0 0 ' + W + ' ' + H); $('c-retardo').innerHTML = s;
  }

  // ───── lectura ─────
  function ms(s) { return num(s * 1000, Math.abs(s * 1000) >= 10 ? 1 : 2); }
  function lecturaComponente(d, hz, amp, etq) {
    var h = T.respuestaSOS(d, hz), m = Math.hypot(h[0], h[1]), ph = Math.atan2(h[1], h[0]) * 180 / Math.PI, tf = T.retardoFase(d, hz), tg = T.retardoGrupo(d, hz);
    var dir = Math.abs(tf) < 5e-6 ? 'sin desplazamiento apreciable' : tf > 0 ? 'la onda de salida llega ' + ms(tf) + ' ms después' : 'la onda de salida se adelanta ' + ms(-tf) + ' ms';
    var db1 = 20 * Math.log10(Math.max(m, 1e-12)), db2 = 2 * db1;
    return '<p><strong>' + etq + ' de ' + sig(hz) + ' Hz</strong> (entrada de ' + amp + ' µV): <em>una pasada (<code>filt</code>)</em> ganancia ' + num(db1, 2) + ' dB (amplitud de salida ' + num(amp * m, 2) + ' µV), fase ' + con(ph, 1) + '°, ' + dir +
      ', retardo de grupo ' + ms(tg) + ' ms. <em>Doble pasada (<code>filtfilt</code>)</em> ganancia ' + num(db2, 2) + ' dB (amplitud ' + num(amp * m * m, 2) + ' µV), fase 0°, desplazamiento 0 ms, retardo de grupo 0 ms.</p>';
  }
  function lectura(S) {
    var p = [], d = S.d, nb = d.nb, rell = T.relleno(d, NREG), n0 = Math.round(T0 * FS), yc = S.yc, yz = S.yz;
    var pip = (st.tipo === 'hp' && st.fc === 0.5 && st.N === 4) ? ' (los del paso alto del pipeline)' : (st.tipo === 'lp' && st.fc === 150 && st.N === 4) ? ' (los del paso bajo del pipeline)' : '';
    p.push('<p><strong>Filtro:</strong> ' + (st.tipo === 'hp' ? 'paso alto' : 'paso bajo') + ' de ' + sig(st.fc) + ' Hz, orden de diseño ' + st.N + pip + '; ' + nb + ' biquads por pasada; relleno de bordes de <code>filtfilt</code>: ' + rell + ' muestras (' + num(1000 * rell / FS, 0) + ' ms). <strong>Escenario:</strong> ' + ESC_NOMBRE[st.esc] + '. <strong>Zona:</strong> ' + ZONAS[st.zona][2] + '.</p>');
    if (st.esc === 'dos' || st.esc === 'borde') {
      p.push(lecturaComponente(d, st.f2, T.A2, 'Componente lenta') + lecturaComponente(d, st.f1, T.A1, 'Componente rápida'));
      p.push('<p>Una pasada no retrasa todas las frecuencias igual (el retardo de grupo cambia con la frecuencia); doble pasada las deja a todas sin desplazamiento, pero con la ganancia elevada al cuadrado.</p>');
    }
    if (st.esc === 'borde') {
      var ic = T.estacionaria(d, { f1: st.f1, f2: st.f2, offset: st.off }, NREG, FS, 'filt', true), iz = T.estacionaria(d, { f1: st.f1, f2: st.f2, offset: st.off }, NREG, FS, 'filtfilt', true);
      function dv(y, ide, a, b) { return num(T.maxDif(y, ide, Math.round(a * FS), Math.round(b * FS)), 2); }
      p.push('<p><strong>Efecto de borde.</strong> La entrada lleva un offset de ' + num(st.off, 0) + ' µV desde la primera muestra (con paso alto, el filtro lo ve como un escalón en t = 0). Desviación máxima de la salida respecto al régimen permanente ideal (offset, ganancia y fase de cada componente), en µV: ' +
        'primer 0,5 s, <code>filt</code> ' + dv(yc, ic, 0, 0.5) + ' y <code>filtfilt</code> ' + dv(yz, iz, 0, 0.5) + '; entre 3 y 4 s, <code>filt</code> ' + dv(yc, ic, 3, 4) + ' y <code>filtfilt</code> ' + dv(yz, iz, 3, 4) +
        '; último 0,5 s, <code>filt</code> ' + dv(yc, ic, TOT - 0.5, TOT) + ' y <code>filtfilt</code> ' + dv(yz, iz, TOT - 0.5, TOT) + '.</p>');
      p.push('<p>Primera muestra de la salida: <code>filt</code> ' + num(yc[0], 2) + ' µV; <code>filtfilt</code> ' + num(yz[0], 2) + ' µV (entrada: ' + num(S.x[0], 2) + ' µV).' + (st.tipo === 'lp' ? ' Con paso bajo el offset es una componente de 0 Hz que pasa intacta: no hay transitorio de offset.' : '') + '</p>');
    }
    if (st.esc === 'escalon') {
      var fin = st.tipo === 'lp' ? T.A_ESCALON : 0, tol = 0.05 * T.A_ESCALON;
      function asent(y) { var t = T.asentamiento(y, n0, fin, tol, FS); return t === null ? 'no se asienta en los ' + num(TOT - T0, 0) + ' s restantes' : num(t, 3) + ' s'; }
      function sobre(y) { return Math.max(0, T.maxAbs(y, n0, NREG) - (st.tipo === 'lp' ? T.A_ESCALON : 0)); }
      var pre = T.maxAbs(yz, n0 - 100, n0), prec = T.maxAbs(yc, n0 - 100, n0);
      p.push('<p><strong>Respuesta al escalón de ' + T.A_ESCALON + ' µV.</strong> Antes del escalón (0,2 s previos): |salida| máxima de <code>filt</code> = ' + num(prec, 3) + ' µV (un filtro causal no puede anticiparse); de <code>filtfilt</code> = ' + num(pre, 3) + ' µV (la respuesta es bilateral: hay precursor).</p>');
      p.push('<p>Tiempo hasta quedar dentro de ±' + num(tol, 2) + ' µV (5 % del escalón) del valor final (' + num(fin, 0) + ' µV): <code>filt</code>, ' + asent(yc) + '; <code>filtfilt</code>, ' + asent(yz) + '. ' + (st.tipo === 'lp'
        ? 'Sobreoscilación por encima de ' + num(T.A_ESCALON, 0) + ' µV: <code>filt</code> ' + num(sobre(yc), 2) + ' µV; <code>filtfilt</code> ' + num(sobre(yz), 2) + ' µV.'
        : 'Pico de |salida| tras el escalón: <code>filt</code> ' + num(T.maxAbs(yc, n0, NREG), 2) + ' µV; <code>filtfilt</code> ' + num(T.maxAbs(yz, n0, NREG), 2) + ' µV.') + '</p>');
    }
    if (st.esc === 'impulso') {
      function pico(y) { var mx = 0, im = n0; for (var i = 0; i < NREG; i++) if (Math.abs(y[i]) > mx) { mx = Math.abs(y[i]); im = i; } return { v: mx, lag: (im - n0) / FS }; }
      var pc = pico(yc), pz = pico(yz), antes = T.maxAbs(yz, n0 - 100, n0), antesC = T.maxAbs(yc, n0 - 100, n0);
      p.push('<p><strong>Respuesta al impulso de ' + T.A_IMPULSO + ' µV.</strong> <code>filt</code>: pico de |salida| ' + num(pc.v, 2) + ' µV, ' + (pc.lag === 0 ? 'en el mismo instante del impulso' : num(pc.lag * 1000, 0) + ' ms después'
        + ' del impulso') + '; |salida| máxima en los 0,2 s previos: ' + num(antesC, 3) + ' µV (causal). <code>filtfilt</code>: pico ' + num(pz.v, 2) + ' µV, ' + (pz.lag === 0 ? 'en el mismo instante del impulso (respuesta simétrica)' : num(pz.lag * 1000, 0) + ' ms respecto al impulso') + '; |salida| máxima en los 0,2 s previos: ' + num(antes, 2) + ' µV (' + num(100 * antes / Math.max(pz.v, 1e-12), 1) + ' % de su pico): oscilación que precede al evento.</p>');
    }
    $('c-lectura').innerHTML = p.join('');
  }

  function render() {
    medir(); clampEstado(); escribirControles();
    var S = calcular();
    dibujarTiempo(S); dibujarRetardo(S); lectura(S);
  }

  // ───── eventos ─────
  function leerRango(id, fn) { $(id).addEventListener('input', function () { var v = Number($(id).value); if (isFinite(v)) { fn(v); render(); } }); }
  $('c-esc').addEventListener('change', function () { st.esc = $('c-esc').value; st.zona = ZONA_DEF[st.esc]; render(); });
  $('c-tipo').addEventListener('change', function () { st.tipo = $('c-tipo').value; st.fc = DEF_TIPO[st.tipo]; render(); });
  $('c-zona').addEventListener('change', function () { st.zona = $('c-zona').value; render(); });
  ['c-ver-in', 'c-ver-filt', 'c-ver-ff'].forEach(function (id) {
    $(id).addEventListener('change', function () { st.ver = { 'in': $('c-ver-in').checked ? 1 : 0, filt: $('c-ver-filt').checked ? 1 : 0, ff: $('c-ver-ff').checked ? 1 : 0 }; render(); });
  });
  leerRango('c-fc', function (v) { st.fc = lmap(v, LIM.fc[0], LIM.fc[1]); });
  leerRango('c-n', function (v) { st.N = v; });
  leerRango('c-f1', function (v) { st.f1 = lmap(v, LIM.f1[0], LIM.f1[1]); });
  leerRango('c-f2', function (v) { st.f2 = lmap(v, LIM.f2[0], LIM.f2[1]); });
  leerRango('c-off', function (v) { st.off = v; });
  $('c-reset').addEventListener('click', function () { st = copiar(); render(); });

  var tm = null;
  window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  st = copiar(); render();
})();
