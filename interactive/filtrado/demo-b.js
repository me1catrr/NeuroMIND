/* Simulación B — una señal sintética antes y después de filtrar (filt causal frente a filtfilt de fase cero).
 * Datos 100 % sintéticos y salidas precalculadas con DSP.jl (scripts/generar_filtrado.jl → data/synthetic/filtrado/demo-b-datos.js).
 * Se aprovecha que los filtros son lineales: la salida de la suma es la suma de las salidas por componente
 * (comprobado en scripts/generar_filtrado.jl y tests/verificar.mjs). Script clásico, sin red ni dependencias. */
(function () {
  'use strict';
  var D = window.NM_DEMO_B, F = window.NMFiltros;
  var raiz = document.getElementById('demo-b');
  if (!D || !F || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };

  var FS = D.fs, N = D.n;
  var ZONAS = {
    todo: [0, 4, 'Todo (0–4 s)'], inicio: [0, 0.6, 'Inicio (0–0,6 s)'], escalon: [1.6, 2.6, 'Alrededor del escalón (1,6–2,6 s)'],
    impulso: [2.8, 3.6, 'Alrededor del impulso (2,8–3,6 s)'], final: [3.4, 4, 'Final (3,4–4 s)']
  };
  var NOMBRES = {
    offset: 'Desplazamiento continuo (30 µV)', deriva: 'Deriva de 0,2 Hz (20 µV)', alfa: 'Alfa de 10 Hz (10 µV)',
    red: 'Red de 50 Hz (5 µV)', escalon: 'Escalón (15 µV en t = 2 s)', impulso: 'Impulso (100 µV en t = 3 s)'
  };
  var ETIQ = {
    'hp0.1': 'Paso alto 0,1 Hz', 'hp0.5': 'Paso alto 0,5 Hz (pipeline)', 'hp2': 'Paso alto 2 Hz',
    'lp20': 'Paso bajo 20 Hz', 'lp150': 'Paso bajo 150 Hz (pipeline)', 'notch50': 'Notch 49,5–50,5 Hz (pipeline)'
  };
  var DEF = { filtro: 'hp0.5', zona: 'escalon', comps: { alfa: 1, red: 1, escalon: 1 }, ver: { in: 1, filt: 1, ff: 1 } };
  var st;

  function copiar() { return JSON.parse(JSON.stringify(DEF)); }
  function num(x, d) {
    var s = Math.abs(x) < 5 * Math.pow(10, -d - 1) ? (0).toFixed(d) : x.toFixed(d);
    return s.replace('-', '−').replace('.', ',');
  }

  // ───── construcción de controles ─────
  (function () {
    var s = $('b-filtro'), h = '';
    D.filtros.forEach(function (f) { h += '<option value="' + f.id + '">' + (ETIQ[f.id] || f.id) + '</option>'; });
    s.innerHTML = h;
    var c = '';
    D.componentes.forEach(function (k) {
      c += '<label><input type="checkbox" data-comp="' + k + '" /> ' + NOMBRES[k] + '</label>';
    });
    $('b-comps').innerHTML = c;
  })();

  function leerControles() {
    st.filtro = $('b-filtro').value; st.zona = $('b-zona').value;
    st.comps = {};
    Array.prototype.forEach.call(raiz.querySelectorAll('[data-comp]'), function (e) { if (e.checked) st.comps[e.getAttribute('data-comp')] = 1; });
    st.ver = { in: $('b-ver-in').checked ? 1 : 0, filt: $('b-ver-filt').checked ? 1 : 0, ff: $('b-ver-ff').checked ? 1 : 0 };
  }
  function escribirControles() {
    $('b-filtro').value = st.filtro; $('b-zona').value = st.zona;
    Array.prototype.forEach.call(raiz.querySelectorAll('[data-comp]'), function (e) { e.checked = !!st.comps[e.getAttribute('data-comp')]; });
    $('b-ver-in').checked = !!st.ver.in; $('b-ver-filt').checked = !!st.ver.filt; $('b-ver-ff').checked = !!st.ver.ff;
  }

  // ───── cálculo ─────
  function suma(get) {
    var y = new Float64Array(N), k, i, v;
    Object.keys(st.comps).forEach(function (c) { v = get(c); for (i = 0; i < N; i++) y[i] += v[i]; });
    return y;
  }
  function señales() {
    var f = D.filtros.filter(function (q) { return q.id === st.filtro; })[0];
    return {
      f: f,
      x: suma(function (c) { return D.entrada[c]; }),
      yc: suma(function (c) { return f.filt[c]; }),
      yz: suma(function (c) { return f.filtfilt[c]; })
    };
  }
  function maxAbs(y, a, b) {
    var m = 0, i0 = Math.max(0, Math.round(a * FS)), i1 = Math.min(N, Math.round(b * FS));
    for (var i = i0; i < i1; i++) m = Math.max(m, Math.abs(y[i]));
    return m;
  }
  function envolvente(y, a, b) {       // amplitud de una sinusoide con el mismo valor eficaz
    var s = 0, i0 = Math.round(a * FS), i1 = Math.round(b * FS);
    for (var i = i0; i < i1; i++) s += y[i] * y[i];
    return Math.sqrt(2 * s / (i1 - i0));
  }

  // ───── dibujo ─────
  var W = 760, L = 58, R = 18;
  function medir() {
    var w = $('b-tiempo').getBoundingClientRect().width;
    W = Math.max(300, Math.round(w || 760)); L = W < 480 ? 46 : 58; R = W < 480 ? 12 : 18;
  }
  function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }
  function trazo(xs, ys) { var d = ''; for (var i = 0; i < xs.length; i++) d += (i ? 'L' : 'M') + xs[i].toFixed(1) + ',' + ys[i].toFixed(1) + ' '; return d; }
  function paso(rango) { var r = rango, p = [1, 2, 5, 10, 20, 50, 100, 200, 500]; for (var i = 0; i < p.length; i++) if (r / p[i] <= 6) return p[i]; return 1000; }

  function dibujarTiempo(S) {
    var Z = ZONAS[st.zona], t0 = Z[0], t1 = Z[1], H = 300, B = 36, T = 12, h = H - B - T;
    var i0 = Math.round(t0 * FS), i1 = Math.min(N, Math.round(t1 * FS)), idx = [], i;
    for (i = i0; i < i1; i++) idx.push(i);
    var m = 1;
    [S.x, st.ver.filt ? S.yc : null, st.ver.ff ? S.yz : null].forEach(function (y, k) {
      if (!y || (k === 0 && !st.ver.in)) return;
      for (var j = 0; j < idx.length; j++) m = Math.max(m, Math.abs(y[idx[j]]));
    });
    var ymax = Math.ceil(m * 1.08), pas = paso(2 * ymax); ymax = Math.ceil(ymax / pas) * pas;
    function X(t) { return L + (t - t0) / (t1 - t0) * (W - L - R); }
    function Y(v) { return T + (ymax - clamp(v, -ymax, ymax)) / (2 * ymax) * h; }
    var s = '', v;
    for (v = -ymax; v <= ymax + 1e-9; v += pas) {
      s += '<line class="rejilla" x1="' + L + '" x2="' + (W - R) + '" y1="' + Y(v) + '" y2="' + Y(v) + '"/><text x="' + (L - 6) + '" y="' + (Y(v) + 4) + '" text-anchor="end" font-size="12">' + num(v, 0) + '</text>';
    }
    var tp = (t1 - t0) > 2 ? 0.5 : (t1 - t0) > 0.8 ? 0.2 : 0.1;
    for (v = Math.ceil(t0 / tp) * tp; v <= t1 + 1e-9; v += tp) {
      s += '<line class="rejilla" x1="' + X(v) + '" x2="' + X(v) + '" y1="' + T + '" y2="' + (H - B) + '"/><text x="' + X(v) + '" y="' + (H - B + 16) + '" text-anchor="middle" font-size="12">' + num(v, 1) + '</text>';
    }
    s += '<text x="' + ((L + W - R) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" font-size="12">' + (W < 480 ? 'Tiempo (s)' : 'Tiempo (s) — zona: ' + Z[2]) + '</text>';
    s += '<text transform="translate(14 ' + (T + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Amplitud (µV)</text>';
    s += '<rect x="' + L + '" y="' + T + '" width="' + (W - L - R) + '" height="' + h + '" fill="none" class="eje"/>';
    [[2, 'escalon', 'escalón'], [3, 'impulso', 'impulso']].forEach(function (e) {
      if (st.comps[e[1]] && e[0] >= t0 && e[0] <= t1) s += '<line class="guia" x1="' + X(e[0]) + '" x2="' + X(e[0]) + '" y1="' + T + '" y2="' + (H - B) + '"/><text x="' + (X(e[0]) + 4) + '" y="' + (T + 12) + '" font-size="12">' + e[2] + '</text>';
    });
    var xs = idx.map(function (j) { return X(j / FS); });
    if (st.ver.in) s += '<path class="curva-0" d="' + trazo(xs, idx.map(function (j) { return Y(S.x[j]); })) + '"/>';
    if (st.ver.filt) s += '<path class="curva-1" d="' + trazo(xs, idx.map(function (j) { return Y(S.yc[j]); })) + '"/>';
    if (st.ver.ff) s += '<path class="curva-2" d="' + trazo(xs, idx.map(function (j) { return Y(S.yz[j]); })) + '"/>';
    $('b-tiempo').setAttribute('viewBox', '0 0 ' + W + ' ' + H); $('b-tiempo').innerHTML = s;
  }

  function dibujarEspectro(S) {
    var H = 260, B = 36, T = 12, h = H - B - T, FMN = 0.25, FMX = FS / 2, YMN = -80, YMX = 40, NFFT = 4096;
    function X(f) { return L + (Math.log10(f) - Math.log10(FMN)) / (Math.log10(FMX) - Math.log10(FMN)) * (W - L - R); }
    function Y(db) { return T + (YMX - clamp(db, YMN, YMX)) / (YMX - YMN) * h; }
    var s = '', db, k;
    for (db = YMN; db <= YMX; db += 20) s += '<line class="rejilla" x1="' + L + '" x2="' + (W - R) + '" y1="' + Y(db) + '" y2="' + Y(db) + '"/><text x="' + (L - 6) + '" y="' + (Y(db) + 4) + '" text-anchor="end" font-size="12">' + num(db, 0) + '</text>';
    [0.25, 1, 10, 100].forEach(function (f) {
      s += '<line class="rejilla" x1="' + X(f) + '" x2="' + X(f) + '" y1="' + T + '" y2="' + (H - B) + '"/>';
      if (X(FMX) - X(f) > 80) s += '<text x="' + X(f) + '" y="' + (H - B + 16) + '" text-anchor="middle" font-size="12">' + num(f, f < 1 ? 2 : 0) + '</text>';
    });
    s += '<text x="' + X(FMX) + '" y="' + (H - B + 16) + '" text-anchor="end" font-size="12">100 · Nyquist 250</text>';
    s += '<text x="' + ((L + W - R) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" font-size="12">' + (W < 480 ? 'Frecuencia (Hz, log)' : 'Frecuencia (Hz, escala logarítmica) — espectro de los 4 s completos') + '</text>';
    s += '<text transform="translate(14 ' + (T + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Amplitud (dB re 1 µV)</text>';
    s += '<rect x="' + L + '" y="' + T + '" width="' + (W - L - R) + '" height="' + h + '" fill="none" class="eje"/>';
    function curva(y, cls) {
      var A = F.espectro(y, NFFT), xs = [], ys = [];
      for (k = 1; k <= NFFT / 2; k++) { var f = k * FS / NFFT; if (f < FMN || f > FMX) continue; xs.push(X(f)); ys.push(Y(20 * Math.log10(Math.max(A[k], 1e-12)))); }
      return '<path class="' + cls + '" d="' + trazo(xs, ys) + '"/>';
    }
    if (st.ver.in) s += curva(S.x, 'curva-0');
    if (st.ver.filt) s += curva(S.yc, 'curva-1');
    if (st.ver.ff) s += curva(S.yz, 'curva-2');
    $('b-espectro').setAttribute('viewBox', '0 0 ' + W + ' ' + H); $('b-espectro').innerHTML = s;
  }

  function lectura(S) {
    var Z = ZONAS[st.zona], p = [], ids = Object.keys(st.comps);
    if (!ids.length) { $('b-lectura').innerHTML = '<p>Activa al menos un componente de la señal.</p>'; return; }
    p.push('<p><strong>Filtro:</strong> ' + (ETIQ[S.f.id] || S.f.id) + ', orden de diseño ' + D.orden + ' (' + S.f.biquads + ' biquads por pasada; relleno de bordes de <code>filtfilt</code>: ' + S.f.relleno_muestras + ' muestras = ' + num(1000 * S.f.relleno_muestras / FS, 0) + ' ms). <strong>Zona:</strong> ' + Z[2] + '.</p>');
    p.push('<p><strong>Máximo |amplitud| en la zona:</strong> entrada ' + num(maxAbs(S.x, Z[0], Z[1]), 2) + ' µV; <code>filt</code> ' + num(maxAbs(S.yc, Z[0], Z[1]), 2) + ' µV; <code>filtfilt</code> ' + num(maxAbs(S.yz, Z[0], Z[1]), 2) + ' µV.</p>');
    // Respuesta del propio evento, aislada del resto de componentes (posible porque los filtros son lineales).
    [[2, 'escalon', 'del escalón'], [3, 'impulso', 'del impulso']].forEach(function (e) {
      if (!st.comps[e[1]] || e[0] - 0.2 < Z[0] - 1e-9 || e[0] > Z[1]) return;
      var a = e[0] - 0.2, b = e[0] - 1 / FS;
      p.push('<p><strong>Respuesta ' + e[2] + ' en los 0,2 s anteriores a t = ' + num(e[0], 1) + ' s</strong> (se aísla el propio evento del resto de componentes): |salida| máxima de <code>filt</code> = ' + num(maxAbs(S.f.filt[e[1]], a, b), 3) + ' µV; de <code>filtfilt</code> = ' + num(maxAbs(S.f.filtfilt[e[1]], a, b), 3) + ' µV. Un filtro causal no puede responder antes de que llegue la entrada; la respuesta bilateral de <code>filtfilt</code> sí puede.</p>');
    });
    if (Z[0] === 0 && Z[1] <= 1) {
      p.push('<p><strong>En los primeros 0,2 s:</strong> |salida| máxima de <code>filt</code> = ' + num(maxAbs(S.yc, 0, 0.2), 2) + ' µV; de <code>filtfilt</code> = ' + num(maxAbs(S.yz, 0, 0.2), 2) + ' µV (entrada: ' + num(maxAbs(S.x, 0, 0.2), 2) + ' µV).</p>');
    }
    if (st.comps.red && S.f.tipo === 'bs') {
      var win = function (y, a, b) { return num(envolvente(y, a, b), 2); };
      p.push('<p><strong>Amplitud de la componente de 50 Hz en la salida</strong> (entrada: 5 µV; solo si «Red» es el único componente activo la cifra es limpia): <code>filt</code> ' + win(S.yc, 0, 0.2) + ' µV (0–0,2 s), ' + win(S.yc, 1, 1.2) + ' µV (1–1,2 s), ' + win(S.yc, 3.8, 4) + ' µV (3,8–4 s); <code>filtfilt</code> ' + win(S.yz, 0, 0.2) + ' µV, ' + win(S.yz, 1, 1.2) + ' µV, ' + win(S.yz, 3.8, 4) + ' µV.</p>');
    }
    $('b-lectura').innerHTML = p.join('');
  }

  function render() {
    medir();
    var S = señales();
    dibujarTiempo(S); dibujarEspectro(S); lectura(S);
  }

  ['b-filtro', 'b-zona', 'b-ver-in', 'b-ver-filt', 'b-ver-ff'].forEach(function (id) { $(id).addEventListener('change', function () { leerControles(); render(); }); });
  Array.prototype.forEach.call(raiz.querySelectorAll('[data-comp]'), function (e) { e.addEventListener('change', function () { leerControles(); render(); }); });
  $('b-reset').addEventListener('click', function () { st = copiar(); escribirControles(); render(); });

  var tm = null;
  window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  st = copiar(); escribirControles(); render();
})();
