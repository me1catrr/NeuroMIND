/* Simulación «Del retardo de conducción al desfase». Datos sintéticos: solo geometría de sinusoides (interactive/esclerosis-multiple/retardo.js).
 * No lee nada del proyecto ni calcula wPLI. Script clásico, sin red. */
(function () {
  'use strict';
  var R = window.NMRetardo, raiz = document.getElementById('demo-retardo');
  if (!R || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };

  var V_MIEL = 55, V_DESM = 7.5;                    // valores del informe (sin cita): 50–60 y 5–10 m/s; se usan los puntos medios
  var DEF = { f: 10, L: 100, v: V_MIEL };
  var st = Object.assign({}, DEF);
  var F_MIN = 1, F_MAX = 40, L_MIN = 10, L_MAX = 300, V_MIN = 0.5, V_MAX = 80;

  function lmap(x, lo, hi) { return lo * Math.pow(hi / lo, x / 1000); }
  function linv(x, lo, hi) { return Math.round(1000 * Math.log(x / lo) / Math.log(hi / lo)); }
  function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }
  function num(x, d) {
    var s = Math.abs(x) < Math.pow(10, -d) / 2 ? (0).toFixed(d) : x.toFixed(d);
    return s.replace('-', '−').replace('.', ',');
  }

  function escribir() {
    st.f = clamp(st.f, F_MIN, F_MAX); st.L = clamp(st.L, L_MIN, L_MAX); st.v = clamp(st.v, V_MIN, V_MAX);
    $('r-f').value = linv(st.f, F_MIN, F_MAX); $('r-f-val').textContent = num(st.f, st.f < 10 ? 1 : 0) + ' Hz';
    $('r-L').value = st.L; $('r-L-val').textContent = num(st.L, 0) + ' mm';
    $('r-v').value = linv(st.v, V_MIN, V_MAX); $('r-v-val').textContent = num(st.v, st.v < 10 ? 1 : 0) + ' m/s';
  }

  var W = 760, Lm = 58, Rm = 18;
  function medir() { var w = $('r-plot').getBoundingClientRect().width; W = Math.max(300, Math.round(w || 760)); Lm = W < 480 ? 46 : 58; Rm = W < 480 ? 12 : 18; }

  function dibujar() {
    var H = 280, B = 36, T = 14, h = H - B - T, T_ms = R.periodoMs(st.f), dur = 3 * T_ms;
    function X(t) { return Lm + t / dur * (W - Lm - Rm); }
    function Y(a) { return T + (1.25 - a) / 2.5 * h; }
    var tau = R.tauMs(st.L, st.v), tauRef = R.tauMs(st.L, V_MIEL), s = '', i, n = 400;
    [-1, 0, 1].forEach(function (a) {
      s += '<line class="rejilla" x1="' + Lm + '" x2="' + (W - Rm) + '" y1="' + Y(a) + '" y2="' + Y(a) + '"/><text x="' + (Lm - 6) + '" y="' + (Y(a) + 4) + '" text-anchor="end" font-size="12">' + num(a, 0) + '</text>';
    });
    var paso = dur > 1500 ? 500 : dur > 600 ? 100 : dur > 250 ? 50 : dur > 120 ? 20 : 10;
    for (var t = 0; t <= dur + 1e-9; t += paso) {
      s += '<line class="rejilla" x1="' + X(t) + '" x2="' + X(t) + '" y1="' + T + '" y2="' + (H - B) + '"/><text x="' + X(t) + '" y="' + (H - B + 16) + '" text-anchor="middle" font-size="12">' + num(t, 0) + '</text>';
    }
    s += '<text x="' + ((Lm + W - Rm) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" font-size="12">Tiempo (ms) — tres periodos de ' + num(st.f, st.f < 10 ? 1 : 0) + ' Hz</text>';
    s += '<text transform="translate(14 ' + (T + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Amplitud (u. arb.)</text>';
    s += '<rect x="' + Lm + '" y="' + T + '" width="' + (W - Lm - Rm) + '" height="' + h + '" fill="none" class="eje"/>';
    function curva(retardo, cls) {
      var d = '';
      for (i = 0; i <= n; i++) { var tt = dur * i / n; d += (i ? 'L' : 'M') + X(tt).toFixed(1) + ',' + Y(Math.sin(2 * Math.PI * st.f * (tt - retardo) / 1000)).toFixed(1) + ' '; }
      return '<path class="' + cls + '" d="' + d + '"/>';
    }
    s += curva(0, 'curva-0') + curva(tauRef, 'curva-1') + curva(tau, 'curva-2');
    var tm = tau % T_ms;                                  // primer cruce ascendente de B, módulo un periodo
    s += '<line class="guia" x1="' + X(tm) + '" x2="' + X(tm) + '" y1="' + T + '" y2="' + (H - B) + '"/><text x="' + (X(tm) + 4) + '" y="' + (T + 12) + '" font-size="12">B cruza 0 ↑</text>';
    s += '<line class="guia" x1="' + X(0) + '" x2="' + X(0) + '" y1="' + T + '" y2="' + (H - B) + '"/>';
    $('r-plot').setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    $('r-plot').innerHTML = s;
  }

  function lectura() {
    var tau = R.tauMs(st.L, st.v), T_ms = R.periodoMs(st.f), d = R.faseDeg(st.f, tau), w = R.envolver(d);
    var tauRef = R.tauMs(st.L, V_MIEL), dRef = R.faseDeg(st.f, tauRef), vc = R.velocidadCritica(st.f, st.L), p = [];
    p.push('<p><strong>Retardo:</strong> con L = ' + num(st.L, 0) + ' mm y v = ' + num(st.v, st.v < 10 ? 1 : 0) + ' m/s, τ = L/v = ' + num(tau, 1) + ' ms. A ' + num(st.f, st.f < 10 ? 1 : 0) + ' Hz el periodo es ' + num(T_ms, 1) + ' ms, así que τ equivale a un desfase de <strong>' + num(d, 1) + '°</strong> (' + num(d / 360, 2) + ' ciclos).</p>');
    var signo = Math.abs(w) < 0.05 ? 'las dos regiones parecen estar en fase' : (w > 0 ? 'B va retrasada respecto a A' : 'B parece ir <em>adelantada</em> respecto a A');
    p.push('<p><strong>Desfase observable</strong> (solo se puede medir entre −180° y 180°): <strong>' + num(w, 1) + '°</strong>; ' + signo + '.</p>');
    if (d > 180) p.push('<p>El retardo supera medio periodo: el desfase se «envuelve» y su signo cambia aunque B siga llegando <em>después</em> que A. A partir de ese punto el desfase medido por sí solo no permite saber cuántos ciclos completos hay de por medio.</p>');
    p.push('<p><strong>Referencia:</strong> con el valor mielinizado del informe (' + V_MIEL + ' m/s) el retardo sería ' + num(tauRef, 1) + ' ms y el desfase ' + num(dRef, 1) + '°. La velocidad a la que el desfase llega a 180° con estos f y L es v* = 2·f·L = ' + num(vc, 2) + ' m/s.</p>');
    $('r-lectura').innerHTML = p.join('');
  }

  function render() { medir(); escribir(); dibujar(); lectura(); }

  $('r-f').addEventListener('input', function () { st.f = lmap(Number($('r-f').value), F_MIN, F_MAX); render(); });
  $('r-L').addEventListener('input', function () { st.L = Number($('r-L').value); render(); });
  $('r-v').addEventListener('input', function () { st.v = lmap(Number($('r-v').value), V_MIN, V_MAX); render(); });
  $('r-miel').addEventListener('click', function () { st.v = V_MIEL; render(); });
  $('r-desm').addEventListener('click', function () { st.v = V_DESM; render(); });
  $('r-reset').addEventListener('click', function () { st = Object.assign({}, DEF); render(); });
  var tm = null; window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  render();
})();
