/* Simulación A — respuesta en frecuencia de un Butterworth digital (una pasada frente a doble pasada).
 * Datos sintéticos: solo usa la fórmula de interactive/filtrado/filtros.js (validada contra DSP.jl en tests/verificar.mjs).
 * Script clásico, sin red ni dependencias; no lee ni escribe nada del proyecto. */
(function () {
  'use strict';
  var F = window.NMFiltros;
  var raiz = document.getElementById('demo-a');
  if (!F || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };

  var FMIN = 0.01;                    // Hz, extremo inferior del eje
  var YMIN = -120, YMAX = 6;          // dB del gráfico de magnitud
  var DEF = { tipo: 'hp', fs: 500, fc: 0.5, f0: 50, bw: 1, N: 4, fp: 0.25 };
  var DEF_TIPO = { hp: { fc: 0.5, fp: 0.25 }, lp: { fc: 150, fp: 100 }, bs: { f0: 50, bw: 1, fp: 49 } };
  var st = Object.assign({}, DEF);

  function nyq() { return st.fs / 2; }
  function lmap(v, lo, hi) { return lo * Math.pow(hi / lo, v / 1000); }
  function linv(x, lo, hi) { return Math.round(1000 * Math.log(x / lo) / Math.log(hi / lo)); }
  function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }
  function rFc() { return [0.05, 0.98 * nyq()]; }
  function rF0() { return [0.5, 0.98 * nyq()]; }
  function rBw() { return [0.1, Math.max(0.2, Math.min(1.9 * st.f0, 2 * (0.995 * nyq() - st.f0)))]; }
  function rFp() { return [FMIN, 0.999 * nyq()]; }

  // Formato con coma decimal y signo menos tipográfico.
  function num(x, d) {
    var s = Math.abs(x) < 1e-12 ? '0' : x.toFixed(d);
    if (/^-0(\.0+)?$/.test(s)) s = s.slice(1);
    return s.replace('-', '−').replace('.', ',');
  }
  function sig(x) {                   // 3 cifras significativas para frecuencias
    var d = x >= 100 ? 0 : x >= 10 ? 1 : x >= 1 ? 2 : x >= 0.1 ? 2 : 3;
    return num(x, d);
  }

  function clampEstado() {
    st.N = clamp(Math.round(st.N), F.N_MIN, F.N_MAX);
    var a = rFc(); st.fc = clamp(st.fc, a[0], a[1]);
    a = rF0(); st.f0 = clamp(st.f0, a[0], a[1]);
    a = rBw(); st.bw = clamp(st.bw, a[0], a[1]);
    a = rFp(); st.fp = clamp(st.fp, a[0], a[1]);
  }

  function escribirControles() {
    var a;
    $('a-tipo').value = st.tipo; $('a-fs').value = String(st.fs);
    a = rFc(); $('a-fc').value = linv(st.fc, a[0], a[1]); $('a-fc-val').textContent = sig(st.fc) + ' Hz';
    a = rF0(); $('a-f0').value = linv(st.f0, a[0], a[1]); $('a-f0-val').textContent = sig(st.f0) + ' Hz';
    a = rBw(); $('a-bw').value = linv(st.bw, a[0], a[1]); $('a-bw-val').textContent = sig(st.bw) + ' Hz';
    $('a-n').value = st.N; $('a-n-val').textContent = String(st.N);
    a = rFp(); $('a-fp').value = linv(st.fp, a[0], a[1]); $('a-fp-val').textContent = sig(st.fp) + ' Hz';
    var bs = st.tipo === 'bs';
    $('a-fc-etq').hidden = bs; $('a-f0-etq').hidden = !bs; $('a-bw-etq').hidden = !bs;
    $('a-lim').textContent = 'Límites con Fs = ' + st.fs + ' Hz (Nyquist = ' + nyq() + ' Hz): ' +
      (bs ? '0 < f0 − B/2 y f0 + B/2 < Nyquist' : '0 < fc < Nyquist') +
      '; orden de diseño entero entre 1 y 8. Los deslizadores no admiten valores fuera de rango.';
  }

  function cfg() {
    return st.tipo === 'bs'
      ? { tipo: 'bs', fs: st.fs, N: st.N, f0: st.f0, bw: st.bw }
      : { tipo: st.tipo, fs: st.fs, N: st.N, fc: st.fc };
  }

  // ───── dibujo ─────
  // El ancho del dibujo es el real del contenedor (1 unidad = 1 píxel): así el texto mantiene su tamaño en ventanas estrechas.
  var W = 760, L = 58, R = 18;
  function medir() {
    var w = $('a-mag').getBoundingClientRect().width;
    W = Math.max(300, Math.round(w || 760)); L = W < 480 ? 46 : 58; R = W < 480 ? 12 : 18;
  }
  function xpix(f) { return L + (Math.log10(f) - Math.log10(FMIN)) / (Math.log10(nyq()) - Math.log10(FMIN)) * (W - L - R); }

  function ejeX(H, B, etiqueta) {
    var s = '', k, f;
    for (k = -2; k <= 3; k++) {
      f = Math.pow(10, k);
      if (f < FMIN || f > nyq()) continue;
      s += '<line class="rejilla" x1="' + xpix(f) + '" x2="' + xpix(f) + '" y1="12" y2="' + (H - B) + '"/>';
      // la etiqueta de la década se omite si choca con la de Nyquist (se funde con ella)
      if (xpix(nyq()) - xpix(f) > 90) s += '<text x="' + xpix(f) + '" y="' + (H - B + 16) + '" text-anchor="middle" font-size="12">' + num(f, 0 > k ? -k : 0) + '</text>';
    }
    s += '<line class="rejilla" x1="' + xpix(nyq()) + '" x2="' + xpix(nyq()) + '" y1="12" y2="' + (H - B) + '"/>';
    var dec = Math.pow(10, Math.floor(Math.log10(nyq() * 0.999)));
    var pre = xpix(nyq()) - xpix(dec) <= 90 && dec >= FMIN ? num(dec, 0) + ' · ' : '';
    s += '<text x="' + xpix(nyq()) + '" y="' + (H - B + 16) + '" text-anchor="end" font-size="12">' + pre + 'Nyquist ' + nyq() + '</text>';
    if (etiqueta) s += '<text x="' + ((L + W - R) / 2) + '" y="' + (H - 4) + '" text-anchor="middle" font-size="12">Frecuencia (Hz, escala logarítmica)</text>';
    return s;
  }

  function marcas(H, B) {
    var s = '', bs = st.tipo === 'bs', top = 22;
    function v(f, txt, dash) {
      if (f <= FMIN || f >= nyq()) return '';
      return '<line class="guia" x1="' + xpix(f) + '" x2="' + xpix(f) + '" y1="12" y2="' + (H - B) + '"' + (dash ? ' style="stroke-dasharray:2 3"' : '') + '/>' +
        '<text x="' + (xpix(f) + 4) + '" y="' + top + '" font-size="12">' + txt + '</text>';
    }
    if (bs) { s += v(st.f0 - st.bw / 2, 'f₁') + v(st.f0 + st.bw / 2, 'f₂'); }
    else { s += v(st.fc, 'fc'); }
    s += '<line x1="' + xpix(st.fp) + '" x2="' + xpix(st.fp) + '" y1="12" y2="' + (H - B) + '" style="stroke:var(--nm-plum);stroke-width:1.5;stroke-dasharray:1 4"/>' +
      '<text x="' + (xpix(st.fp) + 4) + '" y="' + (top + 14) + '" font-size="12" style="fill:var(--nm-t-plum)">f*</text>';
    return s;
  }

  function pasos(xs, ys, corte) {     // construye path; corte(i) true → empieza segmento nuevo
    var d = '';
    for (var i = 0; i < xs.length; i++) d += (i === 0 || (corte && corte(i)) ? 'M' : 'L') + xs[i].toFixed(1) + ',' + ys[i].toFixed(1) + ' ';
    return d;
  }

  function dibujarMag(c, hz) {
    var H = 330, B = 34, T = 12, h = H - B - T;
    function y(db) { return T + (YMAX - clamp(db, YMIN, YMAX)) / (YMAX - YMIN) * h; }
    var s = '', db, i;
    for (db = YMIN; db <= 0; db += 20) {
      s += '<line class="rejilla" x1="' + L + '" x2="' + (W - R) + '" y1="' + y(db) + '" y2="' + y(db) + '"/>' +
        '<text x="' + (L - 6) + '" y="' + (y(db) + 4) + '" text-anchor="end" font-size="12">' + num(db, 0) + '</text>';
    }
    s += '<text transform="translate(14 ' + (T + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Magnitud (dB)</text>';
    s += ejeX(H, B, false);
    [[-3.0103, -3], [-6.0206, 12]].forEach(function (g) {    // la 2.ª etiqueta va bajo su línea para no pisar a la 1.ª
      s += '<line class="guia" x1="' + L + '" x2="' + (W - R) + '" y1="' + y(g[0]) + '" y2="' + y(g[0]) + '"/>' +
        '<text x="' + (W - R - 4) + '" y="' + (y(g[0]) + g[1]) + '" text-anchor="end" font-size="11">' + num(g[0], 2) + ' dB</text>';
    });
    s += '<rect x="' + L + '" y="' + T + '" width="' + (W - L - R) + '" height="' + h + '" fill="none" class="eje"/>';
    var xs = hz.map(xpix);
    s += '<path class="curva-1" d="' + pasos(xs, c.mag1.map(function (m) { return y(F.dB(m)); })) + '"/>';
    s += '<path class="curva-2" d="' + pasos(xs, c.mag2.map(function (m) { return y(F.dB(m)); })) + '"/>';
    s += marcas(H, B);
    var h1 = F.respuesta(cfg(), st.fp), m1 = F.mag(h1);
    s += '<circle cx="' + xpix(st.fp) + '" cy="' + y(F.dB(m1)) + '" r="5" style="fill:var(--nm-orange)"/>' +
      '<circle cx="' + xpix(st.fp) + '" cy="' + y(2 * F.dB(m1)) + '" r="5" style="fill:var(--nm-petrol)"/>';
    $('a-mag').setAttribute('viewBox', '0 0 ' + W + ' ' + H); $('a-mag').innerHTML = s;
  }

  function dibujarFase(c, hz) {
    var H = 210, B = 34, T = 12, h = H - B - T;
    function y(deg) { return T + (180 - deg) / 360 * h; }
    var s = '', k;
    [-180, -90, 0, 90, 180].forEach(function (d) {
      s += '<line class="rejilla" x1="' + L + '" x2="' + (W - R) + '" y1="' + y(d) + '" y2="' + y(d) + '"/>' +
        '<text x="' + (L - 6) + '" y="' + (y(d) + 4) + '" text-anchor="end" font-size="12">' + num(d, 0) + '°</text>';
    });
    s += '<text transform="translate(14 ' + (T + h / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">Fase (grados)</text>';
    s += ejeX(H, B, true);
    s += '<rect x="' + L + '" y="' + T + '" width="' + (W - L - R) + '" height="' + h + '" fill="none" class="eje"/>';
    var xs = hz.map(xpix), deg = c.fase1.map(function (f) { return f * 180 / Math.PI; });
    s += '<path class="curva-1" d="' + pasos(xs, deg.map(y), function (i) { return Math.abs(deg[i] - deg[i - 1]) > 180; }) + '"/>';
    s += '<line class="curva-2" x1="' + L + '" x2="' + (W - R) + '" y1="' + y(0) + '" y2="' + y(0) + '"/>';
    s += '<text x="' + (L + 6) + '" y="' + (y(0) - 6) + '" font-size="12" style="fill:var(--nm-t-petrol)">doble pasada: fase 0° en todas las frecuencias</text>';
    s += marcas(H, B);
    $('a-fase').setAttribute('viewBox', '0 0 ' + W + ' ' + H); $('a-fase').innerHTML = s;
  }

  function lectura() {
    var C = cfg(), v = F.validar(C), p = [], a;
    if (!v.ok) { $('a-lectura').innerHTML = '<p>' + v.mensaje + '</p>'; return; }
    var h = F.respuesta(C, st.fp), m = F.mag(h), ph = Math.atan2(h[1], h[0]) * 180 / Math.PI;
    var d1 = F.dB(m), d2 = 2 * d1;
    function dbTxt(x) { return x <= YMIN ? 'menos de −120 dB (recortado en el gráfico)' : num(x, 2) + ' dB'; }
    p.push('<p><strong>En f* = ' + sig(st.fp) + ' Hz:</strong> una pasada ' + dbTxt(d1) + ', fase ' + num(ph, 1) + '°; doble pasada ' + dbTxt(d2) + ', fase 0°.</p>');
    if (st.tipo === 'bs') {
      var hc = F.mag(F.respuesta(C, st.f0));
      p.push('<p><strong>Bordes de la banda</strong> f₁ = ' + sig(st.f0 - st.bw / 2) + ' Hz y f₂ = ' + sig(st.f0 + st.bw / 2) + ' Hz: −3,01 dB con una pasada y −6,02 dB con doble pasada. En el centro f₀ = ' + sig(st.f0) + ' Hz la ganancia es de ' + num(F.dB(hc), 0) + ' dB con una pasada (' + hc.toExponential(1).replace('.', ',') + ' en amplitud).</p>');
      p.push('<p>Una pasada tiene ' + (2 * st.N) + ' polos (' + st.N + ' secciones de segundo orden), aunque el prototipo de diseño es de orden ' + st.N + '. Con doble pasada la magnitud es |H|², pero la fase es cero.</p>');
    } else {
      p.push('<p><strong>En el corte</strong> fc = ' + sig(st.fc) + ' Hz: −3,01 dB con una pasada y −6,02 dB con doble pasada.</p>');
      p.push('<p>Pendiente asintótica lejos del corte: ' + (20 * st.N) + ' dB/década con una pasada y ' + (40 * st.N) + ' dB/década con doble pasada (orden de diseño ' + st.N + '; una pasada tiene ' + st.N + ' polos).</p>');
    }
    $('a-lectura').innerHTML = p.join('');
  }

  function render() {
    medir();
    clampEstado();
    escribirControles();
    var C = cfg();
    if (!F.validar(C).ok) { lectura(); return; }
    var hz = F.rejillaLog(FMIN, 0.9995 * nyq(), 700), c = F.curva(C, hz);
    dibujarMag(c, hz); dibujarFase(c, hz); lectura();
  }

  // ───── eventos ─────
  function leer(id, fn) { $(id).addEventListener('input', function () { fn(Number($(id).value)); render(); }); }
  $('a-tipo').addEventListener('change', function () {
    st.tipo = $('a-tipo').value;
    Object.assign(st, DEF_TIPO[st.tipo]); render();
  });
  $('a-fs').addEventListener('change', function () { st.fs = Number($('a-fs').value); render(); });
  leer('a-fc', function (v) { var a = rFc(); st.fc = lmap(v, a[0], a[1]); });
  leer('a-f0', function (v) { var a = rF0(); st.f0 = lmap(v, a[0], a[1]); });
  leer('a-bw', function (v) { var a = rBw(); st.bw = lmap(v, a[0], a[1]); });
  leer('a-n', function (v) { st.N = v; });
  leer('a-fp', function (v) { var a = rFp(); st.fp = lmap(v, a[0], a[1]); });
  $('a-reset').addEventListener('click', function () { st = Object.assign({}, DEF); render(); });

  var tm = null;
  window.addEventListener('resize', function () { clearTimeout(tm); tm = setTimeout(render, 120); });
  render();
})();
