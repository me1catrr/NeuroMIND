/* Simulación «PLI y wPLI con épocas sintéticas» (interactive/estado-arte/wpli.js hace el cálculo). Datos sintéticos; sin red.
 * Tres vistas de una misma realización: (a) plano del espectro cruzado S_k, (b) voto del PLI frente a peso del wPLI por época,
 * (c) PLI, wPLI y coherencia² frente al desfase φ (misma semilla para todos los φ). Script clásico. */
(function () {
  'use strict';
  var W = window.NMWpli, raiz = document.getElementById('demo-wpli');
  if (!W || !raiz) return;
  var $ = function (id) { return document.getElementById(id); };

  var DEF = { N: 40, p: 100, fase: 60, sigma: 50, mezcla: 0, semilla: 1 };    // p en %, sigma y mezcla en centésimas
  var st = Object.assign({}, DEF);
  var LIM = { N: [2, 400], p: [0, 100], fase: [0, 180], sigma: [0, 200], mezcla: [0, 100], semilla: [1, 9999] };

  function clamp(x, lo, hi) { return Math.min(hi, Math.max(lo, x)); }
  function num(x, d) {
    var s = Math.abs(x) < Math.pow(10, -d) / 2 ? (0).toFixed(d) : x.toFixed(d);
    return s.replace('-', '−').replace('.', ',');
  }
  function par() { return { N: st.N, p: st.p / 100, fase: st.fase, sigma: st.sigma / 100, mezcla: st.mezcla / 100, semilla: st.semilla }; }

  function escribir() {
    Object.keys(LIM).forEach(function (k) { st[k] = clamp(Math.round(st[k]), LIM[k][0], LIM[k][1]); });
    $('w-N').value = st.N; $('w-N-val').textContent = st.N;
    $('w-p').value = st.p; $('w-p-val').textContent = st.p + ' %';
    $('w-fase').value = st.fase; $('w-fase-val').textContent = st.fase + '°';
    $('w-sigma').value = st.sigma; $('w-sigma-val').textContent = num(st.sigma / 100, 2);
    $('w-mezcla').value = st.mezcla; $('w-mezcla-val').textContent = num(st.mezcla / 100, 2);
    $('w-semilla').value = st.semilla;
  }

  // Elementos SVG como texto (se asignan con innerHTML, como las otras simulaciones del libro)
  function esc(t) { return String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;'); }
  function el(tag, attrs, texto) {
    var a = '';
    for (var k in attrs) a += ' ' + k + '="' + esc(attrs[k]).replace(/"/g, '&quot;') + '"';
    return texto == null ? '<' + tag + a + '/>' : '<' + tag + a + '>' + esc(texto) + '</' + tag + '>';
  }
  function ancho(svg) { var w = svg.getBoundingClientRect().width; return Math.max(300, Math.round(w || 760)); }

  // (a) plano del espectro cruzado
  function dibujarPlano(s) {
    var svg = $('w-plano'), w = ancho(svg), h = 300; svg.setAttribute('viewBox', '0 0 ' + w + ' ' + h); var P = [];
    var cx = w / 2, cy = h / 2, r = Math.min(w, h) / 2 - 28, mx = 1e-12;
    for (var k = 0; k < s.N; k++) mx = Math.max(mx, Math.abs(s.reS[k]), Math.abs(s.imS[k]));
    var f = r / mx;
    P.push(el('line', { 'class': 'eje', x1: cx - r, y1: cy, x2: cx + r, y2: cy }));
    P.push(el('line', { 'class': 'eje', x1: cx, y1: cy - r, x2: cx, y2: cy + r }));
    P.push(el('text', { x: cx + r - 2, y: cy - 6, 'text-anchor': 'end', 'font-size': 12 }, 'Re S'));
    P.push(el('text', { x: cx + 6, y: cy - r + 12, 'font-size': 12 }, 'Im S (↑ X adelanta a Y)'));
    P.push(el('text', { x: cx + 6, y: cy + r - 4, 'font-size': 12 }, '(↓ X retrasa)'));
    P.push(el('text', { x: cx - r, y: cy + r + 18, 'font-size': 11 }, 'escala: ±' + num(mx, 2) + ' (unidades de potencia)'));
    for (var i = 0; i < s.N; i++) {
      var x = cx + s.reS[i] * f, y = cy - s.imS[i] * f;
      var c = s.acoplada[i]
        ? el('circle', { cx: x, cy: y, r: 4, fill: 'var(--nm-petrol)', opacity: 0.85 })
        : el('circle', { cx: x, cy: y, r: 3.5, fill: 'none', stroke: 'var(--nm-muted)', 'stroke-width': 1.3 });
      P.push(c);
    }
    svg.innerHTML = P.join('');
  }

  // (b) voto del PLI frente a peso del wPLI
  function dibujarPesos(s) {
    var svg = $('w-pesos'), w = ancho(svg), h = 256; svg.setAttribute('viewBox', '0 0 ' + w + ' ' + h); var P = [];
    if (s.N > 60) {
      P.push(el('text', { x: w / 2, y: h / 2, 'text-anchor': 'middle', 'font-size': 13 }, 'Con más de 60 épocas no se dibuja cada una: reduce N para verlas.'));
      svg.innerHTML = P.join('');
      return;
    }
    var Lm = 56, Rm = 12, anchoUtil = w - Lm - Rm, bw = anchoUtil / s.N, mx = 1e-12;
    for (var k = 0; k < s.N; k++) mx = Math.max(mx, Math.abs(s.imS[k]));
    var y1 = 74, y2 = 188, a = 44;                 // ejes de las dos filas y semialtura máxima
    P.push(el('text', { x: 4, y: y1 - a - 8, 'font-size': 12 }, 'Voto del PLI: cada época vale +1 o −1'));
    P.push(el('text', { x: 4, y: y2 - a - 8, 'font-size': 12 }, 'Peso del wPLI: |Im S| de cada época'));
    P.push(el('line', { 'class': 'eje', x1: Lm, y1: y1, x2: w - Rm, y2: y1 }));
    P.push(el('line', { 'class': 'eje', x1: Lm, y1: y2, x2: w - Rm, y2: y2 }));
    for (var i = 0; i < s.N; i++) {
      var v = s.imS[i], sg = v > 0 ? 1 : (v < 0 ? -1 : 0), x = Lm + i * bw + bw * 0.1, b = Math.max(1, bw * 0.8);
      var col = v >= 0 ? 'var(--nm-petrol)' : 'var(--nm-orange)';
      if (sg) P.push(el('rect', { x: x, y: sg > 0 ? y1 - a : y1, width: b, height: a, fill: col, opacity: 0.85 }));
      var hh = Math.abs(v) / mx * a;
      if (hh > 0) P.push(el('rect', { x: x, y: v > 0 ? y2 - hh : y2, width: b, height: hh, fill: col, opacity: 0.85 }));
    }
    P.push(el('text', { x: Lm - 4, y: y1 - a + 8, 'text-anchor': 'end', 'font-size': 11 }, '+1'));
    P.push(el('text', { x: Lm - 4, y: y1 + a - 2, 'text-anchor': 'end', 'font-size': 11 }, '−1'));
    P.push(el('text', { x: Lm - 4, y: y2 - a + 8, 'text-anchor': 'end', 'font-size': 11 }, '+' + num(mx, 1)));
    P.push(el('text', { x: Lm - 4, y: y2 + a - 2, 'text-anchor': 'end', 'font-size': 11 }, '−' + num(mx, 1)));
    P.push(el('text', { x: w - Rm, y: h - 6, 'text-anchor': 'end', 'font-size': 11 }, 'épocas 1…' + s.N + ' · azul: Im S > 0; naranja: Im S < 0'));
    svg.innerHTML = P.join('');
  }

  // (c) curvas frente a φ
  function dibujarCurva() {
    var svg = $('w-curva'), w = ancho(svg), h = 250; svg.setAttribute('viewBox', '0 0 ' + w + ' ' + h); var P = [];
    var Lm = 48, Rm = 14, T = 14, B = 34, pw = w - Lm - Rm, ph = h - T - B;
    var X = function (fase) { return Lm + fase / 180 * pw; }, Y = function (v) { return T + (1 - v) * ph; };
    for (var g = 0; g <= 4; g++) {
      var v = g / 4;
      P.push(el('line', { 'class': 'rejilla', x1: Lm, y1: Y(v), x2: w - Rm, y2: Y(v) }));
      P.push(el('text', { x: Lm - 6, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 11 }, num(v, 2)));
    }
    [0, 45, 90, 135, 180].forEach(function (d) {
      P.push(el('text', { x: X(d), y: h - 14, 'text-anchor': 'middle', 'font-size': 11 }, d + '°'));
    });
    P.push(el('text', { x: Lm + pw / 2, y: h - 1, 'text-anchor': 'middle', 'font-size': 11 }, 'desfase φ de la fuente que interactúa'));
    var datos = W.curvaFase(par(), 36);
    function trazo(clave, cls) {
      var d = datos.map(function (q, i) { return (i ? 'L' : 'M') + X(q.fase).toFixed(1) + ' ' + Y(q[clave]).toFixed(1); }).join(' ');
      P.push(el('path', { 'class': cls, d: d }));
    }
    trazo('coherencia2', 'curva-0'); trazo('pli', 'curva-1'); trazo('wpli', 'curva-2');
    P.push(el('line', { 'class': 'guia', x1: X(st.fase), y1: T, x2: X(st.fase), y2: T + ph }));
    svg.innerHTML = P.join('');
  }

  function lectura(s, m) {
    var acopl = 0; for (var k = 0; k < s.N; k++) acopl += s.acoplada[k];
    var t = 'N = ' + s.N + ' épocas (' + acopl + ' acopladas). ' +
      'PLI = ' + num(m.pli, 3) + ' · wPLI = ' + num(m.wpli, 3) + ' · dwPLI² = ' + num(m.dwpli2, 3) + ' · coherencia² = ' + num(m.coherencia2, 3) + '. ';
    if (st.fase < 3 || st.fase > 177) t += 'Con φ ≈ 0° o 180° la interacción real es invisible para PLI y wPLI (Im S ≈ 0): no es que no exista, es que estas métricas no pueden verla. ';
    if (st.mezcla > 0 && st.p === 0) t += 'Solo hay mezcla instantánea: la coherencia² sube y PLI y wPLI no deberían, salvo el ruido de una muestra pequeña. ';
    if (st.N < 15) t += 'Con pocas épocas las cifras fluctúan mucho entre semillas: cambia la semilla para verlo. ';
    $('w-lectura').textContent = t;
  }

  function dibujar() {
    escribir();
    var s = W.simular(par()), m = W.metricas(s.imS, s);
    dibujarPlano(s); dibujarPesos(s); dibujarCurva(); lectura(s, m);
  }

  [['w-N', 'N'], ['w-p', 'p'], ['w-fase', 'fase'], ['w-sigma', 'sigma'], ['w-mezcla', 'mezcla'], ['w-semilla', 'semilla']].forEach(function (c) {
    $(c[0]).addEventListener('input', function () { var v = parseFloat(this.value); if (isFinite(v)) { st[c[1]] = v; dibujar(); } });
  });
  $('w-reset').addEventListener('click', function () { st = Object.assign({}, DEF); dibujar(); });
  $('w-solo-mezcla').addEventListener('click', function () { st = Object.assign({}, DEF, { p: 0, mezcla: 100, sigma: 30 }); dibujar(); });
  $('w-cero').addEventListener('click', function () { st = Object.assign({}, DEF, { fase: 0, sigma: 30 }); dibujar(); });
  window.addEventListener('resize', function () { clearTimeout(dibujar._t); dibujar._t = setTimeout(dibujar, 120); });
  dibujar();
})();
