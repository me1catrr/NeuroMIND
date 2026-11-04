/* Núcleo numérico de la Simulación C («filt frente a filtfilt»): diseño Butterworth en secciones de segundo orden (biquads),
 * `filt` causal y `filtfilt` de fase cero tal como los implementa DSP.jl 0.8.5, y retardo de grupo/fase de una pasada.
 * Script clásico (file:// y Node), sin dependencias ni red. Reutiliza interactive/filtrado/filtros.js (límites y respuesta analítica).
 *
 * Diseño (el del pipeline, src/preprocessing/Filtering.jl → _build: digitalfilter(Highpass|Lowpass(fc/nyq), Butterworth(N))):
 *   polos analógicos del prototipo de corte 1: p_k = exp[ jπ(1/2 + (2k−1)/(2N)) ];  Ωc = tan(π·fc/Fs)  (bilineal con prewarping)
 *   paso bajo:  s = Ωc·p      paso alto:  s = Ωc / p        polo digital:  z = (1 + s) / (1 − s)
 *   ceros digitales: paso bajo en z = −1 (N veces), paso alto en z = +1 (N veces).
 * Secciones, como DSP.jl (SecondOrderSections{:z}(ZeroPoleGain)): cada par de polos conjugados (a1 = −2·Re z, a2 = |z|²) con dos
 * ceros → un biquad con b0 = 1; si N es impar, un polo real con un cero → biquad de primer orden (b2 = a2 = 0) que va PRIMERO
 * de la cascada; después, los pares desde el más alejado de la circunferencia unidad hasta el más próximo. Nº de biquads = ⌈N/2⌉.
 * La ganancia g se fija para que |H| = 1 en la frecuencia de referencia (z = 1 en paso bajo; z = −1 en paso alto); DSP.jl usa la
 * fórmula del producto de polos y ceros, que da lo mismo salvo redondeo (~1e-16).
 *
 * filt: forma directa II transpuesta, estado inicial CERO, y la salida de la cascada se multiplica por g (como _filt! de DSP.jl).
 * filtfilt (DSP.jl 0.8.5, método de biquads): relleno de min(6·nº_biquads, n−1) muestras por extremo con extensión de simetría impar
 *   ext = [2·x[0] − x[pad..1], x, 2·x[n−1] − x[n−2..n−1−pad]],
 * estado inicial ESTACIONARIO de cada biquad (respuesta a un escalón unidad) multiplicado por la primera muestra de la señal
 * extendida, pasada hacia delante, inversión, misma operación con el primer valor de la señal invertida y recorte del relleno.
 *
 * Retardo de grupo (muestras) de un biquad B/A con z⁻¹ = e^{−jω}:  τ = Re[ Σ k·b_k e^{−jωk} / B ] − Re[ Σ k·a_k e^{−jωk} / A ];
 * el del filtro es la suma de sus secciones. Retardo de fase: τ_f = −φ/ω. filtfilt aplica H(z)·H(1/z) = |H|²: fase 0 y retardo 0.
 * Validación: tests/verificar.mjs (sección N) contra DSP.jl (scripts/generar_simulacion_filtfilt.jl → tests/fixtures/ref_filtfilt.json). */
(function (root) {
  'use strict';
  var F = (typeof module === 'object' && module.exports) ? require('./filtros.js') : root.NMFiltros;

  var FS_SENAL = 500;                 // Hz: la simulación usa la frecuencia de muestreo del dataset
  var TIPOS = ['hp', 'lp'];

  /** Polos analógicos del prototipo de corte 1 (semiplano izquierdo). */
  function polosProto(N) {
    var p = [];
    for (var k = 1; k <= N; k++) {
      var a = Math.PI / 2 + (2 * k - 1) * Math.PI / (2 * N);
      p.push([Math.cos(a), Math.sin(a)]);
    }
    return p;
  }

  /** Diseña el filtro. cfg = {tipo:'hp'|'lp', fs, N, fc}. Devuelve {tipo, fs, N, fc, sos:[{b0,b1,b2,a1,a2}], g, nb, polos:[[re,im]...]}. */
  function disenar(cfg) {
    var v = F.validar({ tipo: cfg.tipo, fs: cfg.fs, N: cfg.N, fc: cfg.fc });
    if (TIPOS.indexOf(cfg.tipo) < 0) throw new Error('tipo no admitido en la Simulación C: ' + cfg.tipo);
    if (!v.ok) throw new Error(v.mensaje);
    var N = cfg.N, Oc = Math.tan(Math.PI * cfg.fc / cfg.fs), hp = cfg.tipo === 'hp', zs = hp ? 1 : -1;
    var pares = [], real = null, todos = [], pp = polosProto(N);
    for (var k = 0; k < N; k++) {
      var pr = pp[k][0], pi = pp[k][1], sr, si;
      if (hp) {                                   // s = Ωc / p
        var m = pr * pr + pi * pi; sr = Oc * pr / m; si = -Oc * pi / m;
      } else { sr = Oc * pr; si = Oc * pi; }      // s = Ωc · p
      var dr = 1 - sr, dn = dr * dr + si * si;    // z = (1 + s)/(1 − s)
      var zr = ((1 + sr) * dr - si * si) / dn, zi = (si * dr + (1 + sr) * si) / dn;
      todos.push([zr, zi]);
      if (Math.abs(pi) < 1e-12) real = zr;
      else if (pi > 0) pares.push([zr, zi]);
    }
    // pares desde el más alejado de la circunferencia unidad hasta el más próximo (orden de la cascada de DSP.jl)
    pares.sort(function (a, b) { return Math.abs(Math.hypot(b[0], b[1]) - 1) - Math.abs(Math.hypot(a[0], a[1]) - 1); });
    var sos = [];
    if (N % 2 === 1) sos.push({ b0: 1, b1: -zs, b2: 0, a1: -real, a2: 0 });
    pares.forEach(function (z) { sos.push({ b0: 1, b1: -2 * zs, b2: 1, a1: -2 * z[0], a2: z[0] * z[0] + z[1] * z[1] }); });
    // ganancia: |H| = 1 en la frecuencia de referencia (z = zs·(−1)): paso bajo z = 1, paso alto z = −1
    var zref = hp ? -1 : 1, g = 1;
    sos.forEach(function (s) {
      g /= (s.b0 + s.b1 * zref + s.b2 * zref * zref) / (1 + s.a1 * zref + s.a2 * zref * zref);
    });
    return { tipo: cfg.tipo, fs: cfg.fs, N: N, fc: cfg.fc, sos: sos, g: g, nb: sos.length, polos: todos };
  }

  /** Aplicación de la cascada a x con estado inicial si[2·nb] (se modifica); escribe en out (puede ser x). */
  function aplicar(d, x, out, si) {
    var n = x.length, nb = d.nb, sos = d.sos, g = d.g;
    for (var i = 0; i < n; i++) {
      var yi = x[i];
      for (var f = 0; f < nb; f++) {
        var q = sos[f], xi = yi;
        yi = q.b0 * xi + si[2 * f];
        si[2 * f] = q.a1 * -yi + (q.b1 * xi + si[2 * f + 1]);
        si[2 * f + 1] = q.b2 * xi - q.a2 * yi;
      }
      out[i] = g * yi;
    }
  }

  /** filt causal con estado inicial cero. */
  function filt(d, x) {
    var out = new Float64Array(x.length);
    aplicar(d, x, out, new Float64Array(2 * d.nb));
    return out;
  }

  /** Estado estacionario de cada biquad ante un escalón unidad (filt_stepstate de DSP.jl para SecondOrderSections). */
  function estadoEstacionario(d) {
    var zi = new Float64Array(2 * d.nb), y = 1;
    for (var f = 0; f < d.nb; f++) {
      var q = d.sos[f], den = 1 + q.a1 + q.a2;
      zi[2 * f] = ((q.a1 + q.a2) * -q.b0 + (q.b1 + q.b2)) / den * y;
      zi[2 * f + 1] = (q.a1 * q.b2 + (-q.a2 * (q.b0 + q.b1) + q.b2)) / den * y;
      y *= (q.b0 + q.b1 + q.b2) / den;
    }
    return zi;
  }

  /** Longitud del relleno de filtfilt para una señal de n muestras: min(6·nb, n−1). */
  function relleno(d, n) { return Math.min(6 * d.nb, n - 1); }

  /** filtfilt de fase cero con extensión impar y estado inicial estacionario (DSP.jl 0.8.5). */
  function filtfilt(d, x) {
    if (x.length === 0) return new Float64Array(0);      // caso límite: DSP.jl falla con un registro vacío; aquí se devuelve vacío
    var n = x.length, pad = relleno(d, n), L = n + 2 * pad, ext = new Float64Array(L), zi = estadoEstacionario(d), si = new Float64Array(2 * d.nb), i;
    for (i = 0; i < pad; i++) ext[i] = 2 * x[0] - x[pad - i];
    for (i = 0; i < n; i++) ext[pad + i] = x[i];
    for (i = 1; i <= pad; i++) ext[n + pad + i - 1] = 2 * x[n - 1] - x[n - 1 - i];
    function pasada() {
      for (var k = 0; k < si.length; k++) si[k] = zi[k] * ext[0];
      aplicar(d, ext, ext, si);
    }
    pasada();
    ext.reverse();
    pasada();
    var out = new Float64Array(n);
    for (i = 0; i < n; i++) out[i] = ext[n + pad - 1 - i];
    return out;
  }

  /** Respuesta compleja de una pasada, evaluada con la cascada (independiente de la fórmula analítica de filtros.js). */
  function respuestaSOS(d, hz) {
    var w = 2 * Math.PI * hz / d.fs, c = Math.cos(w), s = -Math.sin(w), hr = d.g, hi = 0;   // e^{−jω} = c + j·s
    var c2 = Math.cos(2 * w), s2 = -Math.sin(2 * w);
    d.sos.forEach(function (q) {
      var br = q.b0 + q.b1 * c + q.b2 * c2, bi = q.b1 * s + q.b2 * s2;
      var ar = 1 + q.a1 * c + q.a2 * c2, ai = q.a1 * s + q.a2 * s2;
      var m = ar * ar + ai * ai, rr = (br * ar + bi * ai) / m, ri = (bi * ar - br * ai) / m;
      var t = hr * rr - hi * ri; hi = hr * ri + hi * rr; hr = t;
    });
    return [hr, hi];
  }

  /** Retardo de grupo de UNA pasada, en segundos (fórmula cerrada por secciones). */
  function retardoGrupo(d, hz) {
    var w = 2 * Math.PI * hz / d.fs, c1 = Math.cos(w), s1 = -Math.sin(w), c2 = Math.cos(2 * w), s2 = -Math.sin(2 * w), suma = 0;
    d.sos.forEach(function (q) {
      var br = q.b0 + q.b1 * c1 + q.b2 * c2, bi = q.b1 * s1 + q.b2 * s2, dbr = q.b1 * c1 + 2 * q.b2 * c2, dbi = q.b1 * s1 + 2 * q.b2 * s2;
      var ar = 1 + q.a1 * c1 + q.a2 * c2, ai = q.a1 * s1 + q.a2 * s2, dar = q.a1 * c1 + 2 * q.a2 * c2, dai = q.a1 * s1 + 2 * q.a2 * s2;
      suma += (dbr * br + dbi * bi) / (br * br + bi * bi) - (dar * ar + dai * ai) / (ar * ar + ai * ai);
    });
    return suma / d.fs;
  }

  /** Fase de una pasada en radianes (valor principal, −π..π). */
  function fase(d, hz) { var h = respuestaSOS(d, hz); return Math.atan2(h[1], h[0]); }

  /** Retardo de fase de una pasada, en segundos: τ_f = −φ/ω con φ ENVUELTA a (−π, π] (positivo = la salida llega después). */
  function retardoFase(d, hz) { return -fase(d, hz) / (2 * Math.PI * hz); }

  /** Un diseño es estable si todos sus polos tienen módulo < 1 (criterio de Jury para 1 + a1 z⁻¹ + a2 z⁻²). */
  function estable(d) {
    return d.sos.every(function (q) { return Math.abs(q.a2) < 1 && Math.abs(q.a1) < 1 + q.a2; });
  }

  // ───── señales sintéticas de la simulación ─────
  var A1 = 10, A2 = 20, A_ESCALON = 15, A_IMPULSO = 100, T_EVENTO = 4;        // µV y s (declarados en el texto del libro)
  var ESCENARIOS = ['dos', 'borde', 'escalon', 'impulso'];

  /** x[n] para el escenario dado. par: {f1, f2, offset}. n muestras a fs; el evento (escalón o impulso) está en T_EVENTO s. */
  function senal(esc, par, n, fs) {
    var x = new Float64Array(n), n0 = Math.round(T_EVENTO * fs), i, w1 = 2 * Math.PI * par.f1 / fs, w2 = 2 * Math.PI * par.f2 / fs;
    if (ESCENARIOS.indexOf(esc) < 0) throw new Error('escenario desconocido: ' + esc);
    for (i = 0; i < n; i++) {
      if (esc === 'dos') x[i] = A1 * Math.sin(w1 * i) + A2 * Math.sin(w2 * i);
      else if (esc === 'borde') x[i] = par.offset + A1 * Math.sin(w1 * i) + A2 * Math.sin(w2 * i);
      else if (esc === 'escalon') x[i] = i >= n0 ? A_ESCALON : 0;
      else x[i] = i === n0 ? A_IMPULSO : 0;
    }
    return x;
  }

  /** Salida en régimen permanente (ideal, registro infinito) de los escenarios «dos» y «borde» para `filt` (con fase) o `filtfilt`
   *  (fase cero, ganancia |H|²): Σ A·|H|ⁿ·sin(ω·i + φ) + offset·H(0), con H(0) = 1 en paso bajo y 0 en paso alto. */
  function estacionaria(d, par, n, fs, metodo, conOffset) {
    var y = new Float64Array(n), dos = metodo === 'filtfilt', comps = [[A1, par.f1], [A2, par.f2]], i;
    comps.forEach(function (c) {
      var h = respuestaSOS(d, c[1]), m = Math.hypot(h[0], h[1]), ph = Math.atan2(h[1], h[0]), w = 2 * Math.PI * c[1] / fs;
      var amp = c[0] * (dos ? m * m : m), f0 = dos ? 0 : ph;
      for (i = 0; i < n; i++) y[i] += amp * Math.sin(w * i + f0);
    });
    if (conOffset && d.tipo === 'lp') for (i = 0; i < n; i++) y[i] += par.offset;
    return y;
  }

  /** máx |a − b| en [i0, i1). */
  function maxDif(a, b, i0, i1) { var m = 0; for (var i = Math.max(0, i0); i < Math.min(a.length, i1); i++) m = Math.max(m, Math.abs(a[i] - b[i])); return m; }

  function maxAbs(y, i0, i1) { var m = 0; for (var i = Math.max(0, i0); i < Math.min(y.length, i1); i++) m = Math.max(m, Math.abs(y[i])); return m; }

  /** Índice de la última muestra con |y − yFinal| > tol, a partir de i0; devuelve el tiempo de asentamiento (s) desde i0,
   *  o null si no se asienta antes del final del registro. */
  function asentamiento(y, i0, yFinal, tol, fs) {
    var ult = -1;
    for (var i = i0; i < y.length; i++) if (Math.abs(y[i] - yFinal) > tol) ult = i;
    if (ult === y.length - 1) return null;
    return (ult < 0 ? 0 : ult + 1 - i0) / fs;
  }

  var api = {
    FS_SENAL: FS_SENAL, TIPOS: TIPOS, ESCENARIOS: ESCENARIOS, A1: A1, A2: A2, A_ESCALON: A_ESCALON, A_IMPULSO: A_IMPULSO, T_EVENTO: T_EVENTO,
    disenar: disenar, filt: filt, filtfilt: filtfilt, relleno: relleno, estadoEstacionario: estadoEstacionario,
    respuestaSOS: respuestaSOS, retardoGrupo: retardoGrupo, retardoFase: retardoFase, fase: fase, estable: estable,
    senal: senal, estacionaria: estacionaria, maxDif: maxDif, maxAbs: maxAbs, asentamiento: asentamiento
  };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.NMFiltfilt = api;
})(typeof window !== 'undefined' ? window : this);
