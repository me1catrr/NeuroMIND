/* Núcleo numérico de las demostraciones docentes: respuesta EXACTA de un Butterworth digital y FFT.
 * Script clásico (no módulo): funciona abriendo el HTML desde disco (file://) y también desde Node.
 * Sin dependencias ni red.
 *
 * Fórmula. Para un diseño digital Butterworth de orden N construido como lo hace DSP.jl (prototipo analógico
 * + transformación de frecuencia + bilineal con prewarping), la respuesta compleja es
 *     H(e^{jω}) = Hp( s ),   Hp(s) = 1 / ∏_k (s − p_k),   p_k = exp[ jπ(1/2 + (2k−1)/(2N)) ],  k = 1..N
 * con Ω = tan(ω/2) y ω = 2πf/Fs, y
 *     paso bajo   (corte Ωc = tan(π fc/Fs)):               s = j Ω / Ωc
 *     paso alto:                                           s = −j Ωc / Ω
 *     banda eliminada (bordes Ω1, Ω2; B = Ω2−Ω1; Ω0² = Ω1·Ω2):   s = j B Ω / (Ω0² − Ω²)
 * La salida de `filtfilt` (doble pasada) tiene H_ef = H(z)·H(1/z) = |H|²: magnitud al cuadrado y fase cero.
 * Esto es la respuesta ideal de un registro infinito; NO incluye los efectos de borde de una señal finita.
 *
 * Validación: tests/verificar.mjs compara esta función con freqresp de DSP.jl (tests/fixtures/ref_freqresp.json).
 */
(function (root) {
  'use strict';

  var FS_VALIDAS = [250, 500, 1000];
  var N_MIN = 1, N_MAX = 8;

  function polos(N) {
    var p = [];
    for (var k = 1; k <= N; k++) {
      var a = Math.PI / 2 + (2 * k - 1) * Math.PI / (2 * N);
      p.push([Math.cos(a), Math.sin(a)]);
    }
    return p;
  }

  // Hp(s) con s = sre + j·sim
  function hProto(N, sre, sim) {
    var p = polos(N), dre = 1, dim = 0;
    for (var k = 0; k < N; k++) {
      var a = sre - p[k][0], b = sim - p[k][1];
      var t = dre * a - dim * b;
      dim = dre * b + dim * a;
      dre = t;
    }
    var m = dre * dre + dim * dim;
    return [dre / m, -dim / m];
  }

  /* cfg: {tipo:'hp'|'lp'|'bs', fs, N, fc}  o  {tipo:'bs', fs, N, f0, bw}.  hz: frecuencia en Hz (0 < hz <= fs/2). */
  function respuesta(cfg, hz) {
    var fs = cfg.fs, N = cfg.N, Om = Math.tan(Math.PI * hz / fs), sim, r;
    if (cfg.tipo === 'lp') {
      sim = Om / Math.tan(Math.PI * cfg.fc / fs);
    } else if (cfg.tipo === 'hp') {
      if (Om === 0) return [0, 0];
      sim = -Math.tan(Math.PI * cfg.fc / fs) / Om;
    } else if (cfg.tipo === 'bs') {
      var f1 = cfg.f0 - cfg.bw / 2, f2 = cfg.f0 + cfg.bw / 2;
      var O1 = Math.tan(Math.PI * f1 / fs), O2 = Math.tan(Math.PI * f2 / fs);
      var den = O1 * O2 - Om * Om;
      if (den === 0) return [0, 0];
      sim = (O2 - O1) * Om / den;
    } else {
      throw new Error('tipo desconocido: ' + cfg.tipo);
    }
    r = hProto(N, 0, sim);
    return r;
  }

  function mag(h) { return Math.sqrt(h[0] * h[0] + h[1] * h[1]); }
  function dB(x) { return 20 * Math.log10(Math.max(x, 1e-12)); }   // suelo de −240 dB: solo para dibujar

  /* Curvas para una rejilla de frecuencias: una pasada (filt) y doble pasada (filtfilt). */
  function curva(cfg, hzs) {
    var m1 = [], m2 = [], fase1 = [];
    for (var i = 0; i < hzs.length; i++) {
      var h = respuesta(cfg, hzs[i]), m = mag(h);
      m1.push(m);
      m2.push(m * m);
      fase1.push(Math.atan2(h[1], h[0]));
    }
    return { hz: hzs, mag1: m1, mag2: m2, fase1: fase1 };
  }

  /* Rejilla logarítmica entre fmin y fmax. */
  function rejillaLog(fmin, fmax, n) {
    var a = Math.log10(fmin), b = Math.log10(fmax), r = [];
    for (var i = 0; i < n; i++) r.push(Math.pow(10, a + (b - a) * i / (n - 1)));
    return r;
  }

  /* Límites válidos de la demostración. Devuelve {ok:boolean, mensaje:string}. */
  function validar(cfg) {
    if (FS_VALIDAS.indexOf(cfg.fs) < 0) return { ok: false, mensaje: 'Fs debe ser 250, 500 o 1000 Hz.' };
    if (!(cfg.N >= N_MIN && cfg.N <= N_MAX) || Math.floor(cfg.N) !== cfg.N)
      return { ok: false, mensaje: 'El orden de diseño debe ser un entero entre ' + N_MIN + ' y ' + N_MAX + '.' };
    var nyq = cfg.fs / 2;
    if (cfg.tipo === 'bs') {
      var f1 = cfg.f0 - cfg.bw / 2, f2 = cfg.f0 + cfg.bw / 2;
      if (!(cfg.bw > 0) || !(f1 > 0) || !(f2 < nyq))
        return { ok: false, mensaje: 'La banda eliminada debe cumplir 0 < f0−B/2 y f0+B/2 < Nyquist (' + nyq + ' Hz).' };
    } else if (!(cfg.fc > 0 && cfg.fc < nyq)) {
      return { ok: false, mensaje: 'La frecuencia de corte debe cumplir 0 < fc < Nyquist (' + nyq + ' Hz).' };
    }
    return { ok: true, mensaje: '' };
  }

  /* ───── FFT radix-2 (in situ) y espectro de amplitud con ventana de Hann ───── */
  function fft(re, im) {
    var n = re.length, i, j, k, m;
    for (i = 1, j = 0; i < n; i++) {
      var bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { var tr = re[i]; re[i] = re[j]; re[j] = tr; var ti = im[i]; im[i] = im[j]; im[j] = ti; }
    }
    for (var len = 2; len <= n; len <<= 1) {
      var ang = -2 * Math.PI / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (i = 0; i < n; i += len) {
        var cr = 1, ci = 0;
        for (k = 0; k < len / 2; k++) {
          m = i + k;
          var ur = re[m], ui = im[m];
          var vr = re[m + len / 2] * cr - im[m + len / 2] * ci;
          var vi = re[m + len / 2] * ci + im[m + len / 2] * cr;
          re[m] = ur + vr; im[m] = ui + vi;
          re[m + len / 2] = ur - vr; im[m + len / 2] = ui - vi;
          var nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
        }
      }
    }
  }

  /* A[k] = 2·|FFT(x·hann)[k]| / Σhann, k = 0..nfft/2 ; ventana de n muestras y relleno de ceros hasta nfft. */
  function espectro(x, nfft) {
    nfft = nfft || 4096;
    var n = x.length, re = new Float64Array(nfft), im = new Float64Array(nfft), sw = 0, i;
    for (i = 0; i < n; i++) {
      var w = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / (n - 1));
      re[i] = x[i] * w; sw += w;
    }
    fft(re, im);
    var A = new Float64Array(nfft / 2 + 1);
    for (i = 0; i <= nfft / 2; i++) A[i] = 2 * Math.sqrt(re[i] * re[i] + im[i] * im[i]) / sw;
    return A;
  }

  var api = {
    FS_VALIDAS: FS_VALIDAS, N_MIN: N_MIN, N_MAX: N_MAX,
    respuesta: respuesta, curva: curva, rejillaLog: rejillaLog, validar: validar, mag: mag, dB: dB,
    espectro: espectro
  };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.NMFiltros = api;
})(typeof window !== 'undefined' ? window : this);
