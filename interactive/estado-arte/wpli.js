/* Núcleo numérico de la simulación PLI / wPLI (capítulo «Estado del arte»). Script clásico (file:// y Node), sin dependencias.
 *
 * Modelo sintético por época k, a una sola frecuencia (señales analíticas complejas):
 *     X_k = a_k·u_k + m·c_k + σ·nx_k
 *     Y_k = a_k·e^{−iφ}·u_k + m·c_k + σ·ny_k
 * con a_k ∈ {0,1} (época «acoplada» con probabilidad p) y u, c, nx, ny gaussianas complejas CN(0,1).
 *   u  = fuente que interactúa (su retardo entre regiones se ve como el desfase φ)
 *   c  = fuente común que ven los dos canales a la vez (mezcla instantánea, caricatura de la conducción de volumen): misma amplitud m
 *        y misma fase en X e Y; m = 0 → no hay fuente común
 *   n  = ruido independiente en cada canal
 * Espectro cruzado S_k = X_k·conj(Y_k). Para una época acoplada sin ruido ni mezcla, S_k = |u|²·e^{iφ}: Im S_k = |u|² sin φ.
 *
 * Métricas (definiciones explícitas; las mismas que ve quien usa la simulación):
 *   PLI    = | (1/N) Σ sign(Im S_k) |
 *   wPLI   = | Σ Im S_k | / Σ |Im S_k|                       (0 si el denominador es 0)
 *   dwPLI² = ((Σ Im S_k)² − Σ Im S_k²) / ((Σ |Im S_k|)² − Σ Im S_k²)   (0 si el denominador es 0; puede ser negativo)
 *   coherencia² = |Σ S_k|² / (Σ|X_k|² · Σ|Y_k|²)
 * Im S_k se pone a 0 cuando es solo redondeo (|Im S| ≤ 1e-12·|X||Y|); ver simular().
 * dwPLI² es el estimador «debiased» del wPLI² de Vinck et al. (2011), como _dwpli_moments en src/connectivity/wPLI.jl.
 * Diferencia con ese código: allí se suma eps() al denominador para evitar 0/0; aquí se devuelve 0 explícitamente.
 *
 * Generador: mulberry32 (aritmética de 32 bits) y Box–Muller; scripts/generar_simulacion_wpli.jl lo reimplementa idéntico
 * y tests/verificar.mjs compara ambos. Es geometría y estadística de datos SINTÉTICOS: no usa datos EEG. */
(function (root) {
  'use strict';

  /** Mulberry32: devuelve una función que da uniformes en [0, 1) a partir de una semilla entera de 32 bits. */
  function mulberry32(semilla) {
    var a = semilla >>> 0;
    return function () {
      a = (a + 0x6D2B79F5) >>> 0;
      var t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /** Gaussiana compleja CN(0,1) = (z1 + i·z2)/√2 con Box–Muller (consume dos uniformes). Devuelve [re, im]. */
  function gaussComp(rand) {
    var u1 = 1 - rand(), u2 = rand();                 // u1 ∈ (0, 1]: evita log(0)
    var r = Math.sqrt(-2 * Math.log(u1)), ang = 2 * Math.PI * u2;
    return [r * Math.cos(ang) * Math.SQRT1_2, r * Math.sin(ang) * Math.SQRT1_2];
  }

  var DEF = { N: 40, p: 1, fase: 60, sigma: 0.5, mezcla: 0, semilla: 1 };

  /** Simula N épocas. par.fase en grados. mezcla = m (amplitud de la fuente común vista por los dos canales; 0 = ninguna).
   *  Por época se consumen siempre 9 uniformes (1 para a_k y 8 para u, c, nx, ny), con cualquier parámetro. */
  function simular(par) {
    var P = Object.assign({}, DEF, par), N = P.N | 0;
    var rand = mulberry32(P.semilla), phi = P.fase * Math.PI / 180, cp = Math.cos(phi), sp = Math.sin(phi);
    var ms = P.mezcla;
    var out = { N: N, xr: new Float64Array(N), xi: new Float64Array(N), yr: new Float64Array(N), yi: new Float64Array(N),
                acoplada: new Uint8Array(N), imS: new Float64Array(N), reS: new Float64Array(N) };
    for (var k = 0; k < N; k++) {
      var a = rand() < P.p ? 1 : 0;
      var u = gaussComp(rand), c = gaussComp(rand), nx = gaussComp(rand), ny = gaussComp(rand);
      // e^{−iφ}·u = (cp − i·sp)(ur + i·ui) = (cp·ur + sp·ui) + i(cp·ui − sp·ur)
      var uyr = cp * u[0] + sp * u[1], uyi = cp * u[1] - sp * u[0];
      var xr = a * u[0] + ms * c[0] + P.sigma * nx[0], xi = a * u[1] + ms * c[1] + P.sigma * nx[1];
      var yr = a * uyr + ms * c[0] + P.sigma * ny[0], yi = a * uyi + ms * c[1] + P.sigma * ny[1];
      out.xr[k] = xr; out.xi[k] = xi; out.yr[k] = yr; out.yi[k] = yi; out.acoplada[k] = a;
      out.reS[k] = xr * yr + xi * yi;                 // Re(X·conj(Y))
      var im = xi * yr - xr * yi;                     // Im(X·conj(Y))
      // Cero numérico: si |Im S| no supera 1e-12·|X||Y|, es redondeo (p. ej. mezcla instantánea pura sin ruido, donde Im S = 0
      // exacto en teoría) y se toma 0; si no, wPLI y PLI medirían el signo de errores de 1e-17, que no tiene sentido.
      out.imS[k] = Math.abs(im) <= 1e-12 * Math.sqrt((xr * xr + xi * xi) * (yr * yr + yi * yi)) ? 0 : im;
    }
    return out;
  }

  /** Métricas a partir de las Im S_k (y, si se dan, de X e Y para la coherencia). */
  function metricas(imS, sim) {
    var N = imS.length, suma = 0, sumaAbs = 0, sumaCuad = 0, signos = 0;
    for (var k = 0; k < N; k++) {
      var v = imS[k];
      suma += v; sumaAbs += Math.abs(v); sumaCuad += v * v;
      signos += v > 0 ? 1 : (v < 0 ? -1 : 0);
    }
    var pli = N ? Math.abs(signos) / N : 0;
    var wpli = sumaAbs > 0 ? Math.abs(suma) / sumaAbs : 0;
    var den = sumaAbs * sumaAbs - sumaCuad;
    var dwpli2 = den > 0 ? (suma * suma - sumaCuad) / den : 0;
    var coh = null;
    if (sim) {
      var sr = 0, si = 0, px = 0, py = 0;
      for (var j = 0; j < N; j++) {
        sr += sim.reS[j]; si += sim.imS[j];
        px += sim.xr[j] * sim.xr[j] + sim.xi[j] * sim.xi[j];
        py += sim.yr[j] * sim.yr[j] + sim.yi[j] * sim.yi[j];
      }
      coh = px * py > 0 ? (sr * sr + si * si) / (px * py) : 0;
    }
    return { pli: pli, wpli: wpli, dwpli2: dwpli2, coherencia2: coh, N: N };
  }

  /** PLI y wPLI frente a φ (0–180°), misma semilla para todos los ángulos: la curva es la de UNA realización, no la esperada. */
  function curvaFase(par, pasos) {
    var out = [];
    for (var i = 0; i <= pasos; i++) {
      var fase = 180 * i / pasos, s = simular(Object.assign({}, par, { fase: fase })), m = metricas(s.imS, s);
      out.push({ fase: fase, pli: m.pli, wpli: m.wpli, coherencia2: m.coherencia2 });
    }
    return out;
  }

  var api = { mulberry32: mulberry32, gaussComp: gaussComp, simular: simular, metricas: metricas, curvaFase: curvaFase, DEF: DEF };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.NMWpli = api;
})(typeof window !== 'undefined' ? window : this);
