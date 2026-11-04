/* Núcleo numérico de la simulación «del retardo de conducción al desfase». Script clásico (file:// y Node), sin dependencias.
 *
 * Una oscilación de frecuencia f que recorre un trayecto de longitud L con velocidad v llega con un retardo
 *     τ = L / v          (con L en mm y v en m/s, τ sale en ms: 1 mm / (1 m/s) = 1 ms)
 * y, respecto a la región de partida, con un desfase
 *     Δφ = 360° · f · τ   (f en Hz, τ en s).
 * Solo se puede observar el desfase envuelto en (−180°, 180°]: si τ supera medio periodo, el desfase observable
 * cambia de signo (un retraso parece un adelanto). Eso ocurre para v < v* = 2·f·L.
 * Es geometría de una sinusoide: no es un modelo de conducción nerviosa ni calcula wPLI. */
(function (root) {
  'use strict';
  function tauMs(Lmm, v) { return Lmm / v; }
  function faseDeg(fHz, tms) { return 360 * fHz * tms / 1000; }
  function envolver(deg) {                         // → (−180, 180]
    var w = deg - 360 * Math.floor((deg + 180) / 360);   // [−180, 180)
    return w === -180 ? 180 : w;
  }
  function velocidadCritica(fHz, Lmm) { return 2 * fHz * Lmm / 1000; }   // m/s: τ = T/2
  function periodoMs(fHz) { return 1000 / fHz; }
  var api = { tauMs: tauMs, faseDeg: faseDeg, envolver: envolver, velocidadCritica: velocidadCritica, periodoMs: periodoMs };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.NMRetardo = api;
})(typeof window !== 'undefined' ? window : this);
