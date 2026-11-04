/* Abre la ficha (<details class="nm-ficha" id="ficha-N">) a la que apunta el enlace o la dirección (#ficha-N).
 * Un <details> cerrado no se abre solo cuando se llega a él por un ancla. Script clásico, sin red ni dependencias. */
(function () {
  'use strict';
  function abrir() {
    var id = decodeURIComponent((window.location.hash || '').replace(/^#/, ''));
    if (!id) return;
    var el = document.getElementById(id);
    if (el && el.tagName === 'DETAILS' && !el.open) {
      el.open = true;
      el.scrollIntoView();
    }
  }
  window.addEventListener('hashchange', abrir);
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', abrir); else abrir();
})();
