// Runtime configuration for HeartGRN viewer.
// This file is loaded BEFORE viewer.js and exposes settings via window.HUB_CONFIG.
//
// Edit this file (then push to git) to change deployment-specific settings.
// All values are public — do NOT put server-side secrets here.

// Auto-detect the URL prefix so the viewer works both as
//   - https://178-105-162-190.nip.io/             → API_BASE ""
//   - https://www.heartcellatlas.org/cardionav/   → API_BASE "/cardionav"
//     (Sanger LB strips /cardionav/ before forwarding to our nginx,
//      so client→Sanger fetches use the prefix; our nginx sees root paths.)
(function () {
  var p = (window.location && window.location.pathname) || "/";
  window.__HUB_API_BASE__ = p.indexOf("/cardionav/") === 0 ? "/cardionav" : "";
})();

window.HUB_CONFIG = {
  // ---- Backend API URL ------------------------------------------------------
  // Same-origin: the browser-side prefix depends on the path the user reached.
  // Local dev override: set API_BASE: "http://127.0.0.1:8766", for direct backend.
  API_BASE: window.__HUB_API_BASE__,

  // ---- Page access gate -----------------------------------------------------
  // Server-side gate (nginx Basic Auth or Cloudflare Access) is the source of
  // truth. The client-side hash gate that used to live here is now disabled —
  // it was offline-brute-forceable from the public JS.
  REQUIRE_PASSWORD: false,
};
