// Runtime configuration for HeartGRN viewer.
// This file is loaded BEFORE viewer.js and exposes settings via window.HUB_CONFIG.
//
// Edit this file (then push to git) to change deployment-specific settings.
// All values are public — do NOT put server-side secrets here.

window.HUB_CONFIG = {
  // ---- Backend API URL ------------------------------------------------------
  // Empty string = same-origin (browser hits the same host that served the
  // page, e.g. https://178-105-162-190.nip.io). Set to a full URL only if
  // backend is on a different host (cross-origin).
  //
  // Local dev:
  //   API_BASE: "http://127.0.0.1:8766",
  API_BASE: "",

  // ---- Page access gate -----------------------------------------------------
  // Server-side gate (nginx Basic Auth or Cloudflare Access) is the source of
  // truth. The client-side hash gate that used to live here is now disabled —
  // it was offline-brute-forceable from the public JS.
  REQUIRE_PASSWORD: false,
};
