// Runtime configuration for HeartGRN viewer.
// This file is loaded BEFORE viewer.js and exposes settings via window.HUB_CONFIG.
//
// Edit this file (then push to git) to change deployment-specific settings.
// All values are public — do NOT put server-side secrets here.

window.HUB_CONFIG = {
  // ---- Backend API URL ------------------------------------------------------
  // Point this to the deployed Fly.io app for production.
  // Default (commented) covers local development.
  //
  // Production (post-deploy):
  //   API_BASE: "https://heartgrn-api.fly.dev",
  // Local dev:
  //   API_BASE: "http://127.0.0.1:8766",
  API_BASE: "https://heartgrn-api.fly.dev",

  // ---- Staging password gate ------------------------------------------------
  // Set REQUIRE_PASSWORD=true and PASSWORD_SHA256 to the SHA-256 of your password.
  // Compute on the command line:    echo -n 'your-password' | shasum -a 256
  //
  // Default password: "heartgrn-staging-2026"
  //   sha256("heartgrn-staging-2026") =
  //     a6eb31f7a98b6b9f81ae341c74b572ecff19fee07ef9c8e42add37ef59e4cf00
  //
  // To remove the password gate entirely, set REQUIRE_PASSWORD=false.
  REQUIRE_PASSWORD: true,
  PASSWORD_SHA256: "a6eb31f7a98b6b9f81ae341c74b572ecff19fee07ef9c8e42add37ef59e4cf00",
};
