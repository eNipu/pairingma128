/* app.js - binds the WebAssembly BLS module to the DOM. */
"use strict";

let wasmInit, wasmKeygen, wasmSign, wasmVerify, wasmAggregate, wasmAggVerify;

/* current keypair + signed messages (full hex kept in memory) */
let state = { sk: null, pk: null, sig: null, sigs: {}, msgs: {} };

const $ = (id) => document.getElementById(id);

/* The falsified report used when "tamper" is ticked: the LDL value is edited
   from 108 to 208 mg/dL before the report reaches the insurer. */
const FORGED_B =
  "PATIENT: Ana Rivera | TEST: LDL cholesterol | VALUE: 208 mg/dL | RESULT: ABNORMAL";

function trunc(hex, n = 20) {
  if (!hex) return "—";
  return hex.length <= n + 1 ? hex : hex.slice(0, n) + "…" + hex.slice(-8);
}

function badge(id, ok) {
  const el = $(id);
  el.textContent = ok ? "✓ valid" : "✗ invalid";
  el.className = "badge " + (ok ? "ok" : "fail");
}

function resetBadge(id, text) {
  const el = $(id);
  el.textContent = text;
  el.className = "badge idle";
}

/* Inline, next to the control that needs it. The demo used window.alert(),
   which blocks the tab and puts the message somewhere the eye is not. */
function hint(id, message) {
  $(id).textContent = message || "";
}

function setStatus(ready, message) {
  const el = $("status");
  el.textContent = message || (ready
    ? "webassembly ready — bls12 optimal-ate pairing, 128-bit"
    : "loading webassembly module");
  el.className = "status " + (ready ? "ready" : "loading");
}

function enableAll() {
  ["btn-keygen", "btn-sign", "btn-verify", "btn-signA", "btn-signB", "btn-agg"]
    .forEach((id) => { $(id).disabled = false; });
}

async function boot() {
  try {
    const Module = await createBlsModule();
    wasmInit       = Module.cwrap("wasm_bls_init", "number", []);
    wasmKeygen     = Module.cwrap("wasm_bls_keygen", "string", []);
    wasmSign       = Module.cwrap("wasm_bls_sign", "string", ["string", "string"]);
    wasmVerify     = Module.cwrap("wasm_bls_verify", "number", ["string", "string", "string"]);
    wasmAggregate  = Module.cwrap("wasm_bls_aggregate", "string", ["string"]);
    wasmAggVerify  = Module.cwrap("wasm_bls_aggregate_verify", "number",
                                  ["string", "string", "string"]);

    const t0 = performance.now();
    wasmInit();
    setStatus(true);
    console.log("BLS init in " + Math.round(performance.now() - t0) + " ms");
    enableAll();
  } catch (e) {
    setStatus(false, "webassembly failed to load: " + e);
  }
}

function clearSigned() {
  ["msg", "msgA", "msgB"].forEach((id) => $(id).removeAttribute("data-signed"));
}

function onKeygen() {
  const kp = wasmKeygen();
  const sep = kp.indexOf("|");
  state.sk = kp.slice(0, sep);
  state.pk = kp.slice(sep + 1);
  state.sig = null;
  state.sigs = {};
  state.msgs = {};
  clearSigned();
  $("out-sk").textContent = trunc(state.sk, 24);
  $("out-pk").textContent = trunc(state.pk, 48);
  $("out-sig").textContent = "not yet signed";
  $("out-agg").textContent = "not yet aggregated";
  resetBadge("out-verify", "awaiting a signature");
  resetBadge("out-aggverify", "awaiting an aggregate");
  hint("hint-sign", "");
  hint("hint-agg", "");
}

function onSign() {
  if (!state.sk) return hint("hint-sign", "Generate a key first.");
  hint("hint-sign", "");
  state.sig = wasmSign(state.sk, $("msg").value);
  $("msg").setAttribute("data-signed", "");
  $("out-sig").textContent = trunc(state.sig, 48);
  resetBadge("out-verify", "not yet verified");
}

function onVerify() {
  if (!state.sig) return hint("hint-sign", "Sign a report first.");
  hint("hint-sign", "");
  badge("out-verify", wasmVerify(state.pk, state.sig, $("msg").value) === 1);
}

function signSlot(which) {
  if (!state.sk) return hint("hint-agg", "Generate a key first.");
  hint("hint-agg", "");
  const input = $(which === "A" ? "msgA" : "msgB");
  state.sigs[which] = wasmSign(state.sk, input.value);
  state.msgs[which] = input.value;
  input.setAttribute("data-signed", "");
}

function onAggregate() {
  if (state.sigs.A == null || state.sigs.B == null) {
    return hint("hint-agg", "Sign both reports first.");
  }
  hint("hint-agg", "");
  const agg = wasmAggregate(state.sigs.A + ";" + state.sigs.B);
  $("out-agg").textContent = trunc(agg, 48);

  const msgB = $("tamper").checked ? FORGED_B : state.msgs.B;
  const ok = wasmAggVerify(
    agg,
    state.pk + ";" + state.pk,
    state.msgs.A + "\n" + msgB
  ) === 1;
  badge("out-aggverify", ok);
}

$("btn-keygen").addEventListener("click", onKeygen);
$("btn-sign").addEventListener("click", onSign);
$("btn-verify").addEventListener("click", onVerify);
$("btn-signA").addEventListener("click", () => signSlot("A"));
$("btn-signB").addEventListener("click", () => signSlot("B"));
$("btn-agg").addEventListener("click", onAggregate);

/* Editing a message after signing invalidates the verdict on screen, so it can
   never show a stale "valid" against text the user has since changed. */
$("msg").addEventListener("input", () => resetBadge("out-verify", "not yet verified"));
["msgA", "msgB"].forEach((id) => $(id).addEventListener("input", () => {
  $(id).removeAttribute("data-signed");
  const which = id === "msgA" ? "A" : "B";
  delete state.sigs[which];
  resetBadge("out-aggverify", "awaiting an aggregate");
}));

boot();
