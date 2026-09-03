/* app.js - binds the WebAssembly BLS module to the DOM. */
"use strict";

let wasmInit, wasmKeygen, wasmSign, wasmVerify, wasmAggregate, wasmAggVerify;

/* current keypair + signed messages (full hex kept in memory) */
let state = { sk: null, pk: null, sigs: [], msgs: [] };

const $ = (id) => document.getElementById(id);

function trunc(hex, n = 20) {
  if (!hex) return "—";
  return hex.length <= n + 1 ? hex : hex.slice(0, n) + "…" + hex.slice(-8);
}

function badge(id, ok) {
  const el = $(id);
  el.textContent = ok ? "✓ valid" : "✗ invalid";
  el.className = "badge " + (ok ? "ok" : "fail");
}

function setStatus(ready) {
  const el = $("status");
  el.textContent = ready ? "WebAssembly ready — BLS12 Optimal-Ate pairing (128-bit)" : "Loading…";
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
    setStatus(false);
    $("status").textContent = "Failed to load WebAssembly: " + e;
  }
}

function onKeygen() {
  const kp = wasmKeygen();
  const sep = kp.indexOf("|");
  state.sk = kp.slice(0, sep);
  state.pk = kp.slice(sep + 1);
  state.sigs = [];
  state.msgs = [];
  $("out-sk").textContent = trunc(state.sk, 24);
  $("out-pk").textContent = trunc(state.pk, 48);
  $("out-sig").textContent = "—";
  $("out-agg").textContent = "—";
  $("out-verify").className = "badge idle"; $("out-verify").textContent = "—";
  $("out-aggverify").className = "badge idle"; $("out-aggverify").textContent = "—";
}

function onSign() {
  if (!state.sk) return alert("Generate a key first.");
  state.sig = wasmSign(state.sk, $("msg").value);
  $("out-sig").textContent = trunc(state.sig, 48);
  $("out-verify").className = "badge idle"; $("out-verify").textContent = "—";
}

function onVerify() {
  if (!state.sig) return alert("Sign a message first.");
  const ok = wasmVerify(state.pk, state.sig, $("msg").value) === 1;
  badge("out-verify", ok);
}

function signSlot(which) {
  if (!state.sk) return alert("Generate a key first.");
  const input = $(which === "A" ? "msgA" : "msgB");
  const sig = wasmSign(state.sk, input.value);
  state.sigs[which] = sig;
  state.msgs[which] = input.value;
  input.style.outline = "1px solid var(--ok)";
}

function onAggregate() {
  if (state.sigs["A"] == null || state.sigs["B"] == null) {
    return alert("Sign both messages first.");
  }
  const agg = wasmAggregate(state.sigs["A"] + ";" + state.sigs["B"]);
  $("out-agg").textContent = trunc(agg, 48);

  const msgB = $("tamper").checked ? "tampered!" : state.msgs["B"];
  const ok = wasmAggVerify(
    agg,
    state.pk + ";" + state.pk,
    state.msgs["A"] + "\n" + msgB
  ) === 1;
  badge("out-aggverify", ok);
}

$("btn-keygen").addEventListener("click", onKeygen);
$("btn-sign").addEventListener("click", onSign);
$("btn-verify").addEventListener("click", onVerify);
$("btn-signA").addEventListener("click", () => signSlot("A"));
$("btn-signB").addEventListener("click", () => signSlot("B"));
$("btn-agg").addEventListener("click", onAggregate);

boot();
