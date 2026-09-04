# Deploying the BLS demo

> Status: **live at <https://al-am.in/bls/>.**

## How it works

The demo is not a deployment of its own. `docs/` is vendored into the portfolio
repo (`eNipu/eNipu.github.io`) as `public/bls/`, and Cloudflare Workers Builds
deploys that repo on every push to `master`. So the demo ships with the site:
no separate Pages project, no DNS record, no API token, no repo secrets.

`docs/` is a self-contained static site (`index.html`, `style.css`, `app.js`,
`bls.js`, `bls.wasm`) that uses only relative URLs, so it works unchanged under
any path prefix. Astro copies `public/` to `dist/` verbatim; the portfolio's
`_headers` sets no CSP on `/bls/*`, which the Emscripten glue needs, and
Workers serves `.wasm` as `application/wasm`.

## Publishing a rebuilt demo

`.wasm` is committed, so there is no build step in CI. After rebuilding with
`examples/bls/build_wasm.sh` (needs emsdk + a WASM GMP):

```sh
rsync -a docs/ ../eNipu.github.io/public/bls/
cd ../eNipu.github.io && git commit -am "Update BLS demo" && git push
```

Workers Builds picks it up. `pairingma128` itself has no deploy workflow.

## Verification

1. `curl -sI https://al-am.in/bls/` → `HTTP 200`
2. `curl -sI https://al-am.in/bls/bls.wasm` → `content-type: application/wasm`
3. In a browser at `/bls/`: status reads "WebAssembly ready", then keygen →
   sign → verify (valid ✓, tampered ✗) → aggregate → verify (valid ✓,
   tampered ✗)
4. The portfolio's "Live demo" link (`src/data/projects.ts`, `/bls/`) resolves
   from both `/` and `/code`

## History

The original plan was a separate Cloudflare Pages project `bls-demo` behind
`bls.al-am.in`, driven by `.github/workflows/deploy-bls-demo.yml`. That workflow
failed on every run — the repo has no `CLOUDFLARE_API_TOKEN` — and the DNS record
was never created, so the URL never resolved. Riding along with the portfolio
removes the credentials, the subdomain and the workflow entirely; the workflow
and `deploy-cloudflare.sh` were deleted.
