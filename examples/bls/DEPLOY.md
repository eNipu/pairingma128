# Deploying the BLS demo to `https://bls.al-am.in/` — plan

> Status: **ready to execute once credentials are provided.** All code and
> deployment scaffolding is already on `master`; only the Cloudflare-side
> actions remain, and they need account access.

## 1. Objective

Serve the BLS-signatures WebAssembly demo (the `docs/` static site in this repo)
at `https://bls.al-am.in/`, linked from the portfolio at `https://al-am.in/`.

## 2. What is already done

| item | state |
|------|-------|
| Refactor + BLS demo merged to `master` and pushed | ✅ |
| `docs/` static site (`index.html`, `app.js`, `style.css`, `bls.js`, `bls.wasm`) | ✅ |
| `examples/bls/deploy-cloudflare.sh` (one-line `wrangler pages deploy`) | ✅ |
| `.github/workflows/deploy-bls-demo.yml` (auto-deploy on push) | ✅ |
| Portfolio entry "BLS signatures + aggregation demo" → `https://bls.al-am.in/` | ✅ live |
| Cloudflare Pages project `bls-demo` created/deployed | ❌ |
| DNS `CNAME bls → bls-demo.pages.dev` | ❌ |

The only reason the URL is down is the two ❌ items above, which require
Cloudflare credentials.

## 3. Approach

Deploy the pre-built static site to a **separate Cloudflare Pages project**
(`bls-demo`), then attach the custom domain `bls.al-am.in`. No build step is
needed — the `.wasm` is committed, so Pages just serves `docs/` as-is.

This matches the portfolio's existing Cloudflare setup and keeps the demo
colocated with its source (no files copied into the portfolio repo).

## 4. Credentials required (one-time, revocable)

A **scoped Cloudflare API token** (dashboard → My Profile → API Tokens → Create):

- **Account → Cloudflare Pages → Edit** (create + deploy the project)
- **Zone → `al-am.in` → DNS → Edit** (create the CNAME)

Plus the **Account ID** (top of any dashboard URL, or Account Home sidebar).
The zone ID can be derived via the API, so it is not needed upfront.

## 5. Steps I will run (with the token)

```sh
# 0. use the provided credentials
export CLOUDFLARE_API_TOKEN=<token>
export CLOUDFLARE_ACCOUNT_ID=<account-id>

# 1. create + deploy the Pages project from docs/
npx wrangler pages deploy docs --project-name=bls-demo

# 2. attach the custom domain (auto-provisions the DNS CNAME when the zone
#    is in the same account; otherwise create it explicitly)
curl -X POST "https://api.cloudflare.com/client/v4/accounts/$CLOUDFLARE_ACCOUNT_ID/pages/projects/bls-demo/domains" \
  -H "Authorization: Bearer $CLOUDFLARE_API_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"name":"bls.al-am.in"}'

# 3. (fallback) explicit DNS record if step 2 did not create it
curl -X POST "https://api.cloudflare.com/client/v4/zones/$ZONE_ID/dns_records" \
  -H "Authorization: Bearer $CLOUDFLARE_API_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"type":"CNAME","name":"bls","content":"bls-demo.pages.dev","proxied":true}'
```

## 6. Verification checklist (I will confirm each)

1. `curl -I https://bls.al-am.in/` → `HTTP 200`
2. `.wasm` served with `Content-Type: application/wasm` (required for
   streaming instantiation; the Emscripten glue falls back to ArrayBuffer
   otherwise, so this is a nice-to-have, not a blocker)
3. Headless Chrome run: page shows `"WebAssembly ready"`, buttons enabled, and
   the full flow passes — keygen → sign → verify (valid ✓ / tampered ✗) →
   aggregate → verify (valid ✓ / tampered ✗)
4. Portfolio link resolves: `https://al-am.in/code` and the homepage both
   surface the demo entry, and clicking "Live demo" reaches the working site

## 7. Alternative (zero Cloudflare steps)

If a scoped token is not desired, the same static site can go live in ~30 s via
**GitHub Pages**: `pairingma128` → Settings → Pages → "Deploy from a branch" →
branch `master`, folder `/docs`. URL: `https://enipu.github.io/pairingma128/`.
The portfolio "Live demo" link would then point there instead.

## 8. Security & cleanup

- Token is **scoped** (Pages + the single zone's DNS) and **revocable**.
- After deployment succeeds, the token can be revoked; the GitHub Action in
  `.github/workflows/deploy-bls-demo.yml` remains available for future deploys
  using repo secrets (`CLOUDFLARE_API_TOKEN`, `CLOUDFLARE_ACCOUNT_ID`).
