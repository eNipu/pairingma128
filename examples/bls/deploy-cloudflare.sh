#!/usr/bin/env bash
# Deploy the BLS demo (docs/) to Cloudflare Pages as a separate project.
# One-time: `npx wrangler login`, then run this. Add a DNS CNAME
#   bls  ->  bls-demo.pages.dev
# in the al-am.in zone to serve it at https://bls.al-am.in/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
npx wrangler pages deploy docs --project-name=bls-demo
