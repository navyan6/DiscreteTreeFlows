#!/bin/bash
# Merge the macOS keychain roots (which include this network's TLS-intercepting
# proxy CA) with certifi's bundle, so Python urllib/Biopython can reach NCBI.
# curl already trusts these; Python does not. Run once per machine.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p .certs
OUT=.certs/system_roots.pem
: > "$OUT"
security find-certificate -a -p /System/Library/Keychains/SystemRootCertificates.keychain >> "$OUT" 2>/dev/null || true
security find-certificate -a -p /Library/Keychains/System.keychain >> "$OUT" 2>/dev/null || true
python3 -c "import certifi; print(open(certifi.where()).read())" >> "$OUT"
echo "wrote $OUT ($(grep -c 'BEGIN CERTIFICATE' "$OUT") certificates)"
echo "Scripts pick this up automatically; for other tools: export SSL_CERT_FILE=$PWD/$OUT"
