#!/usr/bin/env bash
# Deploy the packaged Propedia website tree to the web server — a PURE TRANSFER.
#
# Run `snakemake --cores N site --config mode=full` first: that leaves the finished tree
# under <web_dir>/data (default results/full/web/data), already in the site's layout
# (data/db/… + root csv/tsv/fasta/clusters + the small zips), thanks to package.web_mode_name.
#
# The destination is NEVER hardcoded here (this file is in a public repo). Provide it as
# $PROPEDIA_DEPLOY_DEST or the first argument, in rsync form and reachable from where you
# run this, e.g.:   user@host:/home/liase/www/propedia27-beta/public/data
# Nothing about the lab's servers lives in this file.
#
# Usage:
#   PROPEDIA_DEPLOY_DEST=user@host:/…/public/data  scripts/deploy_site.sh          # real
#   scripts/deploy_site.sh  user@host:/…/public/data                               # real
#   scripts/deploy_site.sh -n  user@host:/…/public/data  [web_data_dir]            # dry run
#
# Excludes the pep-pro complexes zip (propedia.zip is hosted on Zenodo, not the server)
# and the internal .packaged/.zipped/*.missing.txt markers.
set -euo pipefail

DRY=""
if [ "${1:-}" = "-n" ] || [ "${1:-}" = "--dry-run" ]; then DRY="-n"; shift; fi

DEST="${PROPEDIA_DEPLOY_DEST:-}"
if [ -z "${DEST}" ]; then          # not in the env -> take it from the first argument
  DEST="${1:-}"
  [ $# -gt 0 ] && shift
fi
WEB_DATA="${1:-results/full/web/data}"

if [ -z "${DEST}" ]; then
  echo "ERROR: no destination. Set PROPEDIA_DEPLOY_DEST or pass it as an argument," >&2
  echo "       e.g. user@host:/…/propedia27-beta/public/data" >&2
  exit 2
fi
if [ ! -d "${WEB_DATA}" ]; then
  echo "ERROR: web tree not found at '${WEB_DATA}'. Run 'snakemake site' first," >&2
  echo "       or pass the correct web_data_dir as the last argument." >&2
  exit 2
fi

echo "Deploying ${WEB_DATA}/  ->  ${DEST}/  ${DRY:+(dry run)}"
# -a archive, -z compress (helps the text csv/cif/fasta; zips are already compressed but
# rsync skips recompressing incompressible data), --partial resume, --omit-dir-times to
# avoid 'set times' errors on server-owned parent dirs. No --delete: never remove files
# the site added (index.html guards, projects/, etc.).
exec rsync -az --partial --omit-dir-times --info=progress2 ${DRY} \
  --exclude 'propedia.zip' \
  --exclude '.packaged' --exclude '.zipped' --exclude '*.missing.txt' \
  "${WEB_DATA}/" "${DEST}/"
