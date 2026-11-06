#!/usr/bin/env bash
# Update the local ClinVar GRCh38 VCF from NCBI.
#
# Layout (in CLINVAR_DIR):
#   clinvar_YYYYMMDD/clinvar.vcf.gz (+ .tbi, .md5)   one folder per release
#   current -> clinvar_YYYYMMDD                        what the pipeline reads
# Keeps the current release plus one previous release; older ones are deleted.
# A new release is only switched to after its checksum, gzip stream and index check out,
# so `current` always points to a complete, matching VCF + index.
#
# Weekly cron example (Tuesday 03:00):
#   0 3 * * 2 /root/100plus_snpxplorer/definite_standalone/snpXplorer/Exploration/accessory/update_clinvar.sh >/dev/null 2>&1
set -euo pipefail

DIR="${CLINVAR_DIR:-/root/100plus_snpxplorer/definite_standalone/snpXplorer/Data/databases/clinvar}"
URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
LOG="$DIR/update_clinvar.log"

log() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

# release date (YYYYMMDD) from the NCBI .md5 file, which names e.g. .../clinvar_20261101.vcf.gz
release_of() { grep -oE 'clinvar_[0-9]{8}' "$1" | head -1 | cut -d_ -f2; }

# checksum, gzip stream and tabix index of the release in folder $1
verify() {
    local d="$1"
    [ "$(md5sum < "$d/clinvar.vcf.gz" | cut -d' ' -f1)" = "$(cut -d' ' -f1 "$d/clinvar.vcf.gz.md5")" ] || { log "checksum mismatch in $d"; return 1; }
    gzip -t "$d/clinvar.vcf.gz" || { log "corrupt gzip stream in $d"; return 1; }
    # `|| true`: head closes the pipe early, which pipefail would otherwise report as a failure
    local first magic
    first="$(zcat "$d/clinvar.vcf.gz" 2>/dev/null | head -1 || true)"
    [[ "$first" == "##fileformat=VCF"* ]] || { log "not a VCF in $d"; return 1; }
    magic="$(zcat "$d/clinvar.vcf.gz.tbi" 2>/dev/null | head -c 3 || true)"
    [ "$magic" = "TBI" ] || { log "bad tabix index in $d"; return 1; }
}

# point `current` at folder $1 in one atomic step (relative link, so it also works inside Docker)
switch_to() {
    ln -sfn "$1" "$DIR/current.new"
    mv -T "$DIR/current.new" "$DIR/current"
}

# keep the folder `current` points to plus the newest other release; delete the rest
prune() {
    local keep previous
    keep="$(readlink "$DIR/current")"
    previous="$(ls -d "$DIR"/clinvar_[0-9]* 2>/dev/null | xargs -r -n1 basename | grep -vx "$keep" | sort | tail -1 || true)"
    for d in "$DIR"/clinvar_[0-9]*; do
        [ -d "$d" ] || continue
        b="$(basename "$d")"
        if [ "$b" != "$keep" ] && [ "$b" != "$previous" ]; then
            rm -rf -- "$d" && log "removed old release $b"
        fi
    done
}

mkdir -p "$DIR"
cd "$DIR"

# one-time: adopt a plain download already sitting in DIR (files without the release layout)
if [ ! -L current ] && [ -f clinvar.vcf.gz ] && [ -f clinvar.vcf.gz.tbi ] && [ -f clinvar.vcf.gz.md5 ]; then
    rel="$(release_of clinvar.vcf.gz.md5)"
    mkdir -p "clinvar_$rel"
    mv clinvar.vcf.gz clinvar.vcf.gz.tbi clinvar.vcf.gz.md5 "clinvar_$rel/"
    if verify "clinvar_$rel"; then
        switch_to "clinvar_$rel"
        log "adopted existing download as release $rel"
    else
        # put the files back untouched so a failed check never destroys a download
        mv "clinvar_$rel"/clinvar.vcf.gz "clinvar_$rel"/clinvar.vcf.gz.tbi "clinvar_$rel"/clinvar.vcf.gz.md5 .
        rmdir "clinvar_$rel"
        log "existing download failed verification; left it in place and downloading fresh"
    fi
fi

tmp="$(mktemp -d "$DIR/.download.XXXXXX")"
trap 'rm -rf -- "$tmp"' EXIT

curl -sSf --retry 3 -o "$tmp/clinvar.vcf.gz.md5" "$URL/clinvar.vcf.gz.md5"
rel="$(release_of "$tmp/clinvar.vcf.gz.md5")"
[ -n "$rel" ] || { log "could not read the release date from NCBI's md5 file"; exit 1; }

if [ -L current ] && [ "$(readlink current)" = "clinvar_$rel" ]; then
    log "already up to date (release $rel)"
    prune
    exit 0
fi

log "downloading release $rel"
curl -sSf --retry 3 -o "$tmp/clinvar.vcf.gz" "$URL/clinvar.vcf.gz"
curl -sSf --retry 3 -o "$tmp/clinvar.vcf.gz.tbi" "$URL/clinvar.vcf.gz.tbi"
if ! verify "$tmp"; then
    log "release $rel failed verification; keeping $(readlink current 2>/dev/null || echo 'no release')"
    exit 1
fi

rm -rf -- "clinvar_$rel"
mv "$tmp" "clinvar_$rel"
trap - EXIT
chmod 755 "clinvar_$rel"
switch_to "clinvar_$rel"
log "switched to release $rel"
prune
