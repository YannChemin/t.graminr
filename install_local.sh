#!/bin/sh
# Build t.graminr against the GRASS source tree and copy the result into
# every GRASS installation on this machine that carries addons.
#
# This exists because GRASS resolves a module from GRASS_ADDON_BASE
# (~/.grass8/addons) before the source tree's own dist directory. A stale
# copy there silently shadows a freshly built one, so a rebuild alone is
# not enough: `make` and this script go together.
#
# Usage: ./install_local.sh [GRASS source tree]   (default ~/dev/grass)
set -e

PGM=t.graminr
TOPDIR=${1:-$HOME/dev/grass}
DIST=$TOPDIR/dist.x86_64-pc-linux-gnu
CATEGORY=$HOME/dev/grass-addons/src/temporal/$PGM

make -C "$CATEGORY" MODULE_TOPDIR="$TOPDIR"

for base in "$HOME/.grass8/addons" /usr/local/grass86; do
    [ -d "$base/scripts" ] || continue
    if [ ! -w "$base/scripts" ]; then
        echo "$base is not writable; skipping (run as root to install there)"
        continue
    fi
    mkdir -p "$base/etc/$PGM/graminr_lib" "$base/etc/$PGM/params" \
             "$base/docs/html" "$base/docs/man/man1"
    install -m 755 "$DIST/scripts/$PGM"                   "$base/scripts/"
    install -m 644 "$DIST/etc/$PGM/graminr_lib/"*.py      "$base/etc/$PGM/graminr_lib/"
    install -m 644 "$DIST/etc/$PGM/params/"*.csv          "$base/etc/$PGM/params/"
    install -m 644 "$DIST/docs/html/$PGM.html"            "$base/docs/html/"
    install -m 644 "$DIST/docs/man/man1/$PGM.1"           "$base/docs/man/man1/"
    # The installed copy must match what was just built, byte for byte.
    cmp -s "$DIST/scripts/$PGM" "$base/scripts/$PGM" \
        || { echo "$base/scripts/$PGM does not match the build" >&2; exit 1; }
    for f in "$DIST/etc/$PGM/graminr_lib/"*.py; do
        cmp -s "$f" "$base/etc/$PGM/graminr_lib/$(basename "$f")" \
            || { echo "$base/etc/$PGM/graminr_lib/$(basename "$f") differs" >&2; exit 1; }
    done
    echo "installed and verified in $base"
done
