#!/bin/sh
cd "${0%/*}" || exit 1                               # Run from this directory

# Extract the beam tip streamwise displacement (x-component of W on the "right"
# patch) from every written time directory into a two-column data file
# (time, W_x) for plotting with tipDisplacement.gnuplot.
#
# Uses foamDictionary to read the field, so an OpenFOAM environment must be
# sourced. Run automatically at the end of Allrun, or by hand afterwards.

outFile="tipDisplacement.dat"
printf '# time    Wx\n' > "$outFile"

for Wfile in [0-9]*/beam_0/W; do
    [ -f "$Wfile" ] || continue
    timeDir=${Wfile%%/*}

    # Skip non-time directories (e.g. 0_orig, which also matches [0-9]*)
    case "$timeDir" in ''|*[!0-9.]*) continue;; esac

    # "value uniform (x y z)" -> x
    Wx=$(foamDictionary -entry boundaryField/right/value -value "$Wfile" 2>/dev/null \
            | tr -d '()' | awk '{print $2}')
    [ -n "$Wx" ] && printf '%s %s\n' "$timeDir" "$Wx" >> "$outFile"
done

# Sort data rows by time, keeping the header line on top
{ head -n 1 "$outFile"; tail -n +2 "$outFile" | sort -n -k1,1; } > "$outFile.tmp" \
    && mv "$outFile.tmp" "$outFile"

echo "Wrote tip displacement history to $outFile"
