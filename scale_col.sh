if [ $# -ne 3 ]; then
    echo "Usage: $0 file column factor"
    exit 1
fi

file="$1"
col="$2"
factor="$3"

awk -v c="$col" -v f="$factor" '
# If line starts with #, print it as-is
/^#/ { print; next }

# Otherwise, scale the requested column
{
    for (i=1; i<=NF; i++) {
        if (i==c) {
            $i = $i * f
        }
    }
    print
}' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
