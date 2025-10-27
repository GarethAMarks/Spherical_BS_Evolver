#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <output_file> <central_amplitude> <output_position>"
    exit 1
fi

file="$1"
A0="$2"
n="$3"

#extract nth line starting with A0. Note this uses string-matching so A0 must be exact (including trailing zeroes).
line=$(grep -m "$n"  "^$A0" "$file" | sed -n "${n}p")

# Check if a line was found
if [ -n "$line" ]; then
    echo "Found line in $file"
else
    echo "No line starting with '$A0' found in $file."
    exit 1
fi

#obtain needed outputs from output file
omegamin=$(echo -e "$line" | awk -v n=2 '{print $n}')
omegamax=$(echo -e "$line" | awk -v n=3 '{print $n}')
Phi0=$(echo -e "$line" | awk -v n=5 '{print $n}')
M=$(echo -e "$line" | awk -v n=6 '{print $n}')
zeta1=$(echo -e "$line" | awk -v n=7 '{print $n}')

echo "$M"

#get params file for Uli's code and store row labels
params="v3_thinshell/script_params.par"

A0begin="A0begin"
A0end="A0end"
Phi0s="Phi0"
Ms="mass0"
zeta1s="zeta1"
omegamins="omegamin"
omegamaxs="omegamax"

#replace params in params file with those read in
sed -i.bak "/$A0begin/ s/\($A0begin[[:space:]]*\)[^[:space:]]*/\1$A0/" "$params"
sed -i.bak "/$A0end/ s/\($A0end[[:space:]]*\)[^[:space:]]*/\1$A0/" "$params"
sed -i.bak "/$Phi0s/ s/\($Phi0s[[:space:]]*\)[^[:space:]]*/\1$Phi0/" "$params"
sed -i.bak "/$Ms/ s/\($Ms[[:space:]]*\)[^[:space:]]*/\1$M/" "$params"
sed -i.bak "/$zeta1s/ s/\($zeta1s[[:space:]]*\)[^[:space:]]*/\1$zeta1/" "$params"
sed -i.bak "/$omegamins/ s/\($omegamins[[:space:]]*\)[^[:space:]]*/\1$omegamin/" "$params"
sed -i.bak "/$omegamaxs/ s/\($omegamaxs[[:space:]]*\)[^[:space:]]*/\1$omegamax/" "$params"

#nano v3_thinshell/script_params.par

#copy data files 
./v3_thinshell/fshoot "$params"
#cp v3_thinshell/*.dat .
python3 smooth_and_interp.py
