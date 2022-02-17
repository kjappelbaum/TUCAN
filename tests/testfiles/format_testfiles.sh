#!/bin/bash

# Need to use quotation marks oround paths since they can contain spaces.

shopt -s extglob # enables use of patterns when stripping `moldata` to obtain `molname`
shopt -s nullglob # don't execute loop in case there's no matching `moldata`

for moldata in  ./*.mol ./*_metadata.txt ./*_drawing.*; do
    molname=${moldata%@(.mol|_metadata.txt|_drawing.*)}
    if [ ! -d "$molname" ]; then
        echo "Making ${molname} directory."
        mkdir "${molname}"
    fi
    mv "$moldata" "${molname}"
done

for moldir in */ ; do
    molname=${moldir%/}
    molpath=${moldir}${molname}.mol
    if [ ! -f "${molpath}.bak" ]; then
        echo "Converting ${molname} to v3000."
        mv "$molpath" "${molpath}.bak"
        obabel -imol "${molpath}.bak" -omol -O "$molpath" -x3
    fi
done
