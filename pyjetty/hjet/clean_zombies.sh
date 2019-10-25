#!/bin/bash

if [ "x${1}" == "xclean" ]; then
        echo "[i] clean..."
        # find . -name "*.root" -type f -size -1k -exec rm -v {} \;
        rfiles=$(find . -name "*.root" -type f)
        nfiles=$(find . -name "*.root" -type f | wc -l)
        for fn in $rfiles
        do
                echo "[i] file $fn"
                rval=$(python -c "import ROOT as r; f = r.TFile('${fn}'); print(f.IsZombie());")
                if [ "x$rval" == "xTrue" ]; then
                        mkdir -p /tmp/zombies
                        mv -v ${fn} /tmp/zombies
                fi
        done | tqdm --total $nfiles --unit "[i]" >> /dev/null
fi