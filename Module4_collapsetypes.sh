#!/bin/bash

set -e
set -u
set -o pipefail



TAXON=${1}
label=${2}

        touch ${label}.${TAXON}.uniq.fa
        splist=($(cut -f1 ${label}.${TAXON}.sp_list.txt))
        for sp in ${splist[@]}
        do
          grep ${sp} ${label}.${TAXON}.final_clean_relabel.unalign.fa | sed 's/>//g' | seqtk subseq ${label}.${TAXON}.final_clean_relabel.unalign.fa - > ${sp}.fa
          count=($(grep -c ">" ${sp}.fa))
          if [ $count -gt 1 ]
          then
            mafft --thread 4 --retree 2 --reorder ${sp}.fa > ${sp}.mafft.fa
            perl collapsetypes_v4.6.pl -infile=${sp}.mafft.fa -nrdiffs=0
            seqbuddy ${sp}.mafft.fa.unique_haplotypes --clean_seq > ${sp}.uniq.fa
                 #DY seqbuddy has a python error that adds "FTP Error: got more than 8192 bytes" to the end of ${sp}.uniq.fa
                 #DY add this line to remove any instances of "FTP Error: got more than 8192 bytes"
                 sed -i '/FTP Error: got more than 8192 bytes/d' ${sp}.uniq.fa
                 sed -i '/FTP Error: \[Errno 8] nodename nor servname provided, or not known/d'  ${sp}.uniq.fa
            cat ${label}.${TAXON}.uniq.fa ${sp}.uniq.fa > tmp
          else
            cat ${label}.${TAXON}.uniq.fa ${sp}.fa > tmp
          fi
          mv tmp ${label}.${TAXON}.uniq.fa
          rm ${sp}*
        done
