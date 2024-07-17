#!/bin/bash

set -e
set -u
set -o pipefail



TAXON=${1}
CLASS=${2}

      sed -i 's/ /_/g' ${TAXON}.missing_sp_to_delete.txt
      ids=($(cut -f1 MIDORI_${CLASS}.raw_id2acc.txt))
      for i in ${ids[@]}
      do
        sp=($(echo ${i} | cut -f1,2 -d'_'))
        if grep -q ${sp} ${TAXON}.missing_sp_to_delete.txt
        then
          printf '%s\t%s\n' ${i} "delete" >> MIDORI_${CLASS}.taxonomy_status.txt
        else
          printf '%s\t%s\n' ${i} "keep" >> MIDORI_${CLASS}.taxonomy_status.txt
          printf '%s\n' ${i} >> MIDORI_${CLASS}.seqs_to_retain.txt
        fi
      done
      #need to delete seqs by selecting ones to keep
      FILES=($(find . -mindepth 1 -maxdepth 1 -type f -name "MIDORI_*.mafft_edit.fa" | sed 's/\.\///g'))
      for file in ${FILES[@]}
      do
        label=($(echo ${file} | sed 's/.amp_blast.noN.mafft_edit.fa//g'))
        seqtk subseq ${file} MIDORI_${CLASS}.seqs_to_retain.txt > ${label}.${CLASS}.final.fa
      done
