#!/bin/bash

set -e
set -u
set -o pipefail



TAXON=${1}

      splist=($(cut -f1 MIDORI_${TAXON}.missing_sp.txt))
      for sp in ${splist[@]}
      do
        genus=($(echo ${sp} | cut -f1 -d'_'))
        if grep -q ${genus} <(cut -f6 ${TAXON}.combined_taxonomy.txt)
        then
          printf '%s\t%s\t%s\n' ${sp} "present" "genus_match" >> MIDORI.genus_status.txt
          tabtk grep -f7 "${genus}_" ${TAXON}.combined_taxonomy.txt | cut -f1-6 | sort -u | awk -v sp="${sp}" '{print $0 "\t" sp "\t" "genus_match" "\t" sp}' >> genus_match_taxo
        else
          if grep -q ${genus} <(cut -f9 ${TAXON}.combined_taxonomy.txt)
          then
            printf '%s\t%s\t%s\n' ${sp} "absent" "synonym_match" >> MIDORI.genus_status.txt
            species=($(echo ${sp} | cut -f2 -d'_'))
            tabtk grep -f9 "${genus}_" ${TAXON}.combined_taxonomy.txt | cut -f1-6 | sort -u | awk -v species="${species}" -v sp="${sp}" '{print $0 "\t" $6 "_" species "\t" "genus_synonym_match" "\t" sp}' >> genus_match_taxo
          else
            printf '%s\t%s\t%s\n' ${sp} "absent" "no_match" >> MIDORI.genus_status.txt
            printf '%s\n' ${sp} >> ${TAXON}.missing_sp_to_delete.txt
          fi
        fi
      done
