#!/bin/bash

set -e
set -u
set -o pipefail

# Send STDOUT and STDERR to log file
exec > >(tee -a get_taxonomy.`date +%Y-%m-%d`.log)
exec 2> >(tee -a get_taxonomy.`date +%Y-%m-%d`.log >&2)

# run '. ~/.linuxify' # activates GNU versions of grep, sed, awk

# Run

# Steps and associated scripts:

# 1. Process twin-tagged metabarcoding data
#   - *read_preprocessing.sh*
# did not run because we have an alternative bioinformatic pipeline

# 2. Obtain initial taxonomic classification for target taxon
#   - *get_taxonomy.sh*
# The working folder is the folder with the scripts in it.
# usage: bash get_taxonomy.sh <taxonName> <taxonRank> <screenforbio>
# where:
# <taxonName> is the scientific name of the target taxon e.g. Vertebrata
# <taxonRank> is the classification rank of the target taxon e.g. superclass
# <screenforbio> is the path to the screenforbio-mbc directory
cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
bash ~/src/screenforbio-mbc-23GLG/get_taxonomy.sh Vertebrata subphylum ~/src/screenforbio-mbc-23GLG/
# Success.
# Fetching taxonomy of Vertebrata took ~144 hours
# Some of the genus fields have NA in them, even though the Genus species is present (e.g. NA Bufo arabicus), so i manually added in the missing genus names
# output is Vertebrata_ITIS_taxonomy.txt, a copy of which is in screenforbio-mbc-23GLG/archived_files/
mkdir archived_files
cp Vertebrata_ITIS_taxonomy.txt archived_files/Vertebrata_ITIS_taxonomy.txt

# 3. Generate non-redundant curated reference sequence database for target amplicon(s) and fix taxonomic classification
#   - *get_sequences.sh*
# usage: bash get_sequences.sh <extras> <gap-fill> <module> <taxon> <screenforbio>
# where:
# <extras> is 'yes' or 'no', indicating whether to add local FASTA format sequences. if 'yes', files must be in present directory labelled "extra_12S.fa", "extra_16S.fa", "extra_Cytb.fa", "extra_COI.fa", with headers in format Genus_species_uniqueID.
# <gap-fill> is 'no' or a tab-delimited text file of species names to be targeted for gap-filling from NCBI, in format Genus_species.
# <module> is 'one', 'two', 'three' or 'four' indicating whether the script is starting from scratch ('one'), restarting after checking the output of the mafft alignment ('two'), restarting after manual correction of failed taxonomy lookups ('three'), or restarting after manual checks of SATIVA output ('four'). see end of module messages for any requirements for the next module."
# <taxon> is the taxon for which the taxonomy was downloaded with get_taxonomy.sh, e.g. Mammalia or Vertebrata (all outputs should be in present directory).
# <screenforbio> is the path to the screenforbio-mbc directory

# Module 0 - OPTIONAL Create extra_12S.fa and extra_16S.fa files from the Gaoligongshan species' mitogenomes
	# As we only have access to the mitochondrial genomes of these Gaoligongshan species, but no publication rights, the PROTAX reference database does not contain these unpublished sequences.


# Module 1 - Extract subset of raw Midori database for query taxon and loci. Remove sequences with non-binomial species names, reduce subspecies to species labels. Add local sequences (optional). Check for relevant new sequences for list of query species on NCBI (GenBank and RefSeq) (optional). Select amplicon region and remove primers. Remove sequences with ambiguous bases. Align. End of module: optional check of alignments
cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
# copy MIDORI files that i want to process to the working directory: screenforbio-mbc-23GLG/
cp archived_files/MIDORI2_UNIQ_NUC_GB253_lrRNA_RDP.fasta.gz ./MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta.gz; gunzip MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta.gz
cp archived_files/MIDORI2_UNIQ_NUC_GB253_srRNA_RDP.fasta.gz ./MIDORI_UNIQUE_1.4_srRNA_RDP.fasta.gz; gunzip MIDORI_UNIQUE_1.4_srRNA_RDP.fasta.gz

# The new versions of MIDORI have more complex headers, which interfere with the `get_sequences.sh` code.
# V 1.1:  `>AF382008	root;Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo;Homo sapiens`
# V 1.4:  `>AF382008.3.1672.3229	root_1;Eukaryota_2759;Chordata_7711;Mammalia_40674;Primates_9443;Hominidae_9604;Homo_9605;Homo sapiens_9606`
sed -i 's/\.[0-9].*\t/\t/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/ sp__/sp/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/ sp_/sp/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/phylum_//g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/class_//g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/order_//g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/family_//g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/genus_//g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/_[0-9]\{1,9\};/;/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/_[0-9]\{1,9\}$//g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/_;/;/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/_$//g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/_ /-/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/ _/-/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/_/-/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/sp[0-9]\{4,7\} .*/ sp/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40; sed -i 's/sp[0-9]\{1,9\}$/ sp/g' MIDORI_UNIQUE_1.4_lrRNA_RDP.fasta | head -n 40
sed -i 's/\.[0-9].*\t/\t/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/ sp__/sp/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/ sp_/sp/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/phylum_//g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/class_//g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/order_//g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/family_//g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/genus_//g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/_[0-9]\{1,9\};/;/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/_[0-9]\{1,9\}$//g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/_;/;/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/_$//g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/_ /-/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/ _/-/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/_/-/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/sp[0-9]\{4,7\} .*/ sp/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40; sed -i 's/sp[0-9]\{1,9\}$/ sp/g' MIDORI_UNIQUE_1.4_srRNA_RDP.fasta | head -n 40

# there are three 12S_primers files:  12S_primers_kocher.fa, 12S_primers_riaz.fa and 12S_primers_mifish.fa. Duplicate the one that you want to use and change its filename to 12S_primers.fa, and this will be the pair used to pull out 12S amplicons from the Midori reference database
# to use the Kocher primers (12S_primers_kocher.fa), in get_sequences.sh:
     # change the usearch -search_pcr line 178 to
     # usearch -search_pcr ${label}.raw.fa -db ${SCRIPTS}/12S_primers.fa -strand both -maxdiffs 4 -minamp 420 -maxamp 470 -ampout ${label}.amp.fa
     # change the usearch -fastx_truncate line 183 to
     # usearch -fastx_truncate ${label}.amp.fa -stripleft 30 -stripright 28 -fastaout ${label}.amp_only.fa # for Kocher 12S primers
     # change the awk line 188 to
     # cat ${label}.amp.blastn | awk 'BEGIN{FS=OFS}($4>=360){print $1 OFS $7 OFS $8}' > ${label}.amp.blastn.coords # for 12S Kocher primers
# to use the Riaz primers (12S_primers_riaz.fa), in get_sequences.sh:
     # change the usearch -search_pcr line 179 to
     # usearch11 -search_pcr2 ${label}.raw.fa -fwdprimer ACTGGGATTAGATACCCC -revprimer YRGAACAGGCTCCTCTAG -minamp 80 -maxamp 120 -strand both -maxdiffs 4 -fastaout ${label}.amp.fa # for 12S Riaz primers
     # change the usearch -fastx_truncate line 184 to
     # usearch -fastx_truncate ${label}.amp.fa -stripleft 0 -stripright 0 -fastaout ${label}.amp_only.fa # for Riaz 12S primers
     # change the awk line 189 to
     # cat ${label}.amp.blastn | awk 'BEGIN{FS=OFS}($4>=80){print $1 OFS $7 OFS $8}' > ${label}.amp.blastn.coords # for 12S Riaz primers

bash ~/src/screenforbio-mbc-23GLG/get_sequences.sh yes no one Vertebrata ~/src/screenforbio-mbc-23GLG/
# the first no is changed from no to yes to add the extra sequences, which are in the files extra_12S.fa and extra_16S.fa.
# Successful
# Module 1 took 6.98 hours (16Smam and 12SRiaz primers, Midori 1.4, Vertebrata)

# Actions after Module 1 complete
# Module 1 complete. Stopping now for manual inspection of alignments *.mafft.fa inside ./intermediate_files.
# Duplicates have been identified in the "MIDORI_*.amp_blast.noN.mafft.fa" file, do deduplication for data consistency.
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' intermediate_files/MIDORI_srRNA.amp_blast.noN.mafft.fa > MIDORI_srRNA.amp_blast.noN.mafft.singleline.fasta
cat MIDORI_srRNA.amp_blast.noN.mafft.singleline.fasta | tr '\n' '#' | sed 's/#>/\n>/g' | sort | uniq | sed 's/#/\n/g' > MIDORI_srRNA.amp_blast.noN.mafft_edit.fa
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' intermediate_files/MIDORI_lrRNA.amp_blast.noN.mafft.fa > MIDORI_lrRNA.amp_blast.noN.mafft.singleline.fasta
cat MIDORI_lrRNA.amp_blast.noN.mafft.singleline.fasta | tr '\n' '#' | sed 's/#>/\n>/g' | sort | uniq | sed 's/#/\n/g' > MIDORI_lrRNA.amp_blast.noN.mafft_edit.fa
rm -f MIDORI_*.amp_blast.noN.mafft.singleline.fasta
# The reverse matching process alters sequence names by adding the prefix "_R_"; this needs to be corrected.
cat MIDORI_srRNA.amp_blast.noN.mafft_edit.fa | sed 's/^>_R_/>/g' | tr '\n' '#' | sed 's/#>/\n>/g' | sort | sed 's/#/\n/g' > TMP; mv TMP MIDORI_srRNA.amp_blast.noN.mafft_edit.fa
cat MIDORI_lrRNA.amp_blast.noN.mafft_edit.fa | sed 's/^>_R_/>/g' | tr '\n' '#' | sed 's/#>/\n>/g' | sort | sed 's/#/\n/g' > TMP; mv TMP MIDORI_lrRNA.amp_blast.noN.mafft_edit.fa
# Manually check for multiple sequences in the same Accession and remove those that are clearly the result of primer mismatches.
# Haplotype de-redundancy of sequences from the same species
bash collapse_haplotypes_afterModule3.sh srRNA ~/src/screenforbio-mbc-23GLG/
bash collapse_haplotypes_afterModule3.sh lrRNA ~/src/screenforbio-mbc-23GLG/
# MAFFT
mafft --adjustdirection --retree 2 --reorder --thread 4 MIDORI_srRNA.amp_blast.noN.uniq.fasta > MIDORI_srRNA.amp_blast.noN.mafft.fa
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' MIDORI_srRNA.amp_blast.noN.mafft.fa > MIDORI_srRNA.amp_blast.noN.mafft.singleline.fasta
cat MIDORI_srRNA.amp_blast.noN.mafft.singleline.fasta | tr '\n' '#' | sed 's/#>/\n>/g' | sort | sed 's/#/\n/g' > MIDORI_srRNA.amp_blast.noN.mafft_edit.fa
mafft --adjustdirection --retree 2 --reorder --thread 4 MIDORI_lrRNA.amp_blast.noN.mafft_edit.fasta > MIDORI_lrRNA.amp_blast.noN.mafft.fa
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' MIDORI_lrRNA.amp_blast.noN.mafft.fa > MIDORI_lrRNA.amp_blast.noN.mafft.singleline.fasta
cat MIDORI_lrRNA.amp_blast.noN.mafft.singleline.fasta | tr '\n' '#' | sed 's/#>/\n>/g' | sort | sed 's/#/\n/g' > MIDORI_lrRNA.amp_blast.noN.mafft_edit.fa
rm -f MIDORI_*.amp_blast.noN.mafft.singleline.fasta; rm -f MIDORI_*.amp_blast.noN.mafft.fa; rm -f MIDORI_*.amp_blast.noN.mafft_edit.fasta


# Module 2 - Compare sequence species labels in the MIDORI fasta files with the ITIS taxonomy. Non-matching labels are queried against Catalogue of Life to check for known synonyms. Remaining mismatches kept if genus already exists in taxonomy, otherwise flagged for removal. End of module: optional check of flagged species labels.
# Requires a taxon_ITIS_taxonomy.txt file (e.g. Vertebrata_ITIS_taxonomy.txt file)
cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
bash ~/src/screenforbio-mbc-23GLG/get_sequences.sh yes no two Vertebrata ~/src/screenforbio-mbc-23GLG/
# Some species cause the taxize::classification() function inside get_sequences.sh to fail, throwing up the following error:
     # Retrieving data for taxon 'Zygogeomys trichopus'
     #
     # Error in doc_parse_raw(x, encoding = encoding, base_url = base_url, as_html = as_html,  :
     #   CData section not finished
     #  [63]
     # Calls: classification ... read_xml.character -> read_xml.raw -> doc_parse_raw
     # Execution halted
# The developer of taxize has flagged this bug in the Catalogue of Life, to be dealt with in the next release of taxize (for 0.9.8). UPDATE:  the problem appears to be in read_xml2, which is a separate module, which cannot read some non-standard ASCII characters in the names of some taxonomists. UPDATE: The 0.9.8 version of taxize has fixed the problem. The fix below is no longer needed.
     # To ID the species causing this error, open `classification_misbehavers_finder.Rmd`. This script loads "MIDORI_taxon.ITIS_mismatch_sp.txt", uses a tryCatch() loop to run taxize::classification() on each name in this file, and records the species that causes taxize::classification() to throw crashing errors.
     # We call these crashing species 'misbehavers,' and we remove them manually starting at line 532, which is before get_taxonomy_mismatches.R is run
          # remove misbehavers by deleting these species, like this:
          # sed -i '/Hemidactylus adensis/d' MIDORI_${TAXON}.ITIS_mismatch_sp.txt
     # Then add back the removed species starting at line 565, like this (WITH UNDERSCORE!!):
          # sed -i '$a\Hemidactylus_adensis\' MIDORI_${TAXON}.missing_sp.txt
# Success
# There are two output files:
    # Vertebrata.missing_sp_to_delete.txt
    # Vertebrata.combined_taxonomy.txt
    # Check for “NA” in Vertebrata.combined_taxonomy.txt and fix it manually
# Actions after Module 2 complete
# Module 2 complete. Stopping now for manual inspection of failed species lookups (in Vertebrata.missing_sp_to_delete.txt).
# If a failed lookup can be resolved, remove from Vertebrata.missing_sp_to_delete.txt and add taxonomy to a tab-delimited file named Vertebrata.missing_sp_taxonomy.txt with columns for kingdom,phylum,class,order,family,genus,species,status,query - 'status' should be something short and descriptive (_ instead of spaces; eg. mispelling or manual_synonym) and 'query' should be the entry in Vertebrata.missing_sp_to_delete.txt. Vertebrata.missing_sp_taxonomy.txt must not have a header line when the script is restarted.
# If all failed lookups are resolved, delete Vertebrata.missing_sp_to_delete.txt. If some/all failed lookups cannot be resolved, keep the relevant species names in Vertebrata.missing_sp_to_delete.txt. When restarting the script it will check for the presence of this file and act accordingly (sequences for these species will be discarded).
# If no failed lookups can be resolved, do not create Vertebrata.missing_sp_taxonomy.txt, leave Vertebrata.missing_sp_to_delete.txt as it is.
# Restart script when happy.
# Check if the species in the target area are all in Vertebrata.combined_taxonomy.txt, and manually add those that are not in Vertebrata.combined_taxonomy.txt
# Copies of these two files are in the archived/ folder after the run
	cp Vertebrata.missing_sp_to_delete.txt archived_files/
	cp Vertebrata.combined_taxonomy.txt archived_files/


# Module 3 - Discard flagged sequences. Update taxonomy key file for sequences found to be incorrectly labelled in Module 2. Run SATIVA. End of module: optional check of putatively mislabelled sequences
# requires file Vertebrata.combined_taxonomy.txt from Module 2, or there is a copy in archived_files/
cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
bash ~/src/screenforbio-mbc-23GLG/get_sequences.sh yes no three Vertebrata ~/src/screenforbio-mbc-23GLG/
# failed because of unknown reason
# Separate treatment of terrestrial vertebrates and fish
# Tetropoda
grep -E "Amphibia|Aves|Mammalia|Reptilia" Vertebrata.combined_taxonomy.txt > Tetropoda.combined_taxonomy.txt
tabtk isct -1 9 -2 2 Tetropoda.combined_taxonomy.txt MIDORI.raw_id2sp.txt > MIDORI_Tetropoda.raw_id2sp.txt
tabtk isct -1 1 -2 1 MIDORI_Tetropoda.raw_id2sp.txt MIDORI.raw_id2acc.txt > MIDORI_Tetropoda.raw_id2acc.txt
# discard failed lookups
bash Module3_discard_failed_lookups.sh Vertebrata Tetropoda
# get protax & sativa names
join -1 2 -2 9 -o 1.1,1.3,2.7 <(cut -f1 MIDORI_Tetropoda.raw_id2sp.txt | awk 'BEGIN{FS="_"}{print $1 "_" $2 "_" $3 "\t" $1 "_" $2 "\t" $3}' | sort -k2,2) <(sort -k9,9 Tetropoda.combined_taxonomy.txt) > MIDORI_Tetropoda.final_id2acc2sp.txt
awk '{print $1 "\t" $2}' MIDORI_Tetropoda.final_id2acc2sp.txt > MIDORI_Tetropoda.final_rename_seqs_sativa.txt
# Generating master sativa taxonomy file
join -1 3 -2 9 -o 1.2,2.2,2.3,2.4,2.5,2.6,2.9 <(sort -k3,3 MIDORI_Tetropoda.final_id2acc2sp.txt) <(sort -k9,9 Tetropoda.combined_taxonomy.txt) | awk '{print $1 "\t" "Eukaryota;" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7}' | sed 's/_/ /g' | sort -u > Tetropoda.final_taxonomy_sativa.txt
# Working on 12S
#rename alignments for sativa
seqkit replace -p '(.+)$' -r '{kv}' -k MIDORI_Tetropoda.final_rename_seqs_sativa.txt MIDORI_srRNA.Tetropoda.final.fa --keep-key > MIDORI_srRNA.Tetropoda.final_for_sativa.fa
#make sativa taxonomy
join -1 1 -2 1 <(grep ">" MIDORI_srRNA.Tetropoda.final_for_sativa.fa | sed 's/>//g' | sort -k1,1) <(sort -k1,1 Tetropoda.final_taxonomy_sativa.txt) | awk '{print $1 "\t" $2 " " $3}' > MIDORI_srRNA.Tetropoda.final_for_sativa.tax
# sativa
mkdir MIDORI_srRNA.Tetropoda_sativa
python2 ~/src/sativa/sativa.py -s MIDORI_srRNA.Tetropoda.final_for_sativa.fa -t MIDORI_srRNA.Tetropoda.final_for_sativa.tax -x ZOO -T 4 -n MIDORI_srRNA.Tetropoda -o MIDORI_srRNA.Tetropoda_sativa
# Working on 16S
#rename alignments for sativa
seqkit replace -p '(.+)$' -r '{kv}' -k MIDORI_Tetropoda.final_rename_seqs_sativa.txt MIDORI_lrRNA.Tetropoda.final.fa --keep-key > MIDORI_lrRNA.Tetropoda.final_for_sativa.fa
#make sativa taxonomy
join -1 1 -2 1 <(grep ">" MIDORI_lrRNA.Tetropoda.final_for_sativa.fa | sed 's/>//g' | sort -k1,1) <(sort -k1,1 Tetropoda.final_taxonomy_sativa.txt) | awk '{print $1 "\t" $2 " " $3}' > MIDORI_lrRNA.Tetropoda.final_for_sativa.tax
# sativa
mkdir MIDORI_lrRNA.Tetropoda_sativa
python2 ~/src/sativa/sativa.py -s MIDORI_lrRNA.Tetropoda.final_for_sativa.fa -t MIDORI_lrRNA.Tetropoda.final_for_sativa.tax -x ZOO -T 4 -n MIDORI_lrRNA.Tetropoda -o MIDORI_lrRNA.Tetropoda_sativa
# Teleostei
grep -e $'\t'"Teleostei"$'\t' Vertebrata.combined_taxonomy.txt > Teleostei.combined_taxonomy.txt
tabtk isct -1 9 -2 2 Teleostei.combined_taxonomy.txt MIDORI.raw_id2sp.txt > MIDORI_Teleostei.raw_id2sp.txt
tabtk isct -1 1 -2 1 MIDORI_Teleostei.raw_id2sp.txt MIDORI.raw_id2acc.txt > MIDORI_Teleostei.raw_id2acc.txt
# discard failed lookups
bash Module3_discard_failed_lookups.sh Vertebrata Teleostei
# get protax & sativa names
join -1 2 -2 9 -o 1.1,1.3,2.7 <(cut -f1 MIDORI_Teleostei.raw_id2sp.txt | awk 'BEGIN{FS="_"}{print $1 "_" $2 "_" $3 "\t" $1 "_" $2 "\t" $3}' | sort -k2,2) <(sort -k9,9 Teleostei.combined_taxonomy.txt) > MIDORI_Teleostei.final_id2acc2sp.txt
awk '{print $1 "\t" $2}' MIDORI_Teleostei.final_id2acc2sp.txt > MIDORI_Teleostei.final_rename_seqs_sativa.txt
# Generating master sativa taxonomy file
join -1 3 -2 9 -o 1.2,2.2,2.3,2.4,2.5,2.6,2.9 <(sort -k3,3 MIDORI_Teleostei.final_id2acc2sp.txt) <(sort -k9,9 Teleostei.combined_taxonomy.txt) | awk '{print $1 "\t" "Eukaryota;" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7}' | sed 's/_/ /g' | sort -u > Teleostei.final_taxonomy_sativa.txt
# Working on 12S
#rename alignments for sativa
seqkit replace -p '(.+)$' -r '{kv}' -k MIDORI_Teleostei.final_rename_seqs_sativa.txt MIDORI_srRNA.Teleostei.final.fa --keep-key > MIDORI_srRNA.Teleostei.final_for_sativa.fa
#make sativa taxonomy
join -1 1 -2 1 <(grep ">" MIDORI_srRNA.Teleostei.final_for_sativa.fa | sed 's/>//g' | sort -k1,1) <(sort -k1,1 Teleostei.final_taxonomy_sativa.txt) | awk '{print $1 "\t" $2 " " $3}' > MIDORI_srRNA.Teleostei.final_for_sativa.tax
# sativa
mkdir MIDORI_srRNA.Teleostei_sativa
python2 ~/src/sativa/sativa.py -s MIDORI_srRNA.Teleostei.final_for_sativa.fa -t MIDORI_srRNA.Teleostei.final_for_sativa.tax -x ZOO -T 4 -n MIDORI_srRNA.Teleostei -o MIDORI_srRNA.Teleostei_sativa
# SATIVA crashed, it was found that the species name "Lacépède" contains characters outside the ASCII encoding range, so the corresponding sequences "JX867257", "KC768227", "KF356196" were manually deleted from MIDORI_srRNA.Teleostei.final_for_sativa.fa, and then restarted from line197.
# Working on 16S
#rename alignments for sativa
seqkit replace -p '(.+)$' -r '{kv}' -k MIDORI_Teleostei.final_rename_seqs_sativa.txt MIDORI_lrRNA.Teleostei.final.fa --keep-key > MIDORI_lrRNA.Teleostei.final_for_sativa.fa
#make sativa taxonomy
join -1 1 -2 1 <(grep ">" MIDORI_lrRNA.Teleostei.final_for_sativa.fa | sed 's/>//g' | sort -k1,1) <(sort -k1,1 Teleostei.final_taxonomy_sativa.txt) | awk '{print $1 "\t" $2 " " $3}' > MIDORI_lrRNA.Teleostei.final_for_sativa.tax
# check if there's "Lacépède" in MIDORI_lrRNA.Teleostei.final_for_sativa.tax, manually deleted the corresponding sequences from MIDORI_lrRNA.Teleostei.final_for_sativa.fa and rerun from line206
# sativa
mkdir MIDORI_lrRNA.Teleostei_sativa
python2 ~/src/sativa/sativa.py -s MIDORI_lrRNA.Teleostei.final_for_sativa.fa -t MIDORI_lrRNA.Teleostei.final_for_sativa.tax -x ZOO -T 4 -n MIDORI_lrRNA.Teleostei -o MIDORI_lrRNA.Teleostei_sativa
# success
# Actions after Module 3 complete
# Module 3 complete. Stopping now for manual inspection of mislabelled sequences in ./MIDORI_locus_sativa/MIDORI_locus.mis
# For sequences where species-level or genus-level mislabelling can be resolved, make corrections directly in Vertebrata.final_taxonomy_sativa.txt (i.e. replace the taxonomic classification for that sequence with the correct one), this will be used to rename sequences.
# Make higher level changes to the taxonomy at your own risk - untested.
# Restart script when happy.
# archive an original version before changing it with the sativa suggestions
cp Tetropoda.final_taxonomy_sativa.txt ./intermediate_files/Tetropoda.final_taxonomy_sativa_orig.txt
cp Teleostei.final_taxonomy_sativa.txt ./intermediate_files/Teleostei.final_taxonomy_sativa_orig.txt
# then use this R code (delete_seqs_suggested_by_sativa.Rmd) to remove sequences that sativa identifies as incorrect
     # I remove all sequences that sativa identifies as having an incorrect taxonomy at family level and above, as such large errors are most likely to be database errors.
     # I ignore sativa's proposed substitute taxonomies below family level (mostly genus and species), because some unknown, possibly large, proportion of these are due to the *tree* being incorrect, since these are short seqs
#cleanup
mv MIDORI_*.final.fa ./intermediate_files
rm *.final_for_sativa.tax
mv MIDORI_*.amp_blast.noN.mafft_edit.fa ./intermediate_files
rm *.final_rename_seqs_sativa.txt
rm MIDORI*.taxonomy_status.txt
rm MIDORI*.seqs_to_retain.txt
rm MIDORI*.final_id2acc2sp.txt
rm -f *.missing_sp_taxonomy.txt
rm *.missing_sp_to_delete.txt
rm *.consensus_taxonomy.txt
rm *.combined_taxonomy.txt
rm MIDORI.raw_id2sp.txt
rm MIDORI.raw_id2acc.txt
rm MIDORI*.raw_id2sp.txt
rm MIDORI*.raw_id2acc.txt


# Module 4 - Discard flagged sequences. Finalize consensus taxonomy and relabel sequences with correct species label and accession number. Select 1 representative sequence per haplotype per species.
cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
# Tetropoda
	# 12S
	# select accepted
	cut -f1 MIDORI_srRNA.Tetropoda_sativa/MIDORI_srRNA.Tetropoda.mis_to_delete | sort -k1,1 > MIDORI_srRNA.Tetropoda.sativa_flagged.txt
	tabtk isct -c -1 1 -1 1 MIDORI_srRNA.Tetropoda.sativa_flagged.txt <(grep ">" MIDORI_srRNA.Tetropoda.final_for_sativa.fa | sed 's/>//g') | seqtk subseq MIDORI_srRNA.Tetropoda.final_for_sativa.fa - > MIDORI_srRNA.Tetropoda.final_clean.fa
	# rename
	sed 's/;/\t/g' Tetropoda.final_taxonomy_sativa.txt | cut -f1,8 | sed 's/ /_/g' | awk '{print $1 "\t" $2 "_" $1}' > Tetropoda.final_rename_seqs_protax.txt
	seqkit replace -p '(.+)$' -r '{kv}' -k Tetropoda.final_rename_seqs_protax.txt MIDORI_srRNA.Tetropoda.final_clean.fa --keep-key > MIDORI_srRNA.Tetropoda.final_clean_relabel.fa
	#DY there is an upstream bug that i cannot find that creates duplicates of some sequences (e.g. Axis_porcinus_KY117537 gets 81 copies). This really slows down collapsetypes, so i remove these duplicated sequences here
	seqkit rmdup -n MIDORI_srRNA.Tetropoda.final_clean_relabel.fa > MIDORI_srRNA.Tetropoda.final_clean_relabel_rmdup.fa
	mv MIDORI_srRNA.Tetropoda.final_clean_relabel_rmdup.fa MIDORI_srRNA.Tetropoda.final_clean_relabel.fa
	seqbuddy MIDORI_srRNA.Tetropoda.final_clean_relabel.fa --clean_seq > MIDORI_srRNA.Tetropoda.final_clean_relabel.unalign.fa
    	#DY seqbuddy has a python error that adds "FTP Error: got more than 8192 bytes" to the end of ${label}.final_clean_relabel.unalign.fa
    	#DY add this line to remove any instances of "FTP Error: got more than 8192 bytes"
    	sed -i '/FTP Error: got more than 8192 bytes/d' MIDORI_srRNA.Tetropoda.final_clean_relabel.unalign.fa
    	sed -i '/FTP Error: \[Errno 8] nodename nor servname provided, or not known/d'  MIDORI_srRNA.Tetropoda.final_clean_relabel.unalign.fa
	# make final taxonomy
	cut -f2 Tetropoda.final_taxonomy_sativa.txt | cut -f1,2,3,4,5,7 -d";" | sed 's/;/ /g' | sort -u > Tetropoda.final_protax_taxonomy.txt
	# get list of sp
	grep ">" MIDORI_srRNA.Tetropoda.final_clean_relabel.unalign.fa | sed 's/>//g' | awk 'BEGIN{FS="_"}{print $1 "_" $2}' | sort -u > MIDORI_srRNA.Tetropoda.sp_list.txt
	# for each sp in list, if has >1 seq make temp fasta and collapse haplotypes, else print 1 seq
	bash Module4_collapsetypes.sh Tetropoda MIDORI_srRNA
	# end up with mix of upper/lower case and one/multi-line seqs, change to single-line uppercase
	awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' MIDORI_srRNA.Tetropoda.uniq.fa | seqtk seq -l0 - > tmp
	mv tmp Tetropoda.final_database.12S.fa

	# 16S
	# select accepted
	cut -f1 MIDORI_lrRNA.Tetropoda_sativa/MIDORI_lrRNA.Tetropoda.mis_to_delete | sort -k1,1 > MIDORI_lrRNA.Tetropoda.sativa_flagged.txt
	tabtk isct -c -1 1 -1 1 MIDORI_lrRNA.Tetropoda.sativa_flagged.txt <(grep ">" MIDORI_lrRNA.Tetropoda.final_for_sativa.fa | sed 's/>//g') | seqtk subseq MIDORI_lrRNA.Tetropoda.final_for_sativa.fa - > MIDORI_lrRNA.Tetropoda.final_clean.fa
	# rename
	sed 's/;/\t/g' Tetropoda.final_taxonomy_sativa.txt | cut -f1,8 | sed 's/ /_/g' | awk '{print $1 "\t" $2 "_" $1}' > Tetropoda.final_rename_seqs_protax.txt
	seqkit replace -p '(.+)$' -r '{kv}' -k Tetropoda.final_rename_seqs_protax.txt MIDORI_lrRNA.Tetropoda.final_clean.fa --keep-key > MIDORI_lrRNA.Tetropoda.final_clean_relabel.fa
	#DY there is an upstream bug that i cannot find that creates duplicates of some sequences (e.g. Axis_porcinus_KY117537 gets 81 copies). This really slows down collapsetypes, so i remove these duplicated sequences here
	seqkit rmdup -n MIDORI_lrRNA.Tetropoda.final_clean_relabel.fa > MIDORI_lrRNA.Tetropoda.final_clean_relabel_rmdup.fa
	mv MIDORI_lrRNA.Tetropoda.final_clean_relabel_rmdup.fa MIDORI_lrRNA.Tetropoda.final_clean_relabel.fa
	seqbuddy MIDORI_lrRNA.Tetropoda.final_clean_relabel.fa --clean_seq > MIDORI_lrRNA.Tetropoda.final_clean_relabel.unalign.fa
    	#DY seqbuddy has a python error that adds "FTP Error: got more than 8192 bytes" to the end of ${label}.final_clean_relabel.unalign.fa
    	#DY add this line to remove any instances of "FTP Error: got more than 8192 bytes"
    	sed -i '/FTP Error: got more than 8192 bytes/d' MIDORI_lrRNA.Tetropoda.final_clean_relabel.unalign.fa
    	sed -i '/FTP Error: \[Errno 8] nodename nor servname provided, or not known/d'  MIDORI_lrRNA.Tetropoda.final_clean_relabel.unalign.fa
	# make final taxonomy
	cut -f2 Tetropoda.final_taxonomy_sativa.txt | cut -f1,2,3,4,5,7 -d";" | sed 's/;/ /g' | sort -u > Tetropoda.final_protax_taxonomy.txt
	# get list of sp
	grep ">" MIDORI_lrRNA.Tetropoda.final_clean_relabel.unalign.fa | sed 's/>//g' | awk 'BEGIN{FS="_"}{print $1 "_" $2}' | sort -u > MIDORI_lrRNA.Tetropoda.sp_list.txt
	# for each sp in list, if has >1 seq make temp fasta and collapse haplotypes, else print 1 seq
	bash Module4_collapsetypes.sh Tetropoda MIDORI_lrRNA
	# end up with mix of upper/lower case and one/multi-line seqs, change to single-line uppercase
	awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' MIDORI_lrRNA.Tetropoda.uniq.fa | seqtk seq -l0 - > tmp
	mv tmp Tetropoda.final_database.16S.fa

# Teleostei
	# 12S
	# select accepted
	cut -f1 MIDORI_srRNA.Teleostei_sativa/MIDORI_srRNA.Teleostei.mis_to_delete | sort -k1,1 > MIDORI_srRNA.Teleostei.sativa_flagged.txt
	tabtk isct -c -1 1 -1 1 MIDORI_srRNA.Teleostei.sativa_flagged.txt <(grep ">" MIDORI_srRNA.Teleostei.final_for_sativa.fa | sed 's/>//g') | seqtk subseq MIDORI_srRNA.Teleostei.final_for_sativa.fa - > MIDORI_srRNA.Teleostei.final_clean.fa
	# rename
	sed 's/;/\t/g' Teleostei.final_taxonomy_sativa.txt | cut -f1,8 | sed 's/ /_/g' | awk '{print $1 "\t" $2 "_" $1}' > Teleostei.final_rename_seqs_protax.txt
	seqkit replace -p '(.+)$' -r '{kv}' -k Teleostei.final_rename_seqs_protax.txt MIDORI_srRNA.Teleostei.final_clean.fa --keep-key > MIDORI_srRNA.Teleostei.final_clean_relabel.fa
	#DY there is an upstream bug that i cannot find that creates duplicates of some sequences (e.g. Axis_porcinus_KY117537 gets 81 copies). This really slows down collapsetypes, so i remove these duplicated sequences here
	seqkit rmdup -n MIDORI_srRNA.Teleostei.final_clean_relabel.fa > MIDORI_srRNA.Teleostei.final_clean_relabel_rmdup.fa
	mv MIDORI_srRNA.Teleostei.final_clean_relabel_rmdup.fa MIDORI_srRNA.Teleostei.final_clean_relabel.fa
	seqbuddy MIDORI_srRNA.Teleostei.final_clean_relabel.fa --clean_seq > MIDORI_srRNA.Teleostei.final_clean_relabel.unalign.fa
    	#DY seqbuddy has a python error that adds "FTP Error: got more than 8192 bytes" to the end of ${label}.final_clean_relabel.unalign.fa
    	#DY add this line to remove any instances of "FTP Error: got more than 8192 bytes"
    	sed -i '/FTP Error: got more than 8192 bytes/d' MIDORI_srRNA.Teleostei.final_clean_relabel.unalign.fa
    	sed -i '/FTP Error: \[Errno 8] nodename nor servname provided, or not known/d'  MIDORI_srRNA.Teleostei.final_clean_relabel.unalign.fa
	# make final taxonomy
	cut -f2 Teleostei.final_taxonomy_sativa.txt | cut -f1,2,3,4,5,7 -d";" | sed 's/;/ /g' | sort -u > Teleostei.final_protax_taxonomy.txt
	# get list of sp
	grep ">" MIDORI_srRNA.Teleostei.final_clean_relabel.unalign.fa | sed 's/>//g' | awk 'BEGIN{FS="_"}{print $1 "_" $2}' | sort -u > MIDORI_srRNA.Teleostei.sp_list.txt
	# for each sp in list, if has >1 seq make temp fasta and collapse haplotypes, else print 1 seq
	bash Module4_collapsetypes.sh Teleostei MIDORI_srRNA
	# end up with mix of upper/lower case and one/multi-line seqs, change to single-line uppercase
	awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' MIDORI_srRNA.Teleostei.uniq.fa | seqtk seq -l0 - > tmp
	mv tmp Teleostei.final_database.12S.fa

	# 16S
	# select accepted
	cut -f1 MIDORI_lrRNA.Teleostei_sativa/MIDORI_lrRNA.Teleostei.mis_to_delete | sort -k1,1 > MIDORI_lrRNA.Teleostei.sativa_flagged.txt
	tabtk isct -c -1 1 -1 1 MIDORI_lrRNA.Teleostei.sativa_flagged.txt <(grep ">" MIDORI_lrRNA.Teleostei.final_for_sativa.fa | sed 's/>//g') | seqtk subseq MIDORI_lrRNA.Teleostei.final_for_sativa.fa - > MIDORI_lrRNA.Teleostei.final_clean.fa
	# rename
	sed 's/;/\t/g' Teleostei.final_taxonomy_sativa.txt | cut -f1,8 | sed 's/ /_/g' | awk '{print $1 "\t" $2 "_" $1}' > Teleostei.final_rename_seqs_protax.txt
	seqkit replace -p '(.+)$' -r '{kv}' -k Teleostei.final_rename_seqs_protax.txt MIDORI_lrRNA.Teleostei.final_clean.fa --keep-key > MIDORI_lrRNA.Teleostei.final_clean_relabel.fa
	#DY there is an upstream bug that i cannot find that creates duplicates of some sequences (e.g. Axis_porcinus_KY117537 gets 81 copies). This really slows down collapsetypes, so i remove these duplicated sequences here
	seqkit rmdup -n MIDORI_lrRNA.Teleostei.final_clean_relabel.fa > MIDORI_lrRNA.Teleostei.final_clean_relabel_rmdup.fa
	mv MIDORI_lrRNA.Teleostei.final_clean_relabel_rmdup.fa MIDORI_lrRNA.Teleostei.final_clean_relabel.fa
	seqbuddy MIDORI_lrRNA.Teleostei.final_clean_relabel.fa --clean_seq > MIDORI_lrRNA.Teleostei.final_clean_relabel.unalign.fa
    	#DY seqbuddy has a python error that adds "FTP Error: got more than 8192 bytes" to the end of ${label}.final_clean_relabel.unalign.fa
    	#DY add this line to remove any instances of "FTP Error: got more than 8192 bytes"
    	sed -i '/FTP Error: got more than 8192 bytes/d' MIDORI_lrRNA.Teleostei.final_clean_relabel.unalign.fa
    	sed -i '/FTP Error: \[Errno 8] nodename nor servname provided, or not known/d'  MIDORI_lrRNA.Teleostei.final_clean_relabel.unalign.fa
	# make final taxonomy
	cut -f2 Teleostei.final_taxonomy_sativa.txt | cut -f1,2,3,4,5,7 -d";" | sed 's/;/ /g' | sort -u > Teleostei.final_protax_taxonomy.txt
	# get list of sp
	grep ">" MIDORI_lrRNA.Teleostei.final_clean_relabel.unalign.fa | sed 's/>//g' | awk 'BEGIN{FS="_"}{print $1 "_" $2}' | sort -u > MIDORI_lrRNA.Teleostei.sp_list.txt
	# for each sp in list, if has >1 seq make temp fasta and collapse haplotypes, else print 1 seq
	bash Module4_collapsetypes.sh Teleostei MIDORI_lrRNA
	# end up with mix of upper/lower case and one/multi-line seqs, change to single-line uppercase
	awk '{if($0 !~ />/) {print toupper($0)} else {print $0}}' MIDORI_lrRNA.Teleostei.uniq.fa | seqtk seq -l0 - > tmp
	mv tmp Teleostei.final_database.16S.fa

#cleanup
mv *.final_for_sativa.fa ./intermediate_files
mv *.final_clean_relabel.unalign.fa ./intermediate_files
rm *.sativa_flagged.txt
rm *.final_clean.fa
rm *.final_clean_relabel.fa
rm *.sp_list.txt
mv *.final_taxonomy_sativa.txt ./intermediate_files
mv MIDORI_*_sativa/ ./intermediate_files
rm *.final_rename_seqs_protax.txt
rm MIDORI_*.uniq.fa

# Module 4 complete. You have reached the end of get_sequences.sh
#
# Actions after Module 4 complete:
# combining the results of Tetropoda and Teleostei
cat Tetropoda.final_database.16S.fa Teleostei.final_database.16S.fa > Vertebrata.final_database.16S.fa
cat Tetropoda.final_database.12S.fa Teleostei.final_database.12S.fa > Vertebrata.final_database.12S.fa
cat Tetropoda.final_protax_taxonomy.txt Teleostei.final_protax_taxonomy.txt > Vertebrata.final_protax_taxonomy.txt
# make copies
cp Tetropoda.final_database.16S.fa archived_files/Tetropoda.final_database.16S_orig202302.fa
cp Tetropoda.final_database.12S.fa archived_files/Tetropoda.final_database.12S_orig202302.fa
cp Teleostei.final_database.16S.fa archived_files/Teleostei.final_database.16S_orig202302.fa
cp Teleostei.final_database.12S.fa archived_files/Teleostei.final_database.12S_orig202302.fa
# remove any Homo_heidelbergensis sequences
seqkit grep Vertebrata.final_database.16S.fa -r -p Homo_heidelbergensis -v -o Vertebrata.final_database.16S_new.fa
mv Vertebrata.final_database.16S_new.fa Vertebrata.final_database.16S.fa
seqkit grep Vertebrata.final_database.12S.fa -r -p Homo_heidelbergensis -v -o Vertebrata.final_database.12S_new.fa
mv Vertebrata.final_database.12S_new.fa Vertebrata.final_database.12S.fa
rm -f Tetropoda.final_database.*.fa; rm -f Teleostei.final_database.*.fa
mv Tetropoda.final_protax_taxonomy.txt archived_files/; mv Teleostei.final_protax_taxonomy.txt archived_files/
# Final database sequences are in Vertebrata.final_database.locus.fa
# Final taxonomy file is in Vertebrata.final_protax_taxonomy.txt
#
# Next step: train PROTAX models with either:
#   - train_protax.sh for unweighted models
#   - train_weighted_protax.sh for models weighted using a list of expected species

# If sativa is used to change genus and/or species names, it is possible that some of the reference sequences in Vertebrata.final_database.16S.fa and Vertebrata.final_database.12S.fa will have a few sequences without species names (e.g. >_DQ158435) or starting with _TAXCLUSTER (e.g. >__TAXCLUSTER161__Spea_bombifrons_AY523786)
# These are created by sativa:
     # WARNING: Following taxa share >60% indentical sequences und thus considered indistinguishable:
     # Eukaryota;Chordata;Aves;Passeriformes;Turdidae;Turdus;Turdus ruficollis
     # Eukaryota;Chordata;Aves;Passeriformes;Turdidae;Turdus;Turdus eunomus
     # Eukaryota;Chordata;Aves;Passeriformes;Turdidae;Turdus;Turdus atrogularis
     # Eukaryota;Chordata;Aves;Passeriformes;Turdidae;Turdus;Turdus naumanni
     # For the purpose of mislabels identification, they were merged into one taxon:
     # Eukaryota;Chordata;Aves;Passeriformes;Turdidae;Turdus;__TAXCLUSTER141__Turdus ruficollis
# These should be removed because they interfere with PROTAX train (PROTAX needs sequences in the reference dataset to have the format >Ablepharus_kitaibelii_AY308325)
# seqkit grep Vertebrata.final_database.12S.fa -r -p ^_ -v -o Vertebrata.final_database.12S_new.fa
# seqkit grep Vertebrata.final_database.16S.fa -r -p ^_ -v -o Vertebrata.final_database.16S_new.fa
# visually check the new fasta files and then
# The results for MIDORI_1.4 do not have the above mislabels, so ignore it


# 4. Train PROTAX models for target amplicon(s)
#   - *train_protax.sh* (unweighted) or *train_weighted_protax.sh* (weighted)
#   - *check_protax_training.sh* (makes bias-accuracy plots)
#   - *choose between weighted or unweighted.  I use weighted because we have a species list for Ailaoshan

# unweighted
cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
bash ~/src/screenforbio-mbc-23GLG/train_protax.sh Vertebrata.final_protax_taxonomy.txt ~/src/screenforbio-mbc-23GLG
     # usage: bash train_protax.sh taxonomy screenforbio
     # where:
     # taxonomy is the final protax-formatted taxonomy file from get_sequences.sh (e.g. Vertebrata.final_protax_taxonomy.txt)
     # uses fasta files output from module_four of get_sequences.sh:  taxon.final_database.locus.fa (e.g. Vertebrata.final_database.12S.fa)
     # screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)
# End of train_protax.sh
# Success
# This took a total of 3.26 hours (MIDORI 1.4)

# weighted by Gaoligongshan species list:  splist (Tetropoda, Teleostei)
cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
# prepare the species list "splist.csv" of your study area and check that your species are in the Vertebrata.final_protax_taxonomy.txt file

bash ~/src/screenforbio-mbc-23GLG/train_weighted_protax.sh splist.csv Vertebrata.final_protax_taxonomy.txt ~/src/screenforbio-mbc-23GLG/
     # usage: bash train_weighted_protax.sh splist taxonomy screenforbio
     # where:
     # splist is a list of expected species to use in weighting in the format Genus,species (e.g. Homo,sapiens)
     # taxonomy is the final protax-formatted taxonomy file from get_sequences.sh (e.g. Vertebrata.final_protax_taxonomy.txt)
     # screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)
     # note: will take the taxon from the protax taxonomy file name
     # note: assumes curated database FASTA files are in current directory and labelled with format taxon.final_database.locus.fa (e.g. Vertebrata.final_database.12S.fa)

# End of train_weighted_protax.sh
# Success
# This took a total of ~3 hours (MIDORI 1.4)
     # weighted_protax_training.R is hardcoded to run through all four loci (12S, 16S, Cytb, COI), so if one or more loci is missing (i.e. no Vertebrata.final_database.Cytb.fa), the script is halted and generates an error message, but the previous loci do complete successfully.
     # example error message when Cytb is not included. Just ignore this and go to the next step (select mcmc iterations)
          # Working on Cytb in folder ./w_model_Cytb/
          # Working on level1
          # Error in file(file, "rt") : cannot open the connection
          # Calls: read.xdata -> read.table -> file
          # In addition: Warning message:
          # In file(file, "rt") :
          #   cannot open file './w_model_Cytb/train.w1.xdat': No such file or directory
          # Execution halted

# Select an mcmc iteration for each of the four levels for each model and each marker (e.g. ./w_model_16S/w_mcmc1a-d, ./w_model_16S/w_mcmc2a-d, etc.) based on the training plots (labelled ./w_model_16S/weighted_training_plot_16S_level1a_MCMC.pdf, etc). Chains should be well-mixed and acceptance ratio as close to 0.44 as possible. Relabel the selected model as ./w_model_16S/mcmc1 ./w_model_16S/mcmc2 etc.

# For an example of how to choose the model, go to archived_files/protax_training_mcmc_output_16S_example/. There are 4 PDF files (training_plot_16S_level4{a,b,c,d}_MCMC.pdf). To choose amongst these, Panu Somervuo wrote:
     # "In all four cases a-d, the highest log posterior is very similar (around -7912), and also the coefficients corresponding to it (red dot) among a-d are very close to each other (i.e. mislabeling probability around 0.25 , beta1 around 0, beta2 around -40, beta3 around 4 and beta4 around -80, so I would say that all of them would give very similar classification results. Of course, when looking the traceplot of a, it seems that the MCMC has not converged properly, since in the beginning of the plot it is in a different regime. However, the parameter values corresponding to the largest posterior are similar as in b,c,d. I think if taking any one from b,c,or d, they would give very similar (or even identical) classification results."
     # The traceplots of 'mcmc a' are seen to wander by looking at the traceplots themselves and also at the histograms, which are skewed.
     # Acceptance ratio is on the second page of the PDFs

cd ~/src/screenforbio-mbc-23GLG/
. ~/.linuxify; which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed

# unweighted
     # 12S
MOD1CHOSEN12S="mcmc1c"
MOD2CHOSEN12S="mcmc2a"
MOD3CHOSEN12S="mcmc3b"
MOD4CHOSEN12S="mcmc4b"
     # 16S
MOD1CHOSEN16S="mcmc1b"
MOD2CHOSEN16S="mcmc2d"
MOD3CHOSEN16S="mcmc3b"
MOD4CHOSEN16S="mcmc4c"
mv ./model_12S/${MOD1CHOSEN12S} ./model_12S/mcmc1
mv ./model_12S/${MOD2CHOSEN12S} ./model_12S/mcmc2
mv ./model_12S/${MOD3CHOSEN12S} ./model_12S/mcmc3
mv ./model_12S/${MOD4CHOSEN12S} ./model_12S/mcmc4
mv ./model_16S/${MOD1CHOSEN16S} ./model_16S/mcmc1
mv ./model_16S/${MOD2CHOSEN16S} ./model_16S/mcmc2
mv ./model_16S/${MOD3CHOSEN16S} ./model_16S/mcmc3
mv ./model_16S/${MOD4CHOSEN16S} ./model_16S/mcmc4

# weighted
     # 12S
w_MOD1CHOSEN12S="w_mcmc1c"
w_MOD2CHOSEN12S="w_mcmc2b"
w_MOD3CHOSEN12S="w_mcmc3d"
w_MOD4CHOSEN12S="w_mcmc4c"
     # 16S
w_MOD1CHOSEN16S="w_mcmc1b"
w_MOD2CHOSEN16S="w_mcmc2d"
w_MOD3CHOSEN16S="w_mcmc3a"
w_MOD4CHOSEN16S="w_mcmc4d"
mv ./w_model_12S/${w_MOD1CHOSEN12S} ./w_model_12S/w_mcmc1
mv ./w_model_12S/${w_MOD2CHOSEN12S} ./w_model_12S/w_mcmc2
mv ./w_model_12S/${w_MOD3CHOSEN12S} ./w_model_12S/w_mcmc3
mv ./w_model_12S/${w_MOD4CHOSEN12S} ./w_model_12S/w_mcmc4
mv ./w_model_16S/${w_MOD1CHOSEN16S} ./w_model_16S/w_mcmc1
mv ./w_model_16S/${w_MOD2CHOSEN16S} ./w_model_16S/w_mcmc2
mv ./w_model_16S/${w_MOD3CHOSEN16S} ./w_model_16S/w_mcmc3
mv ./w_model_16S/${w_MOD4CHOSEN16S} ./w_model_16S/w_mcmc4

# Next step: Check model training with check_protax_training.sh
# usage: bash check_protax_training.sh modeldir taxon locus screenforbio
# where:
# modeldir is the path to a directory containing the protax model to be checked
# taxon is the taxon for which the model was generated (used for labelling only)
# locus is the locus for which the model was generated (used for labelling only)
# screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)

# unweighted
bash check_protax_training.sh model_12S Vertebrata 12S ~/src/screenforbio-mbc-23GLG/
bash check_protax_training.sh model_16S Vertebrata 16S ~/src/screenforbio-mbc-23GLG/

# weighted
bash check_protax_training.sh w_model_12S Vertebrata 12S ~/src/screenforbio-mbc-23GLG/
bash check_protax_training.sh w_model_16S Vertebrata 16S ~/src/screenforbio-mbc-23GLG/

# Each model check took ~0.11 hours (~7 mins)
# Plots can be found in model_12S/checktrain/unweighted_Vertebrata_12S_biasaccuracy.pdf
# Plots can be found in model_16S/checktrain/unweighted_Vertebrata_16S_biasaccuracy.pdf
# Plots can be found in w_model_12S/checktrain/weighted_Vertebrata_12S_biasaccuracy.pdf
# Plots can be found in w_model_16S/checktrain/weighted_Vertebrata_16S_biasaccuracy.pdf


# 5. Classify query sequences (reads or OTUs) with PROTAX
#   - Process raw data with read_preprocessing.sh (experimental design must follow that described in the manuscript) and classify the output with protax_classify.sh or weighted_protax_classify.sh as appropriate
#   - Classify OTUs with protax_classify_otus.sh or weighted_protax_classify_otus.sh as appropriate

# these are the pathnames to my OTU representative sequences
OTUS12S="12S_otu.fas"
# OTUS16S="16S_otu.fas"
echo ${OTUS12S}
# echo ${OTUS16S}

# unweighted
     # move protax output files to a single folder
# mkdir protaxmodels/
mkdir protaxmodels/
mv model_12S protaxmodels/
mv model_16S protaxmodels/
bash protax_classify_otus.sh ${OTUS12S} 12S protaxmodels ~/src/screenforbio-mbc-23GLG protaxout
# bash protax_classify_otus.sh ${OTUS16S} 16S protaxmodels ~/src/screenforbio-mbc-23GLG protaxout
     # usage: bash protax_classify_otus.sh otus locus protaxdir screenforbio outdir
     # where:
     # otus is the (path to) the OTU fasta to be processed (suffix should be ".fa")
     # locus is the target locus, must be one of: 12S, 16S, CYTB, COI. if you have more than one locus to analyse, run script once for each.
     # protaxdir is the path to a directory containing protax models and clean databases for all 4 loci
     # screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)
     # outdir is the name to give an output directory (inside current) (no slash at end)

# weighted
# move weighted protax output files to a single folder where the weighting is for GLG
mkdir w_protaxmodels_GLG/
mv w_model_12S w_protaxmodels_GLG/
mv w_model_16S w_protaxmodels_GLG/
bash weighted_protax_classify_otus.sh ${OTUS12S} 12S w_protaxmodels_GLG ~/src/screenforbio-mbc-23GLG w_protaxout
# bash weighted_protax_classify_otus.sh ${OTUS16S} 16S w_protaxmodels_GLG ~/src/screenforbio-mbc-23GLG w_protaxout
     # usage: bash weighted_protax_classify_otus.sh otus locus protaxdir screenforbio outdir
     # where:
     # otus is the (path to) the OTU fasta to be processed (suffix should be ".fa")
     # locus is the target locus, must be one of: 12S, 16S, CYTB, COI. if you have more than one locus to analyse, run script once for each.
     # protaxdir is the path to a directory containing weighted protax models and clean databases for all 4 loci
     # screenforbio is the path to the screenforbio-mbc directory (must contain subdirectory protaxscripts)
     # outdir is the name to give the output directory (inside current)

# Success
# Example output from 12S:
#
# Results are in ./w_protaxout_12S/
# Classification for each OTU at each taxonomic level (species, genus, family, order) in files 12S_otu.fas.w_<level>_probs
# e.g. 12S_otu.fas.w_species_probs
# queryID taxID   log(probability)  level   taxon
# OTU1    816     -1.25544          4       Anura,Dicroglossidae,Nanorana,taihangnica
# log(prob) is ln(prob), so to get probability, convert by exp(log(probability))

# Additionally, the best matching hit (for assigned species/genus where available) found with LAST is appended to 12S_otu.fas.w_species_probs in 12S_otu.fas.w_species_probs_sim
# queryID taxID   log(probability) level   taxon                                        bestHit_similarity      bestHit
# OTU1    816     -1.25544         4       Anura,Dicroglossidae,Nanorana,taihangnica    0.979 Nanorana_taihangnica_KJ569109

# 5.1 Use combine_protax_output_tables.Rmd to combine the protax output files. The script is written to process w_protaxout_12S
