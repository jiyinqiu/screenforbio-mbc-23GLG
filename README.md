# ScreenForBio metabarcoding pipeline for Gaoligong

A fork of the ScreenForBio metabarcoding pipeline originally published in:

Axtner, J., Crampton-Platt, A., Hoerig, L.A., Xu, C.C.Y., Yu, D.W., Wilting, A. (2019), An efficient and robust laboratory workflow and tetrapod database for larger scale environmental DNA studies. Giga Sci. 8, 646–17. doi:10.1093/gigascience/giz029

Preprint available from [bioRxiv](https://www.biorxiv.org/content/early/2018/06/12/345082).  
Paper available from [GigaScience](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giz029/5450733).  

This fork is customised to the Gaoligong eDNA dataset, and updated to take advantage of the MIDORI2_UNIQ_NUC_GB253 datasets that are now available. There are some bug fixes and additional R utilities. However, if you have some knowledge of bash and R scripting, then you can easily adapt this pipeline for your paired-end Illumina metabarcoding, from processing of twin-tagged raw reads through to taxonomic classification with *PROTAX*.

The pipeline steps that will be of most use to other projects are those related to *PROTAX*: generating curated reference databases and associated taxonomic classification; training *PROTAX* models (weighted or unweighted); and classification of OTUs with *PROTAX*.

Steps and associated scripts:
1. Process twin-tagged metabarcoding data (not run in this fork, but the script is available)
  - *read_preprocessing.sh*
2. Obtain initial taxonomic classification for target taxon
  - *get_taxonomy.sh*
3. Generate non-redundant curated reference sequence database for target amplicon(s) and fix taxonomic classification
  - *get_sequences.sh*
4. Train PROTAX models for target amplicon(s)
  - *train_protax.sh* (unweighted) or *train_weighted_protax.sh* (weighted)
  - *check_protax_training.sh* (makes bias-accuracy plots)
5. Classify query sequences (reads or OTUs) with PROTAX
  - *protax_classify.sh* or *protax_classify_otus.sh* (unweighted models)
  - *weighted_protax_classify.sh* or *weighted_protax_classify_otus.sh* (weighted models)

**Note:** in some steps the ***screenforbio-mbc*** release associated with the manuscript is specific to the amplicons used in the study - primer sets and relevant settings are hard-coded in *read_preprocessing.sh* and *get_sequences.sh*. This will be generalised in a future release.

### Required software (tested versions)
Pipeline tested on macOS 10.14.4 on a MacBook Pro with i7 CPU (4 cores, 8 virtual cores). If you run with an i5 CPU, you have only 4 virtual cores, and you will need to adjust the requested number of threads in the blastn (-num_threads) and sativa.py (-T) commands in get_sequences.sh

- Homebrew for macOS  
Go to http://brew.sh and follow the instructions for installing Homebrew on macOS

- After Homebrew or Linuxbrew  is installed, run these brew installations
````
brew tap brewsci/bio # a "tap" is a source of "installation formulae" of specialist software, here, bioinformatics
brew install brewsci/bio/seqkit
brew install brewsci/bio/last # v926+
brew install coreutils # cut, join, sort, tr, etc.
brew install seqtk # https://github.com/lh3/seqtk
brew install gnu-sed # GNU sed
brew install grep # GNU grep
brew install gawk # GNU awk
brew install perl # v5.28+
brew install blast # v2.9.0+
brew install mafft # v7.407+
brew install cutadapt # v2.3+
brew install python@3
brew install python@2
brew update; brew upgrade; brew cleanup  # run occasionally to update your software
````

For the rest of the software packages, here are installation instructions:  

- I install all github repositories (repos) in ~/src  
````
mkdir ~/src
````

- Install *linuxify* to prioritise GNU versions of sed, awk, and grep over the macOS versions. GNU grep, GNU sed, and GNU awk are installed with homebrew, but they are given different names (e.g. gsed, gawk, ggrep). However, the scripts use 'sed', 'grep', and 'awk'. To prioritise the GNU versions and to call them as sed, awk, and grep, i use 'Linuxify'  
1. install and run linuxify # https://github.com/fabiomaia/linuxify
````
     cd ~/src
     git clone https://github.com/fabiomaia/linuxify.git
     cd linuxify/
     ./linuxify install # to install all the GNU versions
          # cd ~/src/linuxify/;  ./linuxify uninstall # to remove all the GNU versions and the pathname changes
          # manually remove '. ~/.linuxify' from my ~/.bashrc file
     ls -al ~/src/linuxify/ # should see the file .linuxify
     cp ~/src/linuxify/.linuxify ~/ # cp to root directory
````  
2. to 'linuxify' a terminal session:  run the following at the beginning of a script or a session.
````
     . ~/.linuxify; awk; which sed; which grep
          # awk # should return help page for gawk
          # which sed # should show /usr/local/opt/gnu-sed/libexec/gnubin/sed
          # which grep # should return: '/usr/local/opt/grep/libexec/gnubin/grep'\
````  
3. OPTIONAL if i want to run linuxify automatically with each new shell
````
     # add this to my ~/.bashrc
          . ~/.linuxify
     # then add this to my .bash_profile
          if [ -f ~/.bashrc ]; then
              source ~/.bashrc
          fi
     # https://apple.stackexchange.com/questions/51036/what-is-the-difference-between-bash-profile-and-bashrc
````

bcl2fastq and AdapterRemoval are required for processing of raw reads.  
- bcl2fastq (v2.18)  
     * downloaded from Illumina.com (only for RedHat or CentOS Linux). If you have a local sequencer, it should already be running somewhere. If you send out for sequencing, the provider will have it running.  
- AdapterRemoval (v2.1.7)  
     * first install [miniconda](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg)  
````
     # brew install python@3 # installs python3 if you haven't done so already (see above)
     conda update -n base conda  
     conda install -c bioconda adapterremoval  
````
- usearch (v8.1.1861_i86osx32, v11.0.667_i86osx32)  
````
     # go to https://drive5.com/usearch/
     # register and download the 32-bit usearch11 binary for your OS
     # install usearch binary downloaded from https://drive5.com/usearch/download.html
     cd ~/Downloads
     mv ~/Downloads/usearch11.0.667_i86osx32 /usr/local/bin
     cd /usr/local/bin
     ln -s usearch11.0.667_i86osx32 usearch11  # make symbolic link (symlink)
     chmod 755 usearch11
     usearch11 # should return something like the following
          # usearch v11.0.667_i86osx32, 4.0Gb RAM (17.2Gb total), 8 cores
          # (C) Copyright 2013-18 Robert C. Edgar, all rights reserved.
          # https://drive5.com/usearch

     # download and install usearch 8.1, also from https://drive5.com/usearch/download.html
     mv ~/Downloads/v8.1.1861_i86osx32 /usr/local/bin
     cd /usr/local/bin
     ln -s v8.1.1861_i86osx32 usearch  # make symbolic link (symlink)
     chmod 755 usearch
     usearch # should return something like the following
          # usearch v8.1.1861_i86osx32, 4.0Gb RAM (17.2Gb total), 8 cores
          # (C) Copyright 2013-15 Robert C. Edgar, all rights reserved.
          # http://drive5.com/usearch
````
- RStudio
     * installed from binary downloaded from [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
- R (v3.6.0+)  
     * installed from binary downloaded from [CRAN](https://cran.rstudio.com)
- taxize (v0.9.8) R package.  On the R command line, run:  
     `install.packages("taxize", dep=TRUE)`  
- tabtk (r19)  
````
     cd ~/src/  
     git clone https://github.com/lh3/tabtk.git  
     cd tabtk  
     make  
     mv tabtk /usr/local/bin/tabtk  
````
- sativa (v0.9-57-g8a99328)  Sativa is unstable during installation, in that its built-in raxml source does not always compile. Check that the sativa/raxml/ directory contains one or more raxml binaries (the specific one depends on your configuration) and test the sativa installation with the example dataset (see code below).
````
     cd ~/src; git clone https://github.com/amkozlov/sativa  
     cd sativa
     ./install.sh --no-avx  
     ln -s sativa.py sativa  
     echo 'export PATH="$HOME/src/sativa:${PATH}"' >> ~/.bash_profile  
     # source ~/.bash_profile # if you want to run immediately  
          # to test installation,
               # ls ~/src/sativa/raxml/ # should contain one or more raxml binaries
               # cd ~/src/sativa/example
               # ../sativa.py -s test.phy -t test.tax -x BAC -T 2 # runs fast
````
- seqbuddy (v1.3.0)  
````
     pip3 install buddysuite # installs a bunch of software, requires python3
     # https://github.com/biologyguy/BuddySuite
     buddysuite -setup  
     seqbuddy -h
````
- Entrez Direct (v6.00 and v8.30)  
     * see installation instructions on [NIH](https://www.ncbi.nlm.nih.gov/books/NBK179288/)  

- *get_sequences.sh* also requires MIDORI databases for mitochondrial target genes [Machida *et al.*, 2017](https://www.nature.com/articles/sdata201727). Download relevant MIDORI_UNIQUE FASTAs in RDP format from the [website](http://www.reference-midori.info/download.php). The manuscript used MIDORI_UNIQUE_1.1 versions of COI, Cytb, lrRNA and srRNA. The unzipped FASTAs should be *copied* to the working directory (because the script moves the working MIDORI fasta files to the intermediate_files/ folder after it finishes module one).  There are downloaded versions in the archived_files/ folder

- The 20180221 versions of MIDORI have more complex headers, which interfere with the `get_sequences.sh` code.  
     * *V 1.1*:  `>AF382008	root;Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo;Homo sapiens`  
     * *V 20180221*:  `>AF382008.3.649.1602	root;Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo;Homo sapiens`  

     The filenames will be changed to this format: `MIDORI_UNIQUE_1.2_srRNA_RDP.fasta`, and the extra stuff on the headers will be removed before running *get_sequences.sh*

- *collapsetypes_v4.6.pl* should already be in your ***screenforbio-mbc*** directory. If not, install as follows:  
Download from Douglas Chesters' [sourceforge page](https://sourceforge.net/projects/collapsetypes/).  
````
     chmod 755 ~/Downloads/collapsetypes_v4.6.pl  
     mv ~/Downloads/collapsetypes_v4.6.pl ~/src/screenforbio-mbc-23GLG/  
````

- *PROTAX* scripts are reposted here with the kind permission of Panu Somervuo. These are in the *protaxscripts* subdirectory of ***screenforbio-mbc-23GLG***. This version of *PROTAX* is from [Rodgers *et al.* 2017](https://doi.org/10.1111/1755-0998.12701), and the scripts were originally posted on [Dryad](https://datadryad.org/resource/doi:10.5061/dryad.bj5k0).  

### Usage
All steps in the pipeline are implemented via bash scripts with similar parameter requirements. Each script includes commented usage instructions at the start and displays the same instructions if run without any or an incorrect number of parameters.

For the Ailaoshan fork, the full command history is in `screenforbio_23GLG.sh`. Although formatted as a shell script, it should be run command by command, instead of as a single bash script file, because there are multiple choices to be made during the pipeline.

Some of the bash scripts used within are primarily wrappers for R scripts, all of which are assumed to be in the ***screenforbio-mbc-23GLG*** directory.  

*get_taxonomy.sh* example:

    bash ~/src/screenforbio-mbc-23GLG/get_taxonomy.sh

You will see the following message:

    You are trying to use get_taxonomy.sh but have not provided enough information.
    Please check the following:

    Usage: bash get_taxonomy.sh taxonName taxonRank screenforbio
    Where:
    taxonName is the scientific name of the target taxon e.g. Tetrapoda
    taxonRank is the classification rank of the target taxon e.g. superclass
    screenforbio is the path to the screenforbio-mbc-23GLG directory

To get the ITIS classification for Mammalia:

    bash ~/src/screenforbio-mbc-23GLG/get_taxonomy.sh Mammalia class ~/src/screenforbio-mbc-23GLG/

When the script is running various messages will be printed to both the terminal and a log file. On completion of the script the final messages will indicate the next script that should be run and any actions the user should take beforehand.

Enjoy :-)
