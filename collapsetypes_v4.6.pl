#! /usr/bin/perl

# Name: Collapsetypes.pl
# Description: Perl script for removing identical or near identical DNA sequences from a fasta file. Creates a file with only unique haplotypes / genotypes
# Copyright (C) 2013  Douglas Chesters
# Modified by Enrique Gonzalez-Tortuero
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# contact address dc0357548934@live.co.uk
#
# citation: Chesters, D. (2013) collapsetypes.pl [computer software available at http://sourceforge.net/projects/collapsetypes/]
#
# QUICKSTART:
# -----------
#
# Run with the command:
# $ perl collapsetypes_v4.5.pl -infile=input.fas -ambchar=?nN -nrdiffs=3
#
# Command line options are input file name, the characters you wish to depict ambiguous positions and the number of differences to determine haplotypes.
# * Input is fasta format alignment.
# * Ambiguous positions (-ambchar=<characters>) are those in which the type of nucleotide is unknown, here these are ignored when deciding whether 2 sequences are identical. If you dont have ambiguous characters, it doesnt matter what you use in the command line.
# * The number of differences to determine haplotypes (-nrdiffs=<number>) allows sequencing misreads, i.e if 2 sequences differ by just a single base, they are considered the same haplotype. By default, this program consider that a difference between sequences are the same haplotype.
# * Output is a file with only unique sequences in it.
#
# MORE DETAILED:
# --------------
#
# If you have a Mac or Linux:
# 	* Just open a shell (terminal), add execute permission, and change to the directory containing your input file. If you don't want to add execute permission to this script, you must to change to the folder where the script is.
#
# If you have Windows: 
# 	* You will need to install Perl first. Go to the activePerl website once installed put perl in your dos path (google it), or just place your input file into the /perl/bin/ folder. Then, open a command prompt then cd into the folder with your input file
#
# Type: "perl collapsetypes_v4.5.pl -infile=[input_file_name] -ambchar=[additional_ambiguity_characters] -nrdiffs=[number_of_allowed_diffs]"
# (Eg: perl collapsetypes_v4.5.pl -infile=input.fas -ambchar=?nN -nrdiffs=3)
#
# * Input file must be fasta format ALIGNMENT.
#
# * Additional ambiguity characters I have come across: '?', 'x', 'X', 'n', 'N', '-' in particular, decide what you want to do about gap characters. For example if you have indels, represented by hyphen in the sequences, and you wish that the presence of an indel between 2 sequences would still constitute an identical pair of sequences, then you would include hyphen as an additional ambiguity characters in the command. Conversly, if you think an indel means the 2 sequences are different (this is default), then no hyphen in the command 
# * The number of differences to determine haplotypes (-nrdiffs=<number>) allows sequencing misreads, i.e if 2 sequences differ by just a single base, they are considered the same haplotype. By default, this program consider that a difference between sequences are the same haplotype.
# * output files: 	1) fasta format file with the unique sequences
#			2) haplotype_assignations, which sequences have been printed, with list of identical sequences on the same row (see note below)
#			3) haplotype_assignations2, each unique haplotype is given a number, file lists all sequences and what haplotype number they are, and whether printed
#			4) haplotype_assignations3, lists all sequences which are duplicates, and which other sequences they match. Pattern of matches is not always consistent, see explaination for this, in the section 'this script will sometimes give you unintuitive .....'
#
# The script removes redundent sequences from a input file, printing one of each unique haplotype to the output file other programs that do something similar include:
#	* TCS (Clement, 2000, Mol Ecol)
#	* FaBox (Villesen 2007, mol ecol notes), 
#	* DNASP (Rozas et al 2003)
#	* Collapse (Posada 2006) - this program differs in that the full nucleotide code is accounted for (see algorithm below), and does not require installation (most of the above do require this) and generates some useful additional output files and permits 'near'-identical sequence removal also. Useful if you want to account for some sequencing error. In terms of general behaviour, i have had reports that it is more liberal at sequence removal than Collapse. Also, when used prior to GMYC analysis (Pons et al 06, syst biol), I hear it gives significant improvement to likelihoods the script starts to get very slow as you use large numbers of seqs (>~2000)
#
# The algorithm:
#
# loop 1, read each sequence in turn
# loop2, read each sequence in turn, starting from the index of above sequence +1. So all sequence pairs are compared.
# For each pairwise comparison:
# 	for each occurance of   N,       if corresponding position in other sequence has	[atcgmrwsykvhdb],	other seq gets score+1
#              ..           	Y                               ..                      	[ct]         		..
#              ..           	R                               ..  	                        [ag]			..	
#              ..           	M                               .. 	                        [ac]			..
#              ..           	W                               ..  	                        [at]			..
#              ..           	S                               ..		                [cg]			..
#              ..           	K                               ..		                [tg]			..
#              ..           	V                               ..		                [acgmsr]		..
#              ..           	H                               ..		                [actmyw]		..
#              ..           	D                               .. 	                        [agtrkw]		..
#              ..           	B                               ..    	                        [cgtsky]		..
#
# For all positions except those with [NYRMWSKVHBD], if the 2 sequences are identical put black mark against lowest scoring sequence. If they score the same then black mark the pair member in loop 1 (this allows further comparisons with the other member). Then, print all sequence without black mark to file, unique_haplotypes.
#
#
# This script will sometimes give you unintuitive results, heres why:
#
# 	seq1 CCCGATCGTTCCC
# 	seq2 ATCGATCGTTCCC
# 	seq3 NNNNYYYGTTCCC
#
# When deciding whether seq1 and 3 above are the same haplotype, it can only look at the last 6 bases of each sequence, and indeed they are the same here (also the same as seq2). so as far as the program is concerned seq 3 is the same as seq2 and seq3 is the same as seq1. however the lack of N's in seq 1 and 2 mean theres more bases to compare, and they are not the same at some of these. so the result is seq1 == seq3. seq2 == seq3 even though seq1 does not == seq2.
# So in the haplotype assignations output a sequence composed of largly missing data will probably be assigned to lots of different haplotypes. But will not be printed anyway because it will be low scoring
#
# There should be mp duplicate id's in your input file. Assignations3 file lists all duplicated id's as assignations1 / 2 will ommit if one is same as another which is same as another which is removed. Sequence id's need to be reasonable length. not e.g. 1 character
#
# Its IMPORTANT to consider the characters you use for gap and missing for default i have  "N" , "n" , "x" , "X" , "?" for missing, and assume "-" for gap, treated as different character state, that would mean that, in the example below, seq5 and seq6 are the same haplotype, seq6 and seq7 are different haplotypes
#
# 	seq5 CCCCC
# 	seq6 CCCCN
# 	seq7 -CCCC
#
# Note that N is always used as an ambiguity character, irrespective of what you specify in the command line.
#
#
# Main test dataset:
#
#	>SEQ_c.3
#	?aaaaaaaaaaaaaatt
#	>SEQ_d.1
#	aaaaaaaaaaaaaaata
#	>SEQ_c.1
#	aaaaaaaaaaaaaaatt
#	>SEQ_a.1
#	ttttttttttttttttt
#	>SEQ_c.2
#	aaaaaaaaaaaaaaatt
#	>SEQ_a.2
#	ttttttttttttttytt
#	>SEQ_a.3
#	ttttttttttttttttt
#	>SEQ_b.1
#	ttttttttttttttttc
#	>SEQ_c.3
#	?aaaaaaaaaaaaaatn
#
#
# the three SEQ_a's are all the same, but SEQ_a.2 has an ambig char so should not be printed.
# SEQ_c.3 is a duplicated ID. SEQ_c's all the same, but SEQ_c.3 is the worse scoring of the three.
# SEQ_d.1 the same as second SEQ_c.3 when nonAbmig chars are looked at so the lower SEQ_c.3 belongs to 2 haplotype numbers
#
# output:
#
# 	printed:SEQ_d.1	same_as:SEQ_c.3	
# 	printed:SEQ_c.2	same_as:SEQ_c.3	SEQ_c.1	SEQ_c.3	
# 	printed:SEQ_a.3	same_as:SEQ_a.1	SEQ_a.2	
# 	printed:SEQ_b.1
#
# VERSION HISTORY
# ---------------
#
# 	v4: 	Will now work on files produced across platforms. 
#		Eg you can run the script on a pc, using a fasta file generated on mac
#		Does not replace spaces with underscores on fasta ids
#		Can now handle identical id's in the input file (although these shouldnt really be present), by appending numbers to redundent id's. This feature works ONLY WITH FASTA FORMAT, not with nexus (its more difficult with nexus due to possibilty of sequence interleaving)
#		Tabs as seperators instead of spaces in information files, since id's may contain spaces this would disorganise excel's data to columns feature
# 	v4.1: Takes nexus format (tested with sequential version, but interleaved should work). 
#		Accepts x, N and ? as missing data char
#		Note comments (anything between square brackets) are removed from nexus files. 
#		I know sometimes extra sequence is included in these comments. if you are desperate to have comments reinserted into outfile let me know.
# 	v4.2: '-' added as missing data character. writes output in nexus format
# 	v4.3: big re-write. now based on sequences id's rather than file indices. 
# 	            Another output file (.assign3) is given, 
# 	            SOME INFORMATION OMITTED FROM OTHER FILES SHOULD BE FOUND IN THIS.
# 	v4.4: User specified missing data / gap characters
# 	v4.5: Missing characters specified in command line. nexus format support discontinued. streamlined.
#	v4.6: improvement to UI and structure, by Enrique Gonzalez-Tortuero

## Using core modules
use strict;
use warnings;

## Starting variables of the program
# External parameters
my $command_line_specified_ambiguity_characters;
my $collapse_when_identical;
my $collapsing_cutoff;
my $input_sequences_fasta;
for (my $i = 0; $i < scalar @ARGV; $i++) {
	if ($ARGV[$i] =~ /^-infile=(\S+)/) {
		$input_sequences_fasta = $1;
	}
	if ($ARGV[$i] =~ /^-ambchar=(\S+)/){
		$command_line_specified_ambiguity_characters = $1;
	}
	if ($ARGV[$i] =~ /^-nrdiffs=(\d+)/) {
		$collapsing_cutoff = $1; # This is the cutoff.
		if ($collapsing_cutoff == 0) {
			$collapse_when_identical = 1; # 1 == default (collapse ONLY when identical).
		} else {
			$collapse_when_identical = 0; #  0 == collapse also, when they differ by ... bases
		}
	}
}

# Internal parameters
my @missing_data_characters;
my %standard_ambiguity_characters = 
	(
	N => 'ATCGMRWSYKVHDB',
	Y => 'CT',
	R => 'AG',
	M => 'AC',
	W => 'AT',
	S => 'CG',
	K => 'TG',
	V => 'ACGMSR',
	H => 'ACTMYW',
	D => 'AGTRKW',
	B => 'CGTSKY'
	);
my $file_as_string = "";
my @all_lines = ();
my @all_sequence_ids_array;
my %sequences_from_formatted_file;
my %original_sequences;
my %user_IDs_key;
my $no_sequences = 0;
my $current_seq_id;
my $comparison_seq_id;
my $current_seq;
my $current_seq_copy;
my $comparison_seq;
my $comparison_seq_copy;
my $comp_score;
my $curr_score;
my %unique_hash = ();
my %id_list_array;

## Print header of the program
print "\n********\tcollapsetypes_v4.5\n********\tdc0357548934\@live.co.uk for bug reports & suggestions\n\n";

## Checking external parameters
# Checking input file to avoid strange symbols
unless($input_sequences_fasta =~ /[a-zA-Z0-9]/){
	die "\nFilename you specified is strange! quitting without analysis.\n"
}

# Checking the file to avoid nexus files!
print "Looking for input file: $input_sequences_fasta\n\n";
open(INPUT, "$input_sequences_fasta") || die "i cant find input file. it should be named $input_sequences_fasta\n";

while(my $line= <INPUT>){
	$file_as_string .= $line;
}
close(INPUT);

if ($file_as_string =~ /^.nexus/i){
	die "\n\nFasta format is prefered. Nexus support has been discontinued. Quitting.\n";
}

@all_lines = split /\012\015?|\015\012?/, $file_as_string;
$file_as_string = ""; # this split command allows cross-platform use

# Checking ambiguity characters
if (defined($command_line_specified_ambiguity_characters)) {
	@missing_data_characters = split // ,$command_line_specified_ambiguity_characters;
} else {
	@missing_data_characters = ("N" , "n" , "x" , "X" , "?");
	warn "WARNING: You did not specify ambiguity characters. Using default (N, n, X, x, ?).\n";
	warn "This may come back to bite you later. Please read notes on why correct specification of ambiguity characters is important.\n\n\n";
}

## Printing the haplotyping definitions
print "Processing input file $input_sequences_fasta\nWith the following ambiguity characters: @missing_data_characters\n";
if ($collapse_when_identical == 0) {
	print "Considering all sequences with $collapsing_cutoff differences as the same haplotype.\n"
}

## Creating outputs
open(OUTPUT,">$input_sequences_fasta.unique_haplotypes") || die "cant open output\n";
open(OUT3,">$input_sequences_fasta.haplotype_assignations") || die "cant open output\n";
open(OUT4,">$input_sequences_fasta.haplotype_assignations2") || die "cant open output\n";

print OUT4 "id\thaplotype_number\tprinted_to_outfile\n";

## Reading the FASTA file
open(INPUT, "$input_sequences_fasta") || die "cant open input:$input_sequences_fasta\n";

# Defining variables for this step
my $current_id;
my $entry_counter = 0;
my $file_as_string2 = "";

while(my $line = <INPUT>) {
	$file_as_string .= $line;
}

close(INPUT);

my @all_lines3 = split /\012\015?|\015\012?/, $file_as_string;
$file_as_string2 = join ( "\n" , @all_lines3 );

# This currently has id's AND seqs, but seqs to be removed and put in hash.
@all_sequence_ids_array = split />/, $file_as_string2;
$file_as_string2 = "";

for my $each_line(1 .. $#all_sequence_ids_array) { # as it starts with > and is split at this char, 1st index of array is empty
	my $line = $all_sequence_ids_array[$each_line];
	$entry_counter++;
	my $temporary_sequence_ID = "ID$entry_counter" . "ID";
	if ($line =~ /^(.+)\n/ ) {
		$current_id = $1;
		$line =~ s/^.+\n//;
		$line =~ s/\n//;
		$line =~ s/\r//;
		$line =~ s/\s//g;
		$line =~ s/\t//g;
	} else {
		die "BUG\n$line\n";
	}

	$original_sequences {$temporary_sequence_ID} = $line;
	$user_IDs_key{$temporary_sequence_ID} = $current_id;
	$line =~ tr/a-z/A-Z/; # /g? ..... not req.

	foreach my $missing_data_char(@missing_data_characters) {
		$line =~ s/[$missing_data_char]/N/g;
	}

	$sequences_from_formatted_file{$temporary_sequence_ID} = $line;
	$all_sequence_ids_array[$each_line] = $temporary_sequence_ID;
	$no_sequences++;

}
close(INPUT);
print "$input_sequences_fasta has $entry_counter sequences.\n";

## Pairwise comparisons
for my $comparison_index(1 .. ($no_sequences - 1)) {

	if ($comparison_index % 10 == 0) {
		print "sequence number $comparison_index out of $no_sequences\n"
	}

	my $entry_counter = 0;

	for my $index( ($comparison_index + 1) .. ($no_sequences - 0)) {

		$current_seq_id = $all_sequence_ids_array[$index];
		$comparison_seq_id = $all_sequence_ids_array[$comparison_index];
		$current_seq = $sequences_from_formatted_file{ $current_seq_id };
		$comparison_seq = $sequences_from_formatted_file{ $comparison_seq_id };
		$comparison_seq_copy=$comparison_seq;
		$current_seq_copy=$current_seq;
		$comp_score=0;
		$curr_score=0;
		my $sub_output = score_sequence_pair( $comparison_seq_copy . "___" . $current_seq_copy );
		my @split_sub_output = split(/___/, $sub_output);
		my $comparison_seq_copy = $split_sub_output[0];
		my $current_seq_copy = $split_sub_output[1];
		my $collapsing_cutoff_copy = $collapsing_cutoff;
		if($collapse_when_identical == 0) {
			for my $seq_index22 (0..length($current_seq_copy)) {
				if($collapsing_cutoff_copy >= 1) {
					my $test_char = substr($current_seq_copy,$seq_index22,1);
					my $test_char2 = substr($comparison_seq_copy,$seq_index22,1);
					unless($test_char eq $test_char2){
						substr($current_seq_copy,$seq_index22,1) = $test_char2;
						$collapsing_cutoff_copy--;
	                                }
				}
			}
		}

		if($comparison_seq_copy eq $current_seq_copy) {
			if($comp_score<=$curr_score) {
				$unique_hash{$comparison_seq_id} = 1;
				$id_list_array{$current_seq_id} .= "\t" .  $comparison_seq_id;
			} else {
				$unique_hash{$current_seq_id}=1;
				$id_list_array{$comparison_seq_id} .= "\t" .  $current_seq_id;
			}
		}
		$entry_counter++;
	}
}

# The %unique_hash now contains all unique strings, go through the file again to print these indices out
$current_seq="";
my $count_number_unique = 0;
my $haplotype_number = 1;

open(INPUT, "$input_sequences_fasta") || die "cant open input:$input_sequences_fasta\n";

$entry_counter = 0;

for my $index ( 1 .. ($no_sequences - 0)) {

	$current_seq_id = $all_sequence_ids_array[$index];
	$current_seq = $sequences_from_formatted_file{ $current_seq_id };

	unless(exists $unique_hash{$current_seq_id}) {

		my $original_sequence 	= $original_sequences {$current_seq_id};
		my $original_fastaID 	= $user_IDs_key{$current_seq_id};

		unless(length($original_sequence)>=1 && length($original_fastaID)>=1){
			die "\n\nerror .... quitting\n"
		} 

		print OUTPUT ">$original_fastaID\n$original_sequence\n";
		$count_number_unique++;
		my $current_identicals = $id_list_array{$current_seq_id};
		$current_identicals =~ s/\n//g;
		$current_identicals =~ s/>//g;
		my $length_current_identicals = length($current_identicals);

		if($length_current_identicals >= 1) {
			unless($current_identicals =~ /[a-zA-Z0-9]/) {
				print "no list of identicals for seq:$current_seq_id\n"
			}
			print OUT3 "printed:$original_fastaID\tsame_as:";
			print OUT4 "$original_fastaID\t$haplotype_number\tyes\n";
			my @split_names = split(/\t/,$current_identicals);
			for my $each1(0 .. $#split_names) {
				if(length($split_names[$each1])>1) {
					my $originalID = $user_IDs_key{$split_names[$each1]};
					print OUT3 "$originalID\t";
					print OUT4 "$originalID\t$haplotype_number\tno\n";
				}
			}
			print OUT3 "\n";
		} else {
			print OUT3 "printed:$original_fastaID\n";
			print OUT4 "$original_fastaID\t$haplotype_number\tyes\n"
		}
		$haplotype_number++;
	}
	$entry_counter++;
}

print "\n$count_number_unique unique haplotypes have been printed\n";

close(INPUT);
close(OUTPUT);
close(OUT3);
close(OUT4);

my @allkeys = keys %id_list_array;
@allkeys = sort @allkeys;

open(OUT5,">$input_sequences_fasta.haplotype_assignations3") || die "cant open output\n";

foreach my $key(@allkeys) {
	my $original_fastaID 	= $user_IDs_key{$key};
	print OUT5 "\n$key\tsame_as:\t";

	my @split_names = split(/\t/, $id_list_array{ $key } );
	for my $each1(0 .. $#split_names) {
		if(length($split_names[$each1])>1) {
			my $id = $user_IDs_key{$split_names[$each1]};
			print OUT5 "$id\t";
		}
	}
}

close(OUT5);

print "script completed\n";
exit;

sub score_sequence_pair {

	my $sub_input = shift;
	my @split_sub_input = split(/___/, $sub_input);
	my $comparison_seq_copy = $split_sub_input[0];
	my $current_seq_copy = $split_sub_input[1];

	if(length($comparison_seq_copy) <= 1 || length($current_seq_copy)<=1) {
		print "comp:" , length($comparison_seq_copy) , " curr:" ,length($current_seq_copy) , "\n";
		die "ERROR no sequence input to sub score_sequence_pair\n";
	}

	if(length($comparison_seq_copy) != length($current_seq_copy)) {
		print "comp:" , length($comparison_seq_copy) , " curr:" ,length($current_seq_copy) , "\n";
		print "id's:$current_seq_id vs $comparison_seq_id\n";
		die "ERROR sequence length mismatch. check they are aligned, with no weird characters in seq\n";
	}

	my @standard_ambiguity_characters_array = keys %standard_ambiguity_characters;

	foreach my $character(@standard_ambiguity_characters_array) {
		my $alternative_characters = $standard_ambiguity_characters{$character};

		#### loop through each ambiguity character , N, R, Y, etc. 
		# for each char, find any example in one of the sequences. If found, check the corresponding position in the other sequence of the pair and see if its also ambiguous all sites with an ambiguous base in either one or both of the sequence pair, are replaced by the character Z (a character not used in DNA seqs) then check whether the sequence pair is identical after any such replacements made

		my $test_char;
		my $regex_match_position;
		while($comparison_seq_copy =~ /$character/) {
			$regex_match_position = $-[0];
			$test_char=substr($current_seq_copy,$regex_match_position,1);
			if ($test_char=~ /[$alternative_characters]/){
				$curr_score++
			}
			substr($comparison_seq_copy,$regex_match_position,1)="Z";
			substr($current_seq_copy,$regex_match_position,1)="Z";
		}

		while($current_seq_copy=~/$character/) {
			$regex_match_position = $-[0];
			$test_char=substr($comparison_seq_copy,$regex_match_position,1);
			if ($test_char=~/[$alternative_characters]/i){
				$comp_score++
			}
			substr($comparison_seq_copy,$regex_match_position,1)="Z";
			substr($current_seq_copy,$regex_match_position,1)="Z";
		}
	}
	return( $comparison_seq_copy . "___" . $current_seq_copy );
}