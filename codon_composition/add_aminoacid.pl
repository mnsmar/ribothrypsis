#!/usr/bin/env perl

# Load modules
use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;
use File::Path qw(make_path);

my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Add aminoacid column"],
	[],
	['ifile=s',
		'tab file. Use - for STDIN',
		{required => 1}],
	['codon-col=s',
		'name of codon column'],
	['out-col=s',
		'name of output column'],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

# this is the hash table for the amino acids
my %aacode = (
	TTT => "F", TTC => "F", TTA => "L", TTG => "L",
	TCT => "S", TCC => "S", TCA => "S", TCG => "S",
	TAT => "Y", TAC => "Y", TAA => "STOP", TAG => "STOP",
	TGT => "C", TGC => "C", TGA => "STOP", TGG => "W",
	CTT => "L", CTC => "L", CTA => "L", CTG => "L",
	CCT => "P", CCC => "P", CCA => "P", CCG => "P",
	CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
	CGT => "R", CGC => "R", CGA => "R", CGG => "R",
	ATT => "I", ATC => "I", ATA => "I", ATG => "M",
	ACT => "T", ACC => "T", ACA => "T", ACG => "T",
	AAT => "N", AAC => "N", AAA => "K", AAG => "K",
	AGT => "S", AGC => "S", AGA => "R", AGG => "R",
	GTT => "V", GTC => "V", GTA => "V", GTG => "V",
	GCT => "A", GCC => "A", GCA => "A", GCG => "A",
	GAT => "D", GAC => "D", GAA => "E", GAG => "E",
	GGT => "G", GGC => "G", GGA => "G", GGG => "G",
); 

my $IN = filehandle_for($opt->ifile);
my $header = $IN->getline();
chomp $header;
my @header = split("\t", $header);

my $codoncol;
for (my $i = 0; $i < @header; $i++){
	if ($header[$i] eq $opt->codon_col){
		$codoncol = $i;
	}
}

print $header."\t".$opt->out_col."\n";

while (my $line = $IN->getline) {
	chomp $line;
	my $codon = (split(/\t/, $line))[$codoncol];
	my $aa = $aacode{$codon};
	print $line."\t".$aa."\n";
}


sub filehandle_for {
	my ($file) = @_;

	if ($file eq '-'){
		open(my $IN, "<-");
		return $IN;
	}
	else {
		open(my $IN, "<", $file);
		return $IN
	}
}

