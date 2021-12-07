#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';
use Getopt::Long;

my $docstring=<<END_DOC
runDeepFoldPipeline.pl -datadir=<data directory> -generate_msa=<msa flag> -run_refine=<refine flag> -nr_path=<path to the nr database> -uniclust_path=<path to the uniclust database> -uniref_path=<path to the uniref database> -metaclust_path=<path to the metaclust database>
    run the full DeepFold pipeline:
    [1] scripts/1_predPhiPsi.pl
    [2] scripts/2_runDeepMSA.pl
    [3] scripts/3_runDeepPotential.pl
    [4] scripts/4_runDeepFold.pl
END_DOC
;

####  parse command line argument ####
if (@ARGV<1)
{
    print $docstring;
    exit();
}

my $bindir=dirname(abs_path(__FILE__)); #where script and programs are
GetOptions( 'datadir=s' => \my $datadir         # where the input sequence is 
          , 'generate_msa=s' => \my $msa_flag   # True/False generate MSA by DeepMSA
          , 'run_refine=s' => \my $refine_flag  # True/False run refinement by ModRefiner
	  , 'nr_path=s' => \my $nr_path
	  , 'uniclust_path=s' => \my $uniclust_path
	  , 'uniref_path=s' => \my $uniref_path
	  , 'metaclust_path=s'  => \my $metaclust_path
          );

if($datadir eq ""){
    print $docstring;
    exit();
}
if($msa_flag eq ""){
    $msa_flag="True";
}
if($refine_flag eq ""){
    $refine_flag="False";
}

#### check input file ####
if (!-s "$datadir/seq.txt")
{
    die "FATAL ERROR! Cannot find input sequence $datadir/seq.txt";
}

#### [1] 1_predPhiPsi.pl:       predicted phi/psi angles and secondary structure
#### [2] 2_runDeepMSA.pl:       generate MSA
#### [3] 3_runDeepPotential.pl: predict spatial restraints
#### [4] 4_runDeepFold.pl:      run folding simulations and refinement
while(1)
{
    system("$bindir/scripts/1_predPhiPsi.pl        $datadir $nr_path");
    if($msa_flag eq "True"){
        system("$bindir/scripts/2_runDeepMSA.pl    $datadir $uniclust_path $uniref_path $metaclust_path");
    }
    system("$bindir/scripts/3_runDeepPotential.pl  $datadir");
    system("$bindir/scripts/4_runDeepFold.pl       $datadir $refine_flag");
    last if (-s "$datadir/full_atom_model.pdb");
}

exit(0);
