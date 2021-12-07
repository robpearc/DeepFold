#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $docstring="DeepPotentialmod.pl datadir tag";
if (@ARGV!=2)
{
    print "$docstring\n";
    exit();
}

my $datadir    ="$ARGV[0]";
my $tag        ="$ARGV[1]"; # for tmp folder
my $distancedir=dirname(abs_path(__FILE__));
my $aln        ="$datadir/MSA/protein.aln";

################ working directory ########################
my $work_dir="/tmp/$ENV{USER}/$tag";
system("rm -fr $work_dir");
system("mkdir -p $work_dir");
chdir "$work_dir";

#### DeepPotential needs a3m instead of aln ####
my @suffix_array=qw(
        a3m
        hmsa3m
        hh3a3m
        jaca3m
        hhba3m
        );
my $a3m="";
foreach my $suffix (@suffix_array)
{
    $a3m="$aln"; 
    $a3m=~s/.aln$/.$suffix/mg;
    if (-s "$a3m")
    {
        system("cp $a3m protein.a3m");
        system("head -2 protein.a3m > seq.fasta");
        last;
    }
    elsif (-s "$a3m.gz")
    {
        system("cp $a3m.gz protein.a3m.gz");
        system("gzip -d protein.a3m.gz");
        system("head -2 protein.a3m > seq.fasta");
        last;
    }
}

if ( ! -s "protein.a3m" )
{
    printf "warning: protein.a3m does not exist! use protein.aln instead\n";
    system("cp $aln protein.aln");

    if (! -s "protein.aln")
    {
        printf "error: $aln does not exist!\n";
        exit(1);
    }
    my $sequence=`head -1 protein.aln`;
    chomp($sequence);
    open(FP,">seq.fasta");
    printf FP ">protein\n$sequence\n";
    close(FP);
}

printf "---- run DeepPotential ----\n";
if (-s "$datadir/DeepPotential_20.npz.gz" && -s "$datadir/DeepPotential_pca_20.txt.gz")
{
    #system("cp $datadir/DeepPotential_20.npz.gz .");
    #system("cp $datadir/DeepPotential_pca_20.txt.gz .");
    #system("cp $datadir/DeepPotential_ca_contact.txt.gz .");
    #system("cp $datadir/DeepPotential_cb_contact.txt.gz .");
    #system("gzip -d     DeepPotential_20.npz.gz");
    #system("gzip -d     DeepPotential_pca_20.txt.gz");
    #system("gzip -d     DeepPotential_ca_contact.txt.gz");
    #system("gzip -d     DeepPotential_cb_contact.txt.gz");
    #system("gunzip      DeepPotential_20.npz.gz");
    #system("gunzip      DeepPotential_pca_20.txt.gz");
    #system("gunzip      DeepPotential_ca_contact.txt.gz");
    #system("gunzip      DeepPotential_cb_contact.txt.gz");
}
else
{
    my $cmd="python $distancedir/DeepPotential_04132020/predict.py";
    if ( -s "protein.a3m" )
    {
        system("$cmd protein.a3m DeepPotential");
    }
    else
    {
        system("$cmd protein.aln DeepPotential");
    }
   
    foreach my $target (`ls |grep _`)
    {
        chomp($target);
        system("cat $target|gzip - > $target.gz");
        system("cp $target.gz $datadir/");
    }
}

################# ending procedure ######################
system("sync");
sleep(1);
system("rm -rf $work_dir");
exit(0);
