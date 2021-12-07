#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

#################################################################
# Disclaimer: DeepFold is the software developed at Y Zhang Lab #
# at CCMB, University of Michigan. No any part of this package  #
# could be released outside the Zhang Lab without permission    #
# from the orginal authors. Violation of this rule may result   #
# in lawful consequences.                                       #
#################################################################

#################################################################
# This script will generate spatial restraints by 
# DeepPotential for DeepFold
#     input:
#       seq.txt          (query sequence in FASTA format)
#       MSA/protein.aln  (multiple sequence alignment)
#     output:
#       DeepPotential_* (spatial restraints) 
#################################################################


#################################################################
my @ss;
my $outdir;
my $distancedir=dirname(abs_path(__FILE__))."/../distance/"; #where distance predictors are

#### command line arguments ####
if (@ARGV>0)
{
    my $datadir=abs_path($ARGV[0]);
    $outdir=dirname($datadir);
    @ss=(basename($datadir));
    if (! -d "$datadir") # target list 
    {
        @ss=();
        foreach my $s(`grep -ohP '^\\S+' $datadir`)
        {
            chomp($s);
            push @ss, $s;
        }
    }
}
$outdir=abs_path($outdir);

#################################################################
my @DF=qw(
    DeepPotential
);

foreach my $s(@ss)
{
    my $datadir="$outdir/$s";

    if(!-s "$datadir/MSA/protein.aln")
    {
        print "Error! $datadir/MSA/protein.aln missing. Skip $s\n";
        next;
    }
    print "$datadir/MSA/protein.aln is available. Let's predict spatial restraints.\n";

    my $recorddir="$datadir/record"; #all log files and intermediate job files
    system("mkdir -p $recorddir");

    my $sequence=`head -1 $datadir/MSA/protein.aln`;
    chomp($sequence);
    my $Lch=length $sequence;

    foreach my $F(@DF)
    {
        my $tag="DeepPotential\_$s"; # unique name
        if(-s "$datadir/$F\_20.npz" || -s "$datadir/$F\_20.npz.gz")
        {
            print "Distance prediction is completed. Let's skip\n";
            next;
        }
        my $walltime="24:00:00";
        my $mem="8000mb";
        if ($Lch>300)
        {
            $walltime="48:00:00";
            $mem="12000mb";
        }
        my $jobname="$recorddir/$tag";
        open(FP,">$jobname");
        print FP<<EOF
#!/bin/bash
#SBATCH -e $jobname.err
#SBATCH -o $jobname.out
#SBATCH -t $walltime
#SBATCH --mem=$mem
#SBATCH --job-name=$tag

echo hostname: `hostname`
echo starting time: `date`
echo pwd `pwd`

$distancedir/${F}mod.pl $datadir $tag

echo ending time: `date`
EOF
;
        close(FP);
        system("chmod a+x $jobname");
        system("$jobname");

        #### record submission status ####
=pod
        my $bsub='';
        while(length $bsub ==0)
        {
            $bsub=`sbatch $jobname`;
            chomp($bsub);
            last if (length $bsub);
            sleep(20);
        }
        my $date=`date`;
        chomp($date);
        open(FP,">>$recorddir/note.txt");
        print "$jobname\t at $date $bsub\n";
        print FP "$jobname\t at $date $bsub\n";
        close(FP);
=cut
    }
}
exit();
