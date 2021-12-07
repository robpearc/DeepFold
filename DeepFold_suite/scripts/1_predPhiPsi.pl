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

######## What this program does? ###############################
#
# This program generates input files for DeepFold.
#   input files:
#       seq.txt         (query sequence in FASTA format)
#   output files:
#       seq.dat.ss      (predicted secondary structure)
#       phi.txt         (predicted phi angles)
#       psi.txt         (predicted psi angles)
################################################################


################################################################
my @ss;
my $outdir;
my $CQ_dir=dirname(abs_path(__FILE__)); #location of this script 
my $bindir="$CQ_dir/../bin"; #where DeepFold scripts and programs are
my $nr_path;

#### The argument to the program is the data directory ####
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
    $nr_path="$ARGV[1]";
}
$outdir=abs_path($outdir);

#### command line arguments ended ####
my @ff=qw(
    seq.dat.ss
    phi.txt
    psi.txt
    );

foreach my $s(@ss)
{
    my $datadir="$outdir/$s";
    my $recorddir="$datadir/record"; #to record all log files
    system("mkdir -p $recorddir")   if (!-d "$recorddir");

    #### Prepare input ####
    if (! -s "$datadir/seq.txt")
    {
        if (-s "$datadir/seq.fasta")
        {
            system("cp $datadir/seq.fasta $datadir/seq.txt");
        }
        else
        {
            print "ERROR! $datadir/seq.txt missing.\n";
            next;
        }
    }

    my $inputFasta="$datadir/seq.txt";
    if (!-s "$inputFasta")
    {
        $inputFasta="$datadir/seq.fasta";
        if (!-s "$datadir/seq.fasta")
        {
            print "ERROR! $datadir/seq.txt missing.\n";
            next;
        }
        system("cp $datadir/seq.fasta $inputFasta");
    }
    if (`cat $inputFasta`!~/^>/)
    {
        print "WARNING! rewriting $inputFasta due to missing sequence header\n";
        my $txt=">$s\n";
        foreach my $line(`cat $inputFasta`)
        {
            chomp($line);
            $txt.="$line" if ($line!~/>/);
        }
        open(FP,">$inputFasta");
        print FP "$txt\n";
        close(FP);
    }

    ##### predict ss, phi/psi angles ####
    my $tag="mkin_$s";
    my $flag="False";
    foreach my $f(@ff)
    {
        if (! -s "$datadir/$f")
        {
            $flag="True";
            last;
        }
    }

    if($flag eq "True")
    {
        printf "$tag is missing ----\n";
        my $jobname="$recorddir/$tag";
        open(FP,">$jobname");
        print FP<<EOF
#!/bin/bash
#SBATCH -e $jobname.err
#SBATCH -o $jobname.out
#SBATCH -t 24:00:00
#SBATCH --mem=3500mb
#SBATCH --job-name=$tag

echo hostname: `hostname`   >$recorddir/ware_$tag
echo starting time: `date` >>$recorddir/ware_$tag
echo pwd `pwd`             >>$recorddir/ware_$tag 

$bindir/mkinputmod.pl $tag $datadir $nr_path

echo ending time: `date`   >>$recorddir/ware_$tag
EOF
;
        close(FP);
        system("chmod a+x $jobname");
        system("$jobname");
        #&submitjob($jobname,"$recorddir/note.txt");
    }

}
exit();


sub submitjob
{
    ## jobname is the SBATCH script, note is the file recording the job submission
    my ($jobname,$note) = @_;
    
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
    open(FP,">>$note");
    print "$jobname\t at $date $bsub\n";
    print FP "$jobname\t at $date $bsub\n";
    close(FP);
}

