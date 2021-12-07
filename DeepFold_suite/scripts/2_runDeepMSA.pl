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
# This program generates the DeepMSA MSA
#   input files:
#       seq.txt         (query sequence in FASTA format)
#   output files:
#       MSA/protein.aln     (multiple sequence alignment for restraint prediction)
################################################################


################################################################
my @ss;
my $outdir;
my $CQ_dir=dirname(abs_path(__FILE__)); #location of this script 
my $bindir="$CQ_dir/../bin"; 
my $HHLIB ="$bindir/hhsuite2";
my $uniclust_path;
my $uniref_path;
my $metaclust_path;

######### Needed changes ended #################################
#### command line arguments can overwrite the above variables ####
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
    $uniclust_path="$ARGV[1]";
    $uniref_path="$ARGV[2]";
    $metaclust_path="$ARGV[3]";
}
$outdir=abs_path($outdir);

#### command line arguments ended ####

foreach my $s(@ss)
{
    my $datadir="$outdir/$s";
    my $recorddir="$datadir/record"; #for record all log files
    system("mkdir -p $recorddir")   if (!-d "$recorddir");
    system("mkdir -p $datadir/MSA") if (!-d "$datadir/MSA");

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
    system("cp $inputFasta $datadir/MSA/DeepMSA.fasta");# if (!-s "$datadir/MSA/DeepMSA.fasta");
    system("cp $inputFasta $datadir/MSA/qMSA.fasta");# if (!-s "$datadir/MSA/qMSA.fasta");
    my $Lch=&fasta2len("$inputFasta");
    
    #### Run DeepMSA ####
    my $tag="dMSA_$s"; # unique name
    if (!-s "$datadir/MSA/DeepMSA.aln")
    {
        my $jobname="$recorddir/$tag";
        open(FP,">$jobname");
        print FP<<EOF
#!/bin/bash
#SBATCH -e $jobname.err
#SBATCH -o $jobname.out
#SBATCH -t 48:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=$tag

$HHLIB/scripts/build_MSA.py \\
    -hhblitsdb=$uniclust_path/uniclust30_2017_04 \\
    -jackhmmerdb=$uniref_path/uniref90.fasta \\
    -hmmsearchdb=$metaclust_path/metaclust.fasta \\
    -tmpdir=/tmp/$ENV{USER}/$tag \\
    $datadir/MSA/DeepMSA.fasta
EOF
;
        close(FP);
        system("chmod a+x $jobname");
        system("$jobname");
        #&submitjob($jobname,"$recorddir/note.txt");
    }
    `cp $datadir/MSA/DeepMSA.aln $datadir/MSA/protein.aln`;
}
exit();


sub submitjob
{
    ## jobname is the PBS script, note is the file recording job submission
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

sub fasta2len
{
    my ($filename)=@_;
    my $sequence="";
    foreach my $line(`cat $filename`)
    {   
        if ($line!~/>/)
        {
            chomp($line);
            $sequence.="$line";
        }
    }   
    my $Lch=length $sequence;
    return $Lch; 
}
