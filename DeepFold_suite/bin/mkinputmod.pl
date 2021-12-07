#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $docstring="mkinputmod tag datadir nr_path";
if (@ARGV!=3)
{
    print "$docstring\n";
    exit();
}

# input files of this program
#     seq.txt
#     rmsinp  (it will be generated if without it)
#     seq.dat (it will be generated if without it)
#     exp.dat (it will be generated if without it)
#     

my $tag      ="$ARGV[0]"; # full job name, for tmp folder
my $datadir  ="$ARGV[1]"; # input data folder
my $nr_path  ="$ARGV[2]"; # path to the nr database
my $recorddir="$datadir/record";
my $bindir   =dirname(abs_path(__FILE__)); #folder of main QUARK program

################# directories #############################
my $work_dir="/tmp/$ENV{USER}/$tag";
system("/bin/rm -fr $work_dir");
system("/bin/mkdir -p $work_dir");
chdir "$work_dir";

################# run jobs ###############################
# put your stuff here ------------------>

my %ts=(
     'GLY'=>'G',
     'ALA'=>'A',
     'VAL'=>'V',
     'LEU'=>'L',
     'ILE'=>'I',
     'SER'=>'S',
     'THR'=>'T',
     'CYS'=>'C',
     'MET'=>'M',
     'PRO'=>'P',
     'ASP'=>'D',
     'ASN'=>'N',
     'GLU'=>'E',
     'GLN'=>'Q',
     'LYS'=>'K',
     'ARG'=>'R',
     'HIS'=>'H',
     'PHE'=>'F',
     'TYR'=>'Y',
     'TRP'=>'W',

     'ASX'=>'B',
     'GLX'=>'Z',
     'UNK'=>'X',

     'G'=>'GLY',
     'A'=>'ALA',
     'V'=>'VAL',
     'L'=>'LEU',
     'I'=>'ILE',
     'S'=>'SER',
     'T'=>'THR',
     'C'=>'CYS',
     'M'=>'MET',
     'P'=>'PRO',
     'D'=>'ASP',
     'N'=>'ASN',
     'E'=>'GLU',
     'Q'=>'GLN',
     'K'=>'LYS',
     'R'=>'ARG',
     'H'=>'HIS',
     'F'=>'PHE',
     'Y'=>'TYR',
     'W'=>'TRP',

     'a'=>'CYS',
     'b'=>'CYS',
     'c'=>'CYS',
     'd'=>'CYS',
     'e'=>'CYS',
     'f'=>'CYS',
     'g'=>'CYS',
     'h'=>'CYS',
     'i'=>'CYS',
     'j'=>'CYS',
     'k'=>'CYS',
     'l'=>'CYS',
     'm'=>'CYS',
     'n'=>'CYS',
     'o'=>'CYS',
     'p'=>'CYS',
     'q'=>'CYS',
     'r'=>'CYS',
     's'=>'CYS',
     't'=>'CYS',
     'u'=>'CYS',
     'v'=>'CYS',
     'w'=>'CYS',
     'x'=>'CYS',
     'y'=>'CYS',
     'z'=>'CYS',

     'B'=>'ASX',
     'Z'=>'GLX',
     'X'=>'CYS',
     );

my $bindir_hh="$bindir/hhsuite2/bin"; # QUARK programs
my $db="$nr_path/nr";
my $psipreddir="$bindir/psipred"; #need to move to local


######## copy needed files:
system("cp $datadir/seq.txt $work_dir");
system("cp $datadir/seq.fasta $work_dir/seq.txt") if (!-s "$datadir/seq.txt");

######## prepare 'seq.fasta':
my $sequence="";
foreach my $seqtxt(`cat $datadir/seq.txt`)
{
    next if($seqtxt=~/\>/);
    $seqtxt=~s/\s//mg;
    $seqtxt=~s/\n//mg;
    $sequence=$sequence.$seqtxt;
}
my $Lch=length $sequence;
open(seq,">seq.fasta");
printf seq ">protein\n";
my %seqQ3;
for(my $i=1;$i<=$Lch;$i++){
    my $a=substr($sequence,$i-1,1);
    printf seq "$a";
    $seqQ3{$i}=$ts{$a};
    printf seq "\n" if($i==int($i/80)*80)
}
printf seq "\n";
close(seq);


if(!-s "$datadir/psitmp.chk" || !-s "$datadir/pssm.txt" || !-s "$datadir/blast.out" || !-s "$datadir/mtx")
{
    print "running blast ...\n";
    if(1)
    {
        system("$bindir_hh/blastpgp  -b 1000 -j 3 -h 0.001 -d $db -i seq.fasta -C psitmp.chk -Q pssm.txt > blast.out");
        system("cp psitmp.chk $datadir");
        system("cp pssm.txt $datadir");
        system("cp blast.out $datadir");
        system("cp seq.fasta psitmp.fasta");
        system("echo psitmp.chk > psitmp.pn");
        system("echo psitmp.fasta > psitmp.sn");
        system("$bindir_hh/makemat -P psitmp");
        system("mv psitmp.mtx mtx");
        system("cp mtx $datadir");
    }
}
else
{
    system("cp $datadir/psitmp.chk .");
    system("cp $datadir/blast.out .");
    system("cp $datadir/pssm.txt .");
    system("cp $datadir/mtx .");
}

######## prepare 'protein.ss2':
if(!-s "$datadir/protein.ss2" || !-s "$datadir/protein.horiz")
{
    print "doing psipred\n";
    #`$psipreddir/psipred mtx $psipreddir/weights.dat $psipreddir/weights.dat2 $psipreddir/weights.dat3 $psipreddir/weights.dat4 > protein.ss`;
    system("$psipreddir/bin/psipred mtx $psipreddir/data/weights.dat $psipreddir/data/weights.dat2 $psipreddir/data/weights.dat3 > protein.ss");
    #`$psipreddir/psipass2 $psipreddir/weights_p2.dat 1 1.0 1.0 protein.ss2 protein.ss > protein.horiz`;
    system("$psipreddir/bin/psipass2 $psipreddir/data/weights_p2.dat 1 1.0 1.0 protein.ss2 protein.ss > protein.horiz");
    system("cp protein.ss2 $datadir");
    system("cp protein.horiz $datadir");
}
else{
    system("cp $datadir/protein.ss2 .");
    system("cp $datadir/protein.horiz .");
}
if(!-s "$datadir/seq.dat.ss"){
    ##seq.dat.ss is for simulation
    system("cp $datadir/protein.ss2 $datadir/seq.dat.ss");
}

######## prepare 'seq.dat':
my %sec;
if(!-s "$datadir/seq.dat"){
    print "horiz to seq.dat\n";
    open(psipred,"protein.horiz");  
    open(seq,">$datadir/seq.dat");
    my $j=0;
    while(my $line=<psipred>)
    {
        if($line=~/Conf:\s+(\d+)/)
        {
            my $conf=$1;
            <psipred>=~/Pred:\s+(\S+)/;
            my $pred=$1;
            <psipred>=~/AA:\s+(\S+)/;
            my $aa=$1;
            my $num=length $aa;
            for(my $i=1;$i<=$num;$i++)
            {
                $j++;
                my $conf1=substr($conf,$i-1,1);
                my $pred1=substr($pred,$i-1,1);
                my $aa1=substr($aa,$i-1,1);
                $sec{$j}=1;
                $sec{$j}=2 if($conf1 >=1 && $pred1 eq 'H');
                $sec{$j}=4 if($conf1 >=1 && $pred1 eq 'E');
                printf seq "%5d   %3s%5d%5d\n",$j,$seqQ3{$j},$sec{$j},$conf1;
            }
        }
    }
    close(seq);
    close(psipred);
}
`/bin/cp $datadir/seq.dat .`;

if(!-s "$datadir/sol.txt")
{
    print "predict solve\n";
    system("$bindir/getannfeature 10 psitmp.chk protein.ss2 annfeat10.dat");
    system("$bindir/simple_onetestsas 10 60 $bindir/sastrainres10-60.net annfeat10.dat sol.txt");
    system("cp $work_dir/sol.txt $datadir");
}else{
    system("cp $datadir/sol.txt .");
}

if(!-s "$datadir/phi.txt" || !-s "$datadir/psi.txt")
{
    print "predict phipsi\n";
    #input seq; output predicted phi psi; 
    system("$bindir/getannfeature 10 psitmp.chk protein.ss2 annfeat10.dat");
    system("$bindir/simple_onetestphi 10 50 $bindir/phitrainres10-50.net annfeat10.dat phi.txt");#convert
    system("$bindir/getannfeature 8 psitmp.chk protein.ss2 annfeat8.dat");
    system("$bindir/simple_onetestpsi 8 80 $bindir/psitrainres8-80.net annfeat8.dat psi.txt");
    system("cp phi.txt $datadir");
    system("cp psi.txt $datadir");
}
else
{
    system("cp $datadir/phi.txt .");
    system("cp $datadir/psi.txt .");
}

################# ending procedure ######################
system("sync");
sleep(1);
system("rm -fr $work_dir");

exit();
