#!/usr/bin/perl
$datadir="$ARGV[0]"; 
#Only CA coordinates in $fixname are useful. The length of $fixname can be shorter than the complete length, which will make the refinement of the missing region more flexible.
$libdir="$ARGV[1]";#executable, library file
$refname="$ARGV[2]";#model to be refined
$fixname="$ARGV[3]";#reference model

$rannum=0;#initial random number [0,9999]. Modify this parameter only when you want to generate different refined models.
$strength=100;#strength value in [0,100]. Make this parameter smaller if you want the refined model to be guided more by full-atomic energy terms.;

chdir $datadir;

system("$libdir/mcrefinement $datadir $libdir $refname $fixname $rannum");
system("$libdir/emrefinement $datadir $libdir mc$refname $fixname $strength $rannum");
#system("$libdir/HAAD ./emmc$refname");
#system("$libdir/addchainid ./emmc$refname\.h ./hemmc$refname");
my @rst=`cat emmc$refname`;
open(OUT, ">hemmc$refname");
foreach my $r(@rst)
{
    if($r =~ /^ATOM/ || $r =~ /^TER/)
    {
	print OUT "$r";
    }
}
close(OUT);


exit();
