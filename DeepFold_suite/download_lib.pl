#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $libdir    = "";


GetOptions('libdir:s' => \$libdir);


my $fname=""; # in case user provides a file name
($libdir, $fname)=&get_absolute_path($libdir);


if(!$libdir)
{
    print "\nPlease provide a valid path for the DeepFold sequence database\n\n";
    print "Usage:\n";
    print "./download_lib.pl -libdir libdir\n";    
    printf "====================\n";
    printf "Mandatory arguments:\n";
    printf "====================\n";
    printf "-libdir template_library_directory (full path for saving the library files, such as /home/yourname/ITLIB)\n\n";

    #&get_library_size();
    exit();
}

#&get_library_size();


printf "\nYour settings for the DeepFold library are:\n\n";
printf "-libdir   = $libdir\n";

if(!-d "$libdir")
{
    system("mkdir -p $libdir") == 0 or die "System call failed: $!";
}

if(!-d "$libdir")
{
    print "Failed to create $libdir. Please check the problem with your system\n";
    exit;
}

chdir $libdir;

##NR
if(1)
{
    my $nfiles=`ls $libdir/nr/ |wc -l`;
    chomp($nfiles);

    if($nfiles<30) #avoid re-download
    {
	print "Downloading NR database...\n";
	system("wget -o log -c http://zhanglab.ccmb.med.umich.edu/library/nr.tar.gz") == 0 or die "System call failed: $!";  
	system("tar -xzvf nr.tar.gz")== 0 or die "System call failed: $!";
	system("rm -f nr.tar.gz") == 0 or die "System call failed: $!";
	system("rm -f log") == 0 or die "System call failed: $!";
    }
    else
    {
	print "NR database is already present. Skipping nr download ..\n";  
    }
}

if(1)
{
    
    chdir $libdir;

    print "Downloading DeepMSA library files...\n";

    system("wget -o log -c http://gwdu111.gwdg.de/~compbiol/uniclust/2017_04/uniclust30_2017_04_hhsuite.tar.gz") == 0 or die "System call failed: $!";
    system("tar -zxvf uniclust30_2017_04_hhsuite.tar.gz >log") == 0 or die "System call failed: $!";
    system("rm -f uniclust30_2017_04_hhsuite.tar.gz") == 0 or die "System call failed: $!";
    system("rm -f log") == 0 or die "System call failed: $!";

    system("wget -o log -c ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz") == 0 or die "System call failed: $!";
    system("gzip -d uniref90.fasta.gz >log") == 0 or die "System call failed: $!";
    system("rm -f uniref90.tar.bz2") == 0 or die "System call failed: $!";
    system("rm -f log") == 0 or die "System call failed: $!";

    system("wget -o log -c https://metaclust.mmseqs.org/2017_05/metaclust_2017_05.fasta.gz") == 0 or die "System call failed: $!";
    system("tar -zxvf metaclust_2017_05.fasta.gz >log") == 0 or die "System call failed: $!";
    system("rm -f metaclust_2017_05.fasta.gz") == 0 or die "System call failed: $!";
    system("rm -f log") == 0 or die "System call failed: $!";
}

print "All libraries were downloaded\n";



sub get_library_size
{ 
    my $total=0;
    my $total1=0;

    system("rm -f index.html") == 0 or die "System call failed: $!";
    system("wget -o log http://zhanglab.ccmb.med.umich.edu/library/") == 0 or die "System call failed: $!";
    my $rst=`cat index.html`;

    my @files=(
	       "PDB",
	       "MTX",
	       "DEP",
	       "SIG",
	       "CNT",
	       "dotProfiles",       	       
	       "Bfactor",	       
	       );


    my %size=(
	      "nr"=>7.3,
	      "PDB"=>4.4, 
	      "MTX"=>1.8, 
	      "DEP"=>2.8, 
	      "SIG"=>1.9, 
	      "CNT"=>0.87, 
	      "dotProfiles"=>0.85, 
	      "Bfactor"=>0.95,  
	      "summary"=>2.6, 	      
	      );



    my @keys = keys %size;
    foreach my $k (@keys)
    {
	$total1 += $size{$k};
    }


    printf "\nThe C-QUARK library files to be downloaded include:\n";
    printf "--------------------------------------------\n";
    printf "%-15s     Size          Size\n", "Name";
    printf "%-15s (compressed) (decompressed)\n", "";

    printf "--------------------------------------------\n";
    printf "%-15s %-5.2f GB       %5.2f GB\n", "nr", 4198/1024, $size{"nr"}; 
    $total +=4198; 

    foreach my $f(@files)
    {
	if($rst =~ /$f\.tar\.bz2.+\s+\(~(\d+)MB\)/)
	{
	    printf "%-15s %-5.2f GB       %5.2f GB\n", $f, $1/1024, $size{$f};
	    $total += $1;
	}
    }

    printf "%-15s %-5.2f GB       %5.2f GB\n", "summary", 226/1024, $size{"summary"};
    $total +=226;
    printf "%-15s %-5.2f GB       %5.2f GB\n", "others", 1, 6.5;

    $total=11;
    $total1=50;
    printf "\n%-15s~%-5.2f GB      ~%5.2f GB\n", "TOTAL", $total, $total1;
    printf "--------------------------------------------\n";

    system("rm -f index.html log") == 0 or die "System call failed: $!";
}


sub get_absolute_path
{
    my ($path)=@_;
    my $dir = "";
    my $fname= "";
    if($path ne "")
    {
	`mkdir -p $path` if(!-e $path && !-d $path);
	my $apath = abs_path($path); # for '..'	
	if (defined($apath))
	{
	    $dir = $apath; 
	    if(-e $dir && !-d $dir) #it is not an existing directory, should be a file
	    {
		$dir =~ s/(.*)\/.*$/$1/;
		$fname=$apath;
		$fname =~ s/^.*[\/\\]//;
	    }	
	}
    }
    return ($dir, $fname);
}
