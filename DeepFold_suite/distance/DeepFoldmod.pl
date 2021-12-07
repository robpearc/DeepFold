#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $docstring="DeepPotentialmod.pl datadir tag refine_flag";
if (@ARGV!=3)
{
    print "$docstring\n";
    exit();
}

my $datadir     = "$ARGV[0]";
my $tag         = "$ARGV[1]"; # for tmp folder
my $refine_flag = "$ARGV[2]";
my $distancedir=dirname(abs_path(__FILE__));

################ working directory ########################
my $work_dir="/tmp/$ENV{USER}/$tag";
system("rm -fr $work_dir");
system("mkdir -p $work_dir");
chdir "$work_dir";

printf "---- check DeepPotential restraints ----\n";
if (-s "$datadir/DeepPotential_20.npz.gz" && -s "$datadir/DeepPotential_pca_20.txt.gz")
{
    system("cp $datadir/DeepPotential_20.npz.gz .");
    system("cp $datadir/DeepPotential_pca_20.txt.gz .");
    system("cp $datadir/DeepPotential_ca_contact.txt.gz .");
    system("cp $datadir/DeepPotential_cb_contact.txt.gz .");
    #system("gzip -d     DeepPotential_20.npz.gz");
    #system("gzip -d     DeepPotential_pca_20.txt.gz");
    #system("gzip -d     DeepPotential_ca_contact.txt.gz");
    #system("gzip -d     DeepPotential_cb_contact.txt.gz");
    system("gunzip      DeepPotential_20.npz.gz");
    system("gunzip      DeepPotential_pca_20.txt.gz");
    system("gunzip      DeepPotential_ca_contact.txt.gz");
    system("gunzip      DeepPotential_cb_contact.txt.gz");
}
else
{
    printf "---- No DeepPotential restraints were generated! Exiting now. ----\n";
    system("sync");
    sleep(1);
    system("rm -rf $work_dir");
    exit(0);
}

if(!-s "$datadir/model.pdb")
{
    printf "---- model structure ----\n";
    system("python $distancedir/../bin/npz2txt.py DeepPotential_20.npz");
    system("cp $distancedir/../DeepFold ./");
    system("cp -r $distancedir/../library ./");
    system("cp $datadir/phi.txt ./");
    system("cp $datadir/psi.txt ./");
    system("cp $datadir/seq.dat.ss ./");
    system("cp $datadir/seq.txt ./");
    system("./DeepFold ./ ./library");
    system("cp model.pdb $datadir/");
}
else
{
    system("cp $datadir/model.pdb ./");
}

my $refmodel="model.pdb";
my $sidmodel="faspr_model.pdb";
my $conmodel="full_atom_model.pdb";
if($refine_flag eq "True" && ! -e "$datadir/hemmc$refmodel"){
    system("$distancedir/../bin/refine/EMrefinement.pl $datadir $distancedir/../bin/abs/mybin $refmodel $refmodel ");
    $refmodel = "hemmc$refmodel";
    system("$distancedir/../bin/FASPR/FASPR -i $refmodel -o $sidmodel");
}
else{
    system("$distancedir/../bin/FASPR/FASPR -i $refmodel -o $sidmodel");
}
system("$distancedir/../bin/convert_PDB.pl $sidmodel $conmodel");
system("cp $conmodel $datadir/");


################# ending procedure ######################
system("sync");
sleep(1);
system("rm -rf $work_dir");
exit(0);
