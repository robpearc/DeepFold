#!/usr/bin/perl -w

###################################################################################################
# convert_PDB.pl PDB.old PDB.new
#
# A script used to convert an old PDB file into a standard PDB file format. The heavy atoms will be 
# written in the order of N->CA->C->O->[side chain atoms]->[OXT if applicable], followed by the 
# hydrogen atoms.
#
# Written by Xiaoqiang Huang (xiaoqiah@umich.edu)
# Last updated: 05/02/2020
##################################################################################################

use strict;
use warnings;

if(@ARGV!=2){
  printf "usage: convert_PDB.pl PDB.old PDB.new\n";
  printf "note that only the order of atom in each residue will be modified.\n";
  exit(0);
}

# build a hash-of-hash table to assign the atoms and their indices in each residue type
# please note that only 20 canonical residues are supported
my %HoH=(
"GLY"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA1'=>30,' HA2'=>31,' HA3'=>32,' HA '=>30},
"ALA"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42},
"SER"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' OG '=> 5,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG '=>50,' HG1'=>50},
"CYS"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' SG '=> 5,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG '=>50,' HG1'=>50},
"THR"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' OG1'=> 5,' CG2'=> 6,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB '=>40,
        ' HG1'=>50,'1HG2'=>51,'2HG2'=>52,'3HG2'=>53},
"VAL"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG1'=> 5,' CG2'=> 6,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB '=>40,
        '1HG1'=>50,'2HG1'=>51,'3HG1'=>52,'1HG2'=>53,'2HG2'=>54,'3HG2'=>55},
"ARG"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD '=> 6,' NE '=> 7,' CZ '=> 8,' NH1'=> 9,' NH2'=>10,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG1'=>50,' HG2'=>51,' HG3'=>52,
        ' HD1'=>60,' HD2'=>61,' HD3'=>62,
        ' HE '=>70,
        '1HH1'=>80,'2HH1'=>81,'1HH2'=>82,'2HH2'=>83},
"ASP"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' OD1'=> 6,' OD2'=> 7,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42},
"ASN"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' OD1'=> 6,' ND2'=> 7,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        '1HD2'=>50,'2HD2'=>51},
"GLN"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD '=> 6,' OE1'=> 7,' NE2'=> 8,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG1'=>50,' HG2'=>51,' HG3'=>52,
        '1HE2'=>60,'2HE2'=>61},
"GLU"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD '=> 6,' OE1'=> 7,' OE2'=> 8,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG1'=>50,' HG2'=>51,' HG3'=>52},
"HIS"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' ND1'=> 6,' CD2'=> 7,' CE1'=> 8,' NE2'=> 9,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HD1'=>50,' HD2'=>51,
        ' HE1'=>60,' HE2'=>61},
"ILE"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG1'=> 5,' CG2'=> 6,' CD1'=> 7,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB '=>40,
        '1HG1'=>50,'2HG1'=>51,'3HG1'=>52,'1HG2'=>53,'2HG2'=>54,'3HG2'=>55,
        '1HD1'=>60,'2HD1'=>61,'3HD1'=>62,' HD1'=>63,' HD2'=>64,' HD3'=>65},
"LYS"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD '=> 6,' CE '=> 7,' NZ '=> 8,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG1'=>50,' HG2'=>51,' HG3'=>52,
        ' HD1'=>60,' HD2'=>61,' HD3'=>62,
        ' HE1'=>70,' HE2'=>71,' HE3'=>72,
        ' HZ1'=>80,' HZ2'=>81,' HZ3'=>82},
"LEU"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD1'=> 6,' CD2'=> 7,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG '=>50,
        '1HD1'=>60,'2HD1'=>61,'3HD1'=>62,'1HD2'=>63,'2HD2'=>64,'3HD2'=>65},
"MET"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' SD '=> 6,' CE '=> 7,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG1'=>50,' HG2'=>51,' HG3'=>52,
        ' HE1'=>60,' HE2'=>61,' HE3'=>62},
"PRO"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD '=> 6,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HG1'=>50,' HG2'=>51,' HG3'=>52,
        ' HD1'=>60,' HD2'=>61,' HD3'=>62},
"PHE"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD1'=> 6,' CD2'=> 7,' CE1'=> 8,' CE2'=> 9,' CZ '=>10,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HD1'=>50,' HD2'=>51,' HD3'=>52,
        ' HE1'=>60,' HE2'=>61,' HE3'=>62,
	' HZ '=>70},
"TRP"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD1'=> 6,' CD2'=> 7,' NE1'=> 8,' CE2'=> 9,' CE3'=>10,' CZ2'=>11,' CZ3'=>12,' CH2'=>13,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HD1'=>50,
        ' HE1'=>60,' HE2'=>61,' HE3'=>62,
	' HZ1'=>70,' HZ2'=>71,' HZ3'=>72,
	' HH2'=>80},
"TYR"=>{' N  '=> 0,' CA '=> 1,' C  '=> 2,' O  '=> 3,' CB '=> 4,' CG '=> 5,' CD1'=> 6,' CD2'=> 7,' CE1'=> 8,' CE2'=> 9,' CZ '=>10,' OH '=>11,' OXT'=>15,
        ' HT1'=>20,' HT2'=>21,' HT3'=>22,' HN1'=>20,' HN2'=>21,' HN3'=>22,' H1 '=>20,' H2 '=>21,' H3 '=>22,' H  '=>20,' HN '=>20,
        ' HA '=>30,
        ' HB1'=>40,' HB2'=>41,' HB3'=>42,
        ' HD1'=>50,' HD2'=>51,
        ' HE1'=>60,' HE2'=>61,
	' HH '=>70},
);

#----------------------------------------------------------
# beginning of the script
#----------------------------------------------------------
my $pdbold=$ARGV[0];
my $pdbnew=$ARGV[1];

my $NATOM=100;

open FILE,"<$pdbold";
open FILE2,">$pdbnew";
my $iniChId='UNK';
my $iniRsId='UNKNOWN';
my $firstRs=1;
my $atCounter=0;
my @alllines=();
my (@tails)=();
my $rsDetected=0;
for(my $i=0;$i<$NATOM;$i++){push @tails,'';}
while(<FILE>){
  chomp(my $line=$_);
  my $key=substr $line,0,4;
  if($key ne 'ATOM'){
    #if($key eq 'TER '){$atCounter++;}
    if($rsDetected==0){
      printf FILE2 "$line\n";
      next;
    }
    # output tails in the last residue
    $rsDetected=0;
    for(my $i=0;$i<$NATOM;$i++){
      if($tails[$i] ne ''){
        $atCounter++;
        printf FILE2 "ATOM  ".sprintf("%5d ",$atCounter)."$tails[$i]\n";
      }
    }
    printf FILE2 "$line\n";
    # clear the atom array
    for(my $i=0;$i<$NATOM;$i++){$tails[$i]='';}
    next;
  }
  $rsDetected=1;
  my $chId=substr $line,21,1;
  my $rsId=substr $line,22,4;
  my $rsName=substr $line,17,3;
  my $atName=substr $line,12,4;
  # change histidine name
  $atName='HIS' if($rsName eq 'HSD' or $rsName eq 'HSE');
  # change C-ter OT1 and OT2 into O and OXT, respectively
  $atName=' O  ' if($atName eq ' OT1');
  $atName=' OXT' if($atName eq ' OT2');
  my @buff=split(//,$atName);
  if(@buff!=4){
    printf("ERROR the length of atom name field is greater than 4, skip\n");
    #printf(FILE2 "ERROR the length of atom name field is greater than 4, skip\n");
  }
  if($buff[0] eq 'H'){$atName=$buff[3].$buff[0].$buff[1].$buff[2];}
  elsif($buff[0]=~/\d/){
    for(my $k=2;$k<4;$k++){
      if($buff[$k] eq ' '){
        my $temp=$buff[0];
        $buff[0]=$buff[$k];
	$buff[$k]=$temp;
	last;
      }
    }
    $atName=$buff[0].$buff[1].$buff[2].$buff[3];
  }
  # encounter a new chain, therefore must encounter a new residue
  if($iniChId ne $chId){
    $iniChId=$chId;
    $iniRsId = $rsId;
    for(my $i=0;$i<$NATOM;$i++){
      if($tails[$i] ne ''){
        $atCounter++;
	printf FILE2 "ATOM  ".sprintf("%5d ",$atCounter)."$tails[$i]\n";
      }
    }
    # clear the atom array
    for(my $i=0;$i<$NATOM;$i++){$tails[$i]='';}
  }
  # encounter a new residue
  elsif($iniRsId ne $rsId){
    $iniRsId = $rsId;
    for(my $i=0;$i<$NATOM;$i++){
      if($tails[$i] ne ''){
        $atCounter++;
	printf FILE2 "ATOM  ".sprintf("%5d ",$atCounter)."$tails[$i]\n";
      }
    }
    # clear the atom array
    for(my $i=0;$i<$NATOM;$i++){$tails[$i]='';}
  }
  if(!exists($HoH{$rsName}{$atName})){
    printf "WARNING residue $rsName$rsId does not have atom $atName, skip!\n";
    #printf FILE2 "WARNING residue $rsName$rsId does not have atom $atName, skip!\n";
    next;
  }
  if($tails[$HoH{$rsName}{$atName}] eq ''){
    $tails[$HoH{$rsName}{$atName}]=$atName.substr($line,16);
  }
}
# output the last residue in case there is no 'TER' 
for(my $i=0;$i<$NATOM;$i++){
  if($tails[$i] ne ''){
    $atCounter++;
    printf FILE2 "ATOM  ".sprintf("%5d ",$atCounter)."$tails[$i]\n";
  }
}
# clear the atom array
for(my $i=0;$i<$NATOM;$i++){$tails[$i]='';}
close FILE;
close FILE2;

exit;
