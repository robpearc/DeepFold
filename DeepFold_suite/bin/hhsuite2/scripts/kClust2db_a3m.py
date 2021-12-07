#!/usr/bin/env python
docstring='''
kClust2db.py db.fasta mydb/mydb
    cluster sequences in FASTA file db.fasta using kClust,
    and generate hhblits style database at mydb/mydb

Options:
    -tmpdir=/tmp/$USER/kClust_`date +%N`
        use -tmpdir as temperary folder

    -id=30
        kClust sequence identity cutoff 30%. legal values are: 
        20, 30, 40, 50, 60, 70, 80, 90, 99

    -c=100
        cdhit sequence identity cutoff 100. must be >=40. 
        if set to >100, do not perform cdhit redundancy removal.

    -ncpu=1
        number of CPU threads
'''
import sys,os
import subprocess
import shutil
from string import Template

from HHPaths import bin_dict

id2s_dict= { 20:0.52, 30:1.12, 40:1.73, 50:2.33,
    60:2.93, 70:3.53, 80:4.14, 90:4.74, 99:5.28}
id2s_list=[99,90,80,70,60,50,40,30]

kClust_template=Template("$kClust -i $infile -d $tmpdir/kClust -s $s -M 5000MB")
kClust_mkAln_template=Template("$kClust_mkAln -c '$clustalo --threads=$ncpu -i $$infile -o $$outfile' -d $tmpdir/kClust --no-pseudo-headers|grep -P '^Filename:'|cut -d' ' -f2")
qClust_mkAln_template=Template("$qClust_mkAln -c '$clustalo --threads=$ncpu -i $$infile -o $$outfile' -d $tmpdir/kClust --no-pseudo-headers; $clustalo --threads=$ncpu -i $tmpdir/kClust/ComInpList.fasta -o $tmpdir/kClust/ComOutList.fas --force -b")
reformat_template=Template("$reformat fas a3m $tmpdir/kClust/AllOutList.fas $tmpdir/list.a3m -batch -v 0")
hhblitsdb_template=Template("$hhblitsdb -o $outdb -ia3m $tmpdir/a3m")
cdhit_template=Template("$cdhit -i $infile -o $outfile -c $c -n $n -T $ncpu -M 8000")

# whether to use kClust_mkAln + clustalo or qClust_mkAln + clustalo -b for
# intra-cluster alignment. The latter is slightly faster but the results are
# slightly different.
mkAln=1 # 0 - kClust_mkAln; 1 - qClust_mkAln

def mkdir_if_not_exist(tmpdir):
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

def make_tmpdir(tmpdir):
    ''' create tmp folder '''
    if not tmpdir:
        import random
        tmpdir="/tmp/%s/kClust_%s"%(
            os.getenv("USER"),random.randint(0,10**10))
        while(os.path.isdir(tmpdir)):
            tmpdir="/tmp/%s/kClust_%s"%(
                os.getenv("USER"),random.randint(0,10**10))
    mkdir_if_not_exist(tmpdir)
    sys.stdout.write("created folder %s\n"%tmpdir)
    return tmpdir 

def remove_a3m_gap(infile,outfile,seqname_prefix=''):
    ''' read a3m/fasta format infile, remove gaps and output to outfile. 
    return the number of sequences '''
    fp=open(infile,'rU')
    lines=fp.read().splitlines()
    fp.close()
    txt=''
    Nseq=0
    for line in lines:
        if line.startswith('>'):
            Nseq+=1
            line='>'+seqname_prefix+line[1:]
        else:
            line=line.upper().replace('-','').replace('.','')
        txt+=line+'\n'
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return Nseq

def remove_redundant_cdhit(infile,outfile,cdhit_c,ncpu):
    ''' read fasta format infile, perform redundancy removal, and 
    output to outfile.  return the number of sequences '''
    n=5
    if cdhit_c<0.5:
        n=2
    elif cdhit_c<0.6:
        n=3
    elif cdhit_c<0.7:
        n=4

    cmd=cdhit_template.substitute(
        cdhit=bin_dict["cdhit"],
        infile=infile,
        outfile=outfile,
        c=cdhit_c,
        n=n,
        ncpu=ncpu,
    )
    os.system(cmd)

    if os.path.isfile(outfile):
        infile=outfile.strip()
    fp=open(outfile,'rU')
    Nseq=('\n'+fp.read()).count("\n>")
    return Nseq,infile

def kClust2db(infile,outdb,tmpdir='.',s=1.12,ncpu=1):
    ''' cluster sequences in FASTA file "infile", and generate hhblits
    style database at outdb'''
    sys.stdout.write("#### cluster input fasta ####\n")
    cmd=kClust_template.substitute(dict(
        kClust=bin_dict["kClust"],
        infile=infile,
        tmpdir=tmpdir,
        s=s,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)

    sys.stdout.write("#### alignment within each cluster ####\n")
    cmd=kClust_mkAln_template.substitute(dict(
        kClust_mkAln=bin_dict["kClust_mkAln"],
        clustalo=bin_dict["clustalo"],
        ncpu=ncpu,
        tmpdir=tmpdir,
    ))
    if mkAln:
        cmd=qClust_mkAln_template.substitute(dict(
            qClust_mkAln=bin_dict["qClust_mkAln"],
            clustalo=bin_dict["clustalo"],
            ncpu=ncpu,
            tmpdir=tmpdir,
        ))

    sys.stdout.write(cmd+'\n')
    stdout,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE).communicate()

    sys.stdout.write("#### reformat fas into a3m ####\n")
    a3mdir=os.path.join(tmpdir,"a3m")
    mkdir_if_not_exist(a3mdir)
    if mkAln:
        fp=open(os.path.join(tmpdir,"kClust/AllOutList.fas"),'r')
        lines=fp.read().splitlines()
        fp.close()
    else:
        stdout=stdout.decode()
        fp=open(os.path.join(tmpdir,"kClust/AllOutList.fas"),'w')
        fp.write(stdout)
        fp.close()
        lines=stdout.splitlines()
    fp=open(os.path.join(tmpdir,"list.a3m"),'w')
    for filename in lines:
        basename=os.path.basename(os.path.splitext(filename)[0])
        fp.write(os.path.join(tmpdir,"a3m",basename+".a3m")+'\n')
    fp.close()
    os.system(reformat_template.substitute(dict(
        reformat=bin_dict["reformat"], tmpdir=tmpdir)))

    sys.stdout.write("#### build hhblitsdb ####\n")
    mkdir_if_not_exist(os.path.dirname(outdb))
    cmd=hhblitsdb_template.substitute(dict(
        hhblitsdb=bin_dict["hhblitsdb"],
        ncpu=ncpu,
        outdb=outdb,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return

def decideKclustID(Nseq):
    seqID=id2s_list[-1]
    for i,ids in enumerate(id2s_list):
        if Nseq<(i+1)*500:
            seqID=ids
            break
    return seqID

def kClust2db_main(infile,outdb,seqID,ncpu,tmpdir,cdhit_c=1.0):
    mkdir_if_not_exist(tmpdir)
    fastafile=os.path.join(tmpdir,"db.fasta")
    Nseq=remove_a3m_gap(infile,fastafile)
    if cdhit_c<=1:
        cdhitfile=os.path.join(tmpdir,"cdhit.fasta")
        Nseq,fastafile=remove_redundant_cdhit(fastafile,cdhitfile,cdhit_c,ncpu)
    print(cdhit_c,fastafile,Nseq)

    if seqID==0: # detemine seqID automatically
        seqID=decideKclustID(Nseq)

    s=id2s_dict[seqID]
    print("Nseq=%d. kClust -s %d for seqID %d%%"%(Nseq,s,seqID))

    kClust2db(fastafile,outdb,tmpdir,s,ncpu)
    return

if __name__=="__main__":
    seqID=0
    ncpu=1
    tmpdir=''
    cdhit_c=1.0
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-id="):
            seqID=float(arg[len("-id="):])
            if seqID<=1:
                seqID=100*seqID
            seqID=int(seqID)
            if not seqID in id2s_dict:
                sys.stderr.write("ERROR! Illegal sequence identity cutoff %d\n"%seqID)
            exit()
        elif arg.startswith("-ncpu="):
            ncpu=int(arg[len("-ncpu="):])
        elif arg.startswith("-c="):
            cdhit_c=float(arg[len("-c="):])
            if cdhit_c>1:
                cdhit_c/=100.
            if cdhit_c<0.4:
                sys.stderr.write("ERROR! -c <40: %s\n"%arg)
                exit()
        elif arg.startswith("-tmpdir="):
            tmpdir=os.path.abspath(arg[len("-tmpdir="):])
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! No such option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)!=2:
        sys.stderr.write(docstring)
        exit()

    infile=os.path.abspath(argv[0])
    outdb=os.path.abspath(argv[1])
    tmpdir=make_tmpdir(tmpdir)

    kClust2db_main(infile,outdb,seqID,ncpu,tmpdir,cdhit_c)

    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
