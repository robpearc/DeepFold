#!/usr/bin/env python
docstring='''
build_MSA.py seq.fasta            \\
    -hhblitsdb=uniclust30_2017_10 \\
    -jackhmmerdb=uniref90.fasta   \\
    -hmmsearchdb=metaclust.fasta

    build PSICOV format MSA using hhblits, jackhmmer, or hmmsearch.

options:
    (at least one of -hhblitsdb, -jackhmmerdb or -hmmsearchdb must be set)

    -hhblitsdb=uniclust30_2017_10
        hhsuite databases, to be searched by hhblits.

    -jackhmmerdb=uniref90.fasta
        (decompressed) fasta databases, to be searched by jackhmmer, whose
        search result will be built into a custom hhsuite database to be 
        searched by hhblits, jump-starting from alignment generated by 
        searching hhblitsdb. fasta database must have an ssi index file
        created by esl-sfetch.

    -hmmsearchdb=metaclust.fasta:tara.fasta
        colon deliminated list of decompressed fasta database, to be
        searched hmmsearch, jump-starting from alignment generated by 
        searching either -hhblitsdb or -jackhmmerdb. fasta database must
        have an ssi index file created by esl-sfetch.

    -tmpdir=/tmp/$USER/MSA_`date +%N`
        temporary folder

    -outdir=.
        output folder. default is current folder

    -overwrite={0,1,2,4}
        whether overwrite existing search result.
        0 - do not overwrite any intermediate alignment
        1 - overwrite hhblitsdb search result (.hhbaln and .hhba3m)
        2 - overwrite jackhmmerdb search result (.jacaln and .jaca3m)
        4 - overwrite hmmsearchdb search result (.hmsaln)
        These options are addictive. For example, -overwrite=7 (=1+2+4) for
        overwriting any intermediate alignment (but might still filter
        final alignment if it is too large).

    -ncpu=1
        number of CPU threads. do not use multi-threading by default.

output:
    (filename prefix determined by input filename)
    seq.aln     - final alignment. the only non-optional output
    seq.hhbaln  - (if -hhblitsdb is set) hhblits MSA (PSICOV format)
    seq.hhba3m  - (if -hhblitsdb is set) hhblits MSA (a3m format)
    seq.jacaln - (if -jackhmmerdb is set and -hhblitsdb search does not have
                  enough sequences) jackhmmer + hhblits MSA (PSICOV format)
    seq.jaca3m - (if -jackhmmerdb is set and -hhblitsdb search does not have
                  enough sequences) jackhmmer + hhblits MSA (a3m format)
    seq.hmsaln  - (if -hmmsearchdb is set and neither -hhblitsdb nor
                  -jackhmmerdb search has enough sequences)
                  hhblits + jackhmmer (optional) + hmmsearch output
'''

import sys, os
import shutil
import subprocess
from string import Template,ascii_lowercase

from HHPaths import bin_dict
from kClust2db import kClust2db,id2s_dict

#### search parameters ####
# min query coverage by template
cov_cut=[
    50, # from metapsicov 2.0.3 and deepcontact. for final alignment output
    60, # from metapsicov 1.04. for Nf calculation
    75, # from gremlin. for qhmmbuild searching
]

# target nf, beyond which we do not attempt to build deeper alignment
target_nf=[
    128,  # 64 from gremlin
    129, # 128 from baker casp12
]

# max seqID within MSA
id_cut=[
    99, # from metapsicov and deepcontact
    90, # from gremlin
]

# which file is used by qhmmbuild to build query hmm
# 1  - (recommended) use PSICOV format .aln alignment
# 0  - use hhsuite format .a3m alignment, takes more resources with
#      slightly worse result
aln2hmm=1 

# whether to build custom hhsuite database from qhmmsearch hits
# 0 - directly use filtered and realigned qhmmsearch MSA as output
# 1 - build custom hhsuite database from qhmmsearch hits, and search
#     this custom database with hhblits, similar to jackblits
# 2 - build custom hhsuite database from qhmmsearch hits as well as
#     hhblits/jackblits hits, and search this custom database with hhblits
build_hmmsearch_db=1

# whether to further filter qhmmsearch+hhblits result using more stringent
# id_cut and cov_cut
# 0 - no further filter
# 1 - filter if nf>=target_nf
# 2 - always filter
filter_hmsblits=0

# The threshold above which kClust'ed hhsuite database is build instead of
# non-clustered hhsuite database. For some reason, non-clustered database
# has slightly better alignment quality, but takes much longer time to make.
kClust2db_threshold=1000

# The maximum number of hits to be parsed in a checkpoint alignment
checkali_threshold=30000

#### command templates ####
# Although both gremlin and PconSC2 recommend increasing -maxfilt, it seems
# higher -maxfilt increases alignment depth without significant benefit for
# improving contact prediction accuracy. Therefore, -maxfilt is still
# default value here.
# $infile      - input fasta or a3m
# $db          - hhsuite database
# $ncpu        - number of cpu cores
# $outprefix   - outputs are $outprefix.a3m $outprefix.log $outprefix.aln
# $id_cut      - percentage max seqID among aligned sequences
# $cov_cut     - percentage min cov of query by templates
hhblits_template=Template(bin_dict["hhblits"]+ \
" -i $infile -diff inf -d $db -cpu $ncpu -oa3m $outprefix.a3m "+ \
" -id $id_cut -cov $cov_cut -o $outprefix.log -n 3 -diff inf; "+ \
"grep -v '^>' $outprefix.a3m|sed 's/[a-z]//g' > $outprefix.aln")

# While both hhfilter and rmRedundantSeq can control the max seqID and min
# coverage within an alignment, they are slightly different: hhfilter
# consider seqID normalized by full query length, while rmRedundantSeq
# consider seqID normalized by aligned region only.
# $prefix      - input in $prefix.a3m (hhfilter_template) or
#                in $prefix.aln (alnfilter_template)
#                output are $prefix.$cov_cut.a3m and $prefix.$cov_cut.aln
# $id_cut      - percentage max seqID among aligned sequences
# $cov_cut     - percantage min cov of query by templates
hhfilter_template=Template(bin_dict["hhfilter"]                         + \
" -i $prefix.a3m -o $prefix.$cov_cut.a3m -id $id_cut -cov $cov_cut; "   + \
"grep -v '^>' $prefix.$cov_cut.a3m|sed 's/[a-z]//g' > $prefix.$cov_cut.aln")

alnfilter_template=Template(bin_dict["rmRedundantSeq"]+ \
" $id_cut $cov_cut $prefix.aln > $prefix.$cov_cut.aln")

# $calNf       - calNf executable
# $infile      - PSICOV format input
# $target_nf   - max Nf to consider
calNf_template=Template(bin_dict["calNf"]+" $infile 0.8 0 $target_nf")

# qjackhmmer is almost identical to jackhmmer, except that the output
# alignment only include positions covered by query. Here, the number of
# sequences in ".tbl" and ".fseqs" files is usually larger than in ".first"
# file, because tbl and fseqs contain all hits with -E <= 10, while
# ".first" only contain --incE <= 1e-3. Short queries might hit a long-
# multidomain template with significant -E but not significant --incE. This
# neccesitate trim_eslsfetch for removing spurious fseqs hits and template
# regions that are neither aligned nor close to an aligned region.
# $ncpu        - number of cpu cores
# $outprefix   - outputs are $outprefix.first $outprefix.tbl 
#                $outprefix.fseqs $outprefix.out
# $infile      - input query fasta
# $db          - jackhmmer database
qjackhmmer_template=Template(bin_dict["qjackhmmer"]+" --cpu $ncpu"   + \
" -N 3 -E 10 --incE 1e-3 -A $outprefix.first --tblout $outprefix.tbl"+ \
" -o $outprefix.out $infile $db; "+bin_dict["eslsfetch"]             + \
" -f $db $outprefix.tbl|sed 's/*//g' > $outprefix.fseqs")

# hhblitsdb.pl can take much longer time than even jackhmmer if there
# are many hits. In this case, the hits are first clustered by kClust
# before being used by hhblitsdb.pl.
# $ncpu        - number of cpu cores
# $db          - output hhblits database
# $a3mdir      - input a3m folder
hhblitsdb_template=Template(
bin_dict["hhblitsdb"]+" -cpu $ncpu -o $db -ia3m $a3mdir")

# qhmmbuild is almost identical to hmmbuild, except that the output
# alignment only contains the first sequence in query alignment, so
# as to save disk space. qhmmbuild can be replaced by hmmbuild.
# $infile      - input a3m (qhmmbuild_a3m_template) or 
#                input aln (qhmmbuild_aln_template)
# $outprefix   - outputs are $outprefix.afq $outprefix.hmm
qhmmbuild_aln_template=Template(
"sed = $infile |sed 'N;s/\\n/\\t/'|sed 's/^/>/g'|sed 's/\\t/\\n/g'|"+ \
bin_dict["qhmmbuild"]                                               + \
" -n aln --amino -O $outprefix.afq --informat afa $outprefix.hmm -")

qhmmbuild_a3m_template=Template(
bin_dict["reformat"]+" a3m a2m $infile -|grep -vP "                 + \
"'^(Reformat|Using |Removed |Sequence |inserting |WARNING: )'|"     + \
bin_dict["qhmmbuild"]                                               + \
" -n a3m --amino -O $outprefix.afq --informat afa $outprefix.hmm -")

# This only check number of match states that are in query. Match states
# corresponding to gap, i.e. [-], in query are not considered.
# $infile      - input fasta alignment
num_match_state_template=Template(
    bin_dict["fasta2aln"]+" $infile|head -1|grep -ohP '[A-Z]'|wc -l")

# qhmmsearch is almost identical to hmmsearch, except that the output
# alignment exclude insertions, which saves lots of disk space.
# $outprefix   - outputs are $outfile.match $outfile.tbl $outprefix.out
# $infile      - input hmm
# $db          - hmmsearch database
# $ncpu        - number of cpu cores
qhmmsearch_template=Template(bin_dict["qhmmsearch"]+" --cpu $ncpu"   + \
" -A $outprefix.match --incT 27 -T 27 --incdomT 27 -o $outprefix.out"+ \
" --tblout $outprefix.tbl $infile $db")

qhmmsearch_eslsfetch_template=Template(bin_dict["qhmmsearch"]+" --cpu"+ \
" $ncpu -E 10 --incE 1e-3 -A $outprefix.match --tblout $outprefix.tbl"+ \
" -o $outprefix.out $infile $db; "+bin_dict["eslsfetch"]+ \
" -f $db $outprefix.tbl|sed 's/*//g' > $outprefix.fseqs")


# The HMM built by qhmmbuild is usually shorter than the query. This is
# because some positions are considered match states while other positions
# are considered insertion states. This result in the qhmmsearch alignment
# being shorter than the query. realignMSA program re-align the qhmmsearch
# alignment to query so that the former has the same length as query.
# $input_match - input fasta alignment
# $input_afq   - input match states for query
# $outfile     - output alignment
# $cov_cut1    - coverage of match state positions
# $cov_cut2    - coverage of full query. $cov_cut2<=$cov_cut1
realignMSA_template=Template("cat $input_match|sed 's/[*JUZBOX]/-/g'|"+ \
bin_dict["fastaCov"]+" $cov_cut1|"                                    + \
bin_dict["realignMSA"]+" $input_afq - |"                              + \
bin_dict["fastaCov"]+" $cov_cut2|"                                    + \
bin_dict["fasta2aln"]+" -|tr '[:lower:]' '[:upper:]' > $outfile")

# $id_cut      - max seqID at aligned region
# $cov_cut     - min cov to query by new template
# $infile1     - old alignment from hhblits or jachmmer
# $infile2     - new alignment from hmmsearch
# $outprefix   - $outprefix.nonredundant is alignment non-redundant
#                to $infile1 and non-redundant within the alignment.
#                $outprefix.aln combines $infile1 and $outprefix.nr
rmRedundantSeq_template=Template(bin_dict["rmRedundantSeq"]      + \
" $id_cut $cov_cut $infile1 $infile2 > $outprefix.nonredundant;" + \
" cat $infile1 $outprefix.nonredundant > $outprefix.aln")

#### check databases ####

def check_db(db_dict):
    ''' check if databases are legal '''
    if not db_dict["hhblitsdb"]:
        sys.stderr.write("Warning! Not setting -hhblitsdb is unrecommended\n")
        if not db_dict["jackhmmerdb"]:
            if db_dict["hmmsearchdb"]:
                sys.stderr.write("WARNING! Using single query "+ \
                    "sequence for hmmsearch generates worse result\n")
            else:
                sys.stderr.write("Please at least set -hhblitsdb\n")
                sys.stderr.write("ERROR! No database to search\n")
                exit()

    if db_dict["hhblitsdb"]:
        a3m_db=db_dict["hhblitsdb"]+"_a3m_db"
        if not os.path.isfile(a3m_db):
            sys.stderr.write(
                "ERROR! Cannot locate %s for\n -hhblitsdb=%s\n"%(
                a3m_db,db_dict["hhblitsdb"]))
            exit()

    if db_dict["jackhmmerdb"]:
        for db in db_dict["jackhmmerdb"].split(':'):
            if not os.path.isfile(db):
                sys.stderr.write("ERROR! No such jackhmmerdb file %s\n"%db)
                exit()

            if db.endswith(".gz") or db.endswith(".bz2"):
                sys.stderr.write(
                    "ERROR! jackhmmer cannot use compressed file %s\n"%db)
                exit()

            ssi_db=db+".ssi"
            if not os.path.isfile(ssi_db):
                sys.stderr.write(
                    "ERROR! ssi index file missing for %s.\n"%db+ \
                    "You can create the ssi index file by:\n"   + \
                    "%s --index %s\n"%(bin_dict["eslsfetch"],db))
                exit()

    if db_dict["hmmsearchdb"]:
        for db in db_dict["hmmsearchdb"].split(':'):
            if not os.path.isfile(db):
                sys.stderr.write("ERROR! No such hmmsearchdb file %s\n"%db)
                exit()

            if not build_hmmsearch_db:
                continue

            if db.endswith(".gz") or db.endswith(".bz2"):
                sys.stderr.write(
                    "ERROR! esl-sfetch cannot use compressed file %s\n"%db)
                exit()

            ssi_db=db+".ssi"
            if not os.path.isfile(ssi_db):
                sys.stderr.write(
                    "ERROR! ssi index file missing for %s.\n"%db+ \
                    "You can create the ssi index file by:\n"   + \
                    "%s --index %s\n"%(bin_dict["eslsfetch"],db))
                exit()
    return
    
#### parse query ####

def read_one_sequence(query_fasta="seq.fasta"):
    ''' check if input is legal single sequence fasta and read the sequence '''
    fp=open(query_fasta,'rU')
    txt=fp.read()
    fp.close()
    if ('\n'+txt).count('\n>')!=1:
        sys.stderr.write("ERROR! Input is not single sequence fasta.")
        exit()
    sequence=''
    for line in txt.splitlines():
        if not line.startswith('>'):
            sequence+=line.strip()
    sequence=sequence.upper().replace(' ','').replace('\t','')
    illegal_residues=set(sequence)-set("ABCDEFGHIKLMNOPQRSTUVWXYZ")
    if illegal_residues:
        sys.stderr.write("ERROR! %s contains illegal residues %s\n"%(
            query_fasta,' '.join(illegal_residues)))
        exit()
    return sequence

#### make tmp folder ####

def mkdir_if_not_exist(tmpdir):
    ''' create folder if not exists '''
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    
def make_tmpdir(tmpdir):
    ''' create tmp folder '''
    if not tmpdir:
        import random
        tmpdir="/tmp/%s/MSA_%s"%(
            os.getenv("USER"),random.randint(0,10**10))
        while(os.path.isdir(tmpdir)):
            tmpdir="/tmp/%s/MSA_%s"%(
                os.getenv("USER"),random.randint(0,10**10))
    mkdir_if_not_exist(tmpdir)
    sys.stdout.write("created folder %s\n"%tmpdir)
    return tmpdir 

#### search sequence database ####

def getNf(prefix):
    ''' return Nf on -cov 60 filtered MSA. input file is prefix.a3m.
    output files are prefix.60.a3m and prefix.60.aln
    '''
    cov   =cov_cut[1]
    infile="%s.%d.aln"%(prefix,cov)

    #### filter input alignment ####
    if not os.path.isfile(infile):
        if os.path.isfile(prefix+".a3m"):
            cmd=hhfilter_template.substitute(
                prefix =prefix,
                id_cut =id_cut[0],
                cov_cut=cov,
            )
        else:
            cmd=alnfilter_template.substitute(
                prefix=prefix,
                id_cut=id_cut[0],
                cov_cut=cov,
            )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

    #### calculate Nf ####
    cmd=calNf_template.substitute(
        infile   =infile,
        target_nf=target_nf[-1],
    )
    stdout,stderr=subprocess.Popen(cmd,
        shell=True,stdout=subprocess.PIPE).communicate()
    return float(stdout)

def run_hhblits(query_fasta,db,ncpu,hhblits_prefix):
    ''' run hhblits with -cov 50 and return Nf for -cov 60'''
    cmd=hhblits_template.substitute(
        infile   =query_fasta,
        db       =db,
        ncpu     =ncpu,
        outprefix=hhblits_prefix,# outputs are $outprefix.a3m 
                                 # $outprefix.log $outprefix.aln
        id_cut   =id_cut[0],     # 99 in metapsicov
        cov_cut  =cov_cut[0],    # 50 in metapsicov 2.0.3
    )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return getNf(hhblits_prefix)

def fasta2a3msplit(txt,a3mdir='.'):
    ''' split multiple sequence FASTA txt into single sequence a3m files
    under "a3mdir"'''
    mkdir_if_not_exist(a3mdir)
    seqnum=0
    for block in ('\n'+txt).split('\n>'):
        lines=block.splitlines()
        if len(lines)<2:
            continue
        name=lines[0].split()[0].replace('|','_')
        sequence=''.join(lines[1:])
        fp=open(os.path.join(a3mdir,name+".a3m"),'w')
        fp.write('>'+lines[0]+'\n'+sequence+'\n')
        fp.close()
        seqnum+=1
    return seqnum

def trim_eslsfetch(fseqs_file,first_file,L=0,seqname_prefix='',
    max_seqnum=0):
    ''' 
    trim fseqs_file according to first_file so that for a template in
    fseqs_file, the N and C termini are trimmed such that template 
    sequence flanking the aligned region is <=L at each side.

    fseqs_file - fasta fetched by esl-efetch result
    first_file - alignment output by qjackhmmer
    L          - length of query. If L is 0, its value is inferred from the
                 first sequence in first_file. This first sequence will be
                 assumed to be the query. If L is not 0, all sequences in
                 first_file will be assumed to be templates
    seqname_prefix - append this to the begining of each sequence name. This
                 is mainly for avoiding name conflict, i.e. sequences from
                 different databases sharing the same name.
    max_seqnum - maximum number of hits to be parsed in first_file
                 default is 0, which means parsing all sequences
    '''
    trim_txt=''
    L=0
    trim_dict=dict() # key is sequence name, value (min,max) position to keep

    fp=open(first_file,'rU')
    seqnum=0
    for block in ('\n'+fp.read()).split("\n>"):
        lines=block.splitlines()
        if len(lines)<2:
            continue
        if not L:
            sequence=''.join(lines[1:])
            L=len(sequence)
            trim_dict[lines[0].split()[0]]=(1,L)
        else:
            subseq_name=lines[0].split()[0]
            name,pos_range=subseq_name.split('/')
            min_pos,max_pos=map(int,pos_range.split('-'))
            min_pos=max([0,min_pos-1-L]) # starting from 0
            max_pos+=L                   # starting from 1
            if name in trim_dict:
                min_pos=min([min_pos,trim_dict[name][0]])
                max_pos=max([max_pos,trim_dict[name][1]])
            trim_dict[name]=(min_pos,max_pos)
        seqnum+=1
        if max_seqnum>0 and seqnum>=max_seqnum:
            break
    fp.close()

    fp=open(fseqs_file,'rU')
    for block in ('\n'+fp.read()).split("\n>"):
        lines=block.splitlines()
        if len(lines)<2:
            continue
        sequence=''.join(lines[1:])
        name=lines[0].split()[0]
        if name in trim_dict:
            trim_txt+='>'+seqname_prefix
            # hhmake "helpfully" assume a3m with sequences name ending in
            # "_consensus" must have a master sequence. see "hhalignment.C"
            if "_consensus" in name:
                trim_txt+=name.replace("_consensus","_consen")+'\n'
            else:
                trim_txt+=name+'\n'
            trim_txt+=sequence[trim_dict[name][0]:trim_dict[name][1]]+'\n'
    fp.close()
    return trim_txt

def run_jackblits(query_fasta,db_list,ncpu,hhblits_prefix,jackblits_prefix):
    ''' run jackhmmer, hhblitsdb and hhblits '''
    fp=open(query_fasta,'rU')
    txt=fp.read()
    fp.close()

    #### run jackhmmer ####
    for d,db in enumerate(db_list):
        # outputs are $outprefix.first $outprefix.tbl $outprefix.fseqs
        outprefix=jackblits_prefix+'.'+str(d)
        cmd=qjackhmmer_template.substitute(
            ncpu     =ncpu,
            outprefix=outprefix,
            infile   =query_fasta,
            db       =db,
        )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

        # parse jackhmmer hits
        txt+=trim_eslsfetch(outprefix+".fseqs",outprefix+".first",
            seqname_prefix="jac%d_"%d,
            max_seqnum=checkali_threshold, # avoid excessive number of hits
        )
    
    #### build hhsuite database ####
    db=jackblits_prefix+"-mydb/mydb"
    if txt.count('\n>')>kClust2db_threshold:
        ### cluster at 30% seqID and build database ###
        fp=open(jackblits_prefix+".fseqs",'w')
        fp.write(txt)
        fp.close()
        kClust2db(jackblits_prefix+".fseqs",db,s=id2s_dict[30],
            ncpu=ncpu,tmpdir=os.path.dirname(jackblits_prefix))
    else:
        ### split jackhmmer hits into a3m ###
        a3mdir=jackblits_prefix+"-mya3m"
        seqnum=fasta2a3msplit(txt,a3mdir)
    
        ### build single sequence profile database ###
        mkdir_if_not_exist(jackblits_prefix+"-mydb")
        cmd=hhblitsdb_template.substitute(
            ncpu  =ncpu,
            db    =db,
            a3mdir=a3mdir,
        )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

    #### hhblits search  ####
    if os.path.isfile(hhblits_prefix+".a3m"):
        query_fasta=hhblits_prefix+".a3m"
    else:
        sys.stderr.write("WARNING! Using single sequence for jackblits\n")
    return run_hhblits(query_fasta,db,ncpu,jackblits_prefix)

def run_hmmsearch(query_fasta,sequence,hhblits_prefix,db_list,ncpu,
    hmmsearch_prefix): # build_hmmsearch_db==0 in search_metaclust
    L=len(sequence)

    #### check number of matched states in query #####
    cmd=num_match_state_template.substitute(
        infile=hmmsearch_prefix+".afq")
    sys.stdout.write(cmd+'\n')
    stdout,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE).communicate()
    Lmatch=int(stdout)
    target_Lmatch=(1.*cov_cut[-2]/cov_cut[-1])

    #### re-generate hmm if too few match states ####
    if (Lmatch < (target_Lmatch*L)):
        sys.stdout.write("WARNING! HMM length %.2f L < %.2f L\n"%(
            1.*Lmatch/L,target_Lmatch))
        cov=cov_cut[-1]
        if aln2hmm:
            cmd=';'.join([
                alnfilter_template.substitute(
                    prefix =hhblits_prefix,
                    cov_cut=cov,
                    id_cut =100,
                ),
                qhmmbuild_aln_template.substitute(
                    infile   =hhblits_prefix+".%d.aln"%cov,
                    outprefix=hmmsearch_prefix,
                ),
            ])
        else:
            cmd=';'.join([
                hhfilter_template.substitute(
                    prefix =hhblits_prefix,
                    cov_cut=cov,
                    id_cut =100,
                ),
                qhmmbuild_a3m_template.substitute(
                    infile   =hhblits_prefix+".%d.a3m"%cov,
                    outprefix=hmmsearch_prefix,
                ),
            ])
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

    #### search hmm against metaclust ####
    outfile=hmmsearch_prefix+".nonredundant"
    for d,db in enumerate(db_list):
        outprefix=hmmsearch_prefix+".%d"%d
        cmd=qhmmsearch_template.substitute(
            outprefix=outprefix,
            infile   =hmmsearch_prefix+".hmm",
            ncpu     =ncpu,
            db       =db,
        )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

    ### cat hmmsearch aln from multiple databases ###
    cmd="cat "+' '.join(["%s.%d.match"%(hmmsearch_prefix,d
        ) for d in range(len(db_list))])+" > %s.match"%hmmsearch_prefix
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
   
    ### realign MSA to query hmm ###
    cmd=realignMSA_template.substitute(
        input_match=hmmsearch_prefix+".match",
        input_afq  =hmmsearch_prefix+".afq",
        outfile    =hmmsearch_prefix+".redundant",
        cov_cut1   =cov_cut[-1],
        cov_cut2   =cov_cut[-2],
    )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)

    ### remove redundancy in resulting msa ###
    cmd=rmRedundantSeq_template.substitute(
        id_cut   =id_cut[-1], # typically the same or more stringent
                              # than that in hhblits
        cov_cut  =cov_cut[-2],# less stringenet than cov_cut1
        infile1  =hhblits_prefix+".aln",
        infile2  =hmmsearch_prefix+".redundant",
        outprefix=hmmsearch_prefix,
    )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return getNf(hmmsearch_prefix)

def run_hmsblits(query_fasta,sequence,hhblits_prefix,
    db_list,ncpu,hmmsearch_prefix): # build_hmmsearch_db in [1,2]
    #### read in query profile ####
    if build_hmmsearch_db==2 and os.path.isfile(hhblits_prefix+".a3m"):
        txt=''
        fp=open(hhblits_prefix+".a3m",'rU')
        for line in fp.read().splitlines():
            if line.startswith('>'):
                txt+=line+'\n'
            else:
                txt+=line.replace('-','').replace('.','').upper()+'\n'
        fp.close()
    else:
        fp=open(query_fasta,'rU')
        txt=fp.read()
        fp.close()

    #### search hmm against metaclust ####
    outfile=hmmsearch_prefix+".nonredundant"
    for d,db in enumerate(db_list):
        outprefix=hmmsearch_prefix+".%d"%d
        cmd=qhmmsearch_eslsfetch_template.substitute(
            outprefix=outprefix,
            infile   =hmmsearch_prefix+".hmm",
            ncpu     =ncpu,
            db       =db,
        )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)
        txt+=trim_eslsfetch(outprefix+".fseqs",outprefix+".match",
            L=len(sequence), seqname_prefix="hms%d_"%d,
            max_seqnum=checkali_threshold, # avoid excessive number of hits
        )

    #### build custom hhsuite database ####
    db=hmmsearch_prefix+"-mydb/mydb"
    if txt.count('\n>')>kClust2db_threshold:
        ### cluster at 30% seqID and build database ###
        fp=open(hmmsearch_prefix+".fseqs",'w')
        fp.write(txt)
        fp.close()
        kClust2db(hmmsearch_prefix+".fseqs",db,s=id2s_dict[30],
            ncpu=ncpu,tmpdir=os.path.dirname(hmmsearch_prefix))
    else:
        ### split hmmsearch hits into a3m ###
        a3mdir=hmmsearch_prefix+"-mya3m"
        seqnum=fasta2a3msplit(txt,a3mdir)
    
        ### build single sequence profile database ###
        mkdir_if_not_exist(hmmsearch_prefix+"-mydb")
        cmd=hhblitsdb_template.substitute(
            ncpu  =ncpu,
            db    =db,
            a3mdir=a3mdir,
        )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

    #### hhblits search  ####
    if not os.path.isfile(hhblits_prefix+".a3m") or \
       not os.path.isfile(hhblits_prefix+".aln"):
        sys.stderr.write("WARNING! Using single sequence for"+ \
            " hmmsearch+hhblits\n")
        return run_hhblits(query_fasta,db,ncpu,hmmsearch_prefix)
    
    hmsblits_nf=run_hhblits(hhblits_prefix+".a3m",db,ncpu,hmmsearch_prefix)

    ### remove redundancy in resulting msa ###
    if (filter_hmsblits==1 and hmsblits_nf>=target_nf[0]) or \
        filter_hmsblits==2:
        ### filter $hmmsearch_prefix.aln ###
        cmd=rmRedundantSeq_template.substitute(
            id_cut   =id_cut[-1], # typically the same or more stringent
                                  # than that in hhblits
            cov_cut  =cov_cut[0],
            infile1  =hhblits_prefix+".aln",
            infile2  =hmmsearch_prefix+".aln",
            outprefix=hmmsearch_prefix,
        )
        sys.stdout.write(cmd+'\n')
        os.system(cmd)
        
        ### update $hmmsearch_prefix.a3m accordingly ###
        txt=''
        fp=open(hmmsearch_prefix+".aln",'rU')
        aln_list=fp.read().splitlines()
        fp.close()

        fp_hhba3m=open(hhblits_prefix+".a3m",'rU')
        fp_hmsa3m=open(hmmsearch_prefix+".a3m",'rU')
        for block in ('\n'+fp_hhba3m.read()+fp_hmsa3m.read()).split('\n>'):
            lines=block.splitlines()
            if len(lines)!=2:
                continue
            header,sequence=lines
            match_seq=sequence.translate(None,ascii_lowercase)
            if match_seq in aln_list:
                txt+=">%s\n%s\n"%(header,sequence)
                for l in range(len(aln_list)):
                    if match_seq==aln_list[l]:
                        del aln_list[l]
                        break
        fp_hhba3m.close()
        fp_hmsa3m.close()

        fp=open(hmmsearch_prefix+".a3m",'w')
        fp.write(txt)
        fp.close()
        hmsblits_nf=getNf(hmmsearch_prefix)
    return hmsblits_nf

def search_metaclust(query_fasta,sequence,hhblits_prefix,
    db_list,ncpu,hmmsearch_prefix):
    ''' jump start hmmsearch using hhblits or jack_hhblits a3m'''

    #### convert input to hmm ####
    if aln2hmm:
        input_aln=hhblits_prefix+".aln"
        if not os.path.isfile(input_aln):
            sys.stderr.write("WARNING! No hhblits alignment.\n"+ \
                "Using query sequence as input for hmmsearch\n")
            fp=open(input_aln,'w')
            fp.write(sequence+'\n')
            fp.close()

        cmd=qhmmbuild_aln_template.substitute(
            infile=input_aln,
            outprefix=hmmsearch_prefix,
        )
    else:
        input_a3m=hhblits_prefix+".a3m"
        if not os.path.isfile(input_a3m):
            sys.stderr.write("WARNING! No hhblits alignment.\n"+ \
                "Using query sequence as input for hmmsearch\n")
            fp=open(input_a3m,'w')
            fp.write(">seq\n"+sequence+'\n')
            fp.close()

        cmd=qhmmbuild_a3m_template.substitute(
            infile   =input_a3m,
            outprefix=hmmsearch_prefix,
        )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)

    if not build_hmmsearch_db:
        return run_hmmsearch(query_fasta,sequence,hhblits_prefix,
            db_list,ncpu,hmmsearch_prefix)
    else:
        return run_hmsblits(query_fasta,sequence,hhblits_prefix,
            db_list,ncpu,hmmsearch_prefix)

def build_MSA(prefix, sequence, tmpdir, db_dict, ncpu=1,
    overwrite_dict=dict(hmmsearch=False,jackhmmer=False,hhblits=False)):
    ''' sequentially attempt to build MSA by hhblits, jackhmmer+hhblits,
    and hmmsearch. '''
    nf=hhb_nf=jack_nf=hms_nf=0
    #### preparing query ####
    query_fasta=os.path.join(tmpdir,"seq.fasta")
    fp=open(query_fasta,'w')
    fp.write(">seq\n%s\n"%sequence)
    fp.close()

    #### run hhblits ####
    hhblits_prefix=os.path.join(tmpdir,"hhblits")
    if db_dict["hhblitsdb"]:
        if overwrite_dict["hhblits"] or not os.path.isfile(
            prefix+".hhbaln") or not os.path.isfile(prefix+".hhba3m"):
            # generates hhblits_prefix.a3m hhblits_prefix.aln
            # hhblits_prefix.60.a3m hhblits_prefix.60.aln
            hhb_nf=run_hhblits(query_fasta,db_dict["hhblitsdb"],
                ncpu,hhblits_prefix)
            shutil.copyfile(hhblits_prefix+".aln",prefix+".hhbaln")
            shutil.copyfile(hhblits_prefix+".a3m",prefix+".hhba3m")
        else:
            shutil.copyfile(prefix+".hhbaln",hhblits_prefix+".aln")
            shutil.copyfile(prefix+".hhba3m",hhblits_prefix+".a3m")
            hhb_nf=getNf(hhblits_prefix)
            sys.stdout.write("%s and %s exists, skip hhblitsdb\n"%(
                prefix+".hhbaln",prefix+".hhba3m"))
        
        nf=hhb_nf
        if hhb_nf>=target_nf[0]:
            shutil.copyfile(hhblits_prefix+".aln",prefix+".aln")
            sys.stdout.write("Final MSA by hhblits with Nf >=%.1f\n"%nf)
            return nf

    #### run jack_hhblits ####
    jackblits_prefix=os.path.join(tmpdir,"jackblits")
    if db_dict["jackhmmerdb"]:
        if overwrite_dict["jackhmmer"] or not os.path.isfile(
            prefix+".jacaln") or not os.path.isfile(prefix+".jaca3m"):
            # generates jackblits_prefix.a3m jackblits_prefix.aln
            # jackblits_prefix.60.a3m jackblits_prefix.60.aln
            jack_nf=run_jackblits(query_fasta,
                db_dict["jackhmmerdb"].split(':'),
                ncpu,hhblits_prefix,jackblits_prefix)
            shutil.copyfile(jackblits_prefix+".aln",prefix+".jacaln")
            shutil.copyfile(jackblits_prefix+".a3m",prefix+".jaca3m")
        else:
            shutil.copyfile(prefix+".jacaln",jackblits_prefix+".aln")
            shutil.copyfile(prefix+".jaca3m",jackblits_prefix+".a3m")
            jack_nf=getNf(jackblits_prefix)
            sys.stdout.write("%s and %s exists, skip jackhmmerdb\n"%(
                prefix+".jacaln",prefix+".jaca3m"))
        
        nf=max([jack_nf,hhb_nf])
        if jack_nf>=target_nf[0]:
            shutil.copyfile(jackblits_prefix+".aln",prefix+".aln")
            sys.stdout.write("Final MSA by jackhmmer with Nf >=%.1f\n"%nf)
            return nf

    #### run hmmsearch ####
    hmmsearch_prefix=os.path.join(tmpdir,"hmmsearch")
    if db_dict["hmmsearchdb"]:
        if overwrite_dict["hmmsearch"] or \
            not os.path.isfile(prefix+".hmsaln") or \
            (build_hmmsearch_db and not os.path.isfile(prefix+".hmsa3m")):
            # generates hmmsearch_prefix.afq hmmsearch_prefix.hmm
            # hmmsearch_prefix.redundant hmmsearch_prefix.nonredundant
            # hmmsearch_prefix.aln
            hms_nf=search_metaclust(query_fasta,sequence,
                hhblits_prefix if jack_nf<hhb_nf else jackblits_prefix,
                db_dict["hmmsearchdb"].split(':'),ncpu,hmmsearch_prefix)
            shutil.copyfile(hmmsearch_prefix+".aln",prefix+".hmsaln")
            if os.path.isfile(hmmsearch_prefix+".a3m"):
                shutil.copyfile(hmmsearch_prefix+".a3m",prefix+".hmsa3m")
        else:
            shutil.copyfile(prefix+".hmsaln",hmmsearch_prefix+".aln")
            if os.path.isfile(prefix+".hmsa3m"):
                shutil.copyfile(prefix+".hmsa3m",hmmsearch_prefix+".a3m")
            sys.stdout.write("%s exists, skip hmmsearchdb\n"%(
                prefix+".hmsaln"))
            hms_nf=getNf(hmmsearch_prefix)

        if hms_nf>nf: # hmmsearch replaces jackblits and hhblits result
            nf=hms_nf
            shutil.copyfile(hmmsearch_prefix+".aln",prefix+".aln")
            sys.stdout.write("Final MSA by hmmsearch with Nf >=%.1f\n"%nf)
            return nf

    if hhb_nf>jack_nf:
        shutil.copyfile(hhblits_prefix+".aln",prefix+".aln")
        sys.stdout.write("hhblits MSA has %.1f Nf. Output anyway.\n"%hhb_nf)
    else:
        shutil.copyfile(jackblits_prefix+".aln",prefix+".aln")
        sys.stdout.write("jackhmmer MSA has %.1f Nf. Output anyway.\n"%jack_nf)
    return nf

def parse_overwrite_option(overwrite=0):
    ''' whether overwrite existing search result.
    0 - do not overwrite any alignment
    1 - overwrite hhblitsdb search result (.hhbaln and .hhba3m)
    2 - overwrite jackhmmerdb search result (.jacaln and .jaca3m)
    4 - overwrite hmmsearchdb search result (.hmsaln)
    These options are addictive, e.g., -overwrite=7 (=1+2+4) for
    overwriting any alignment. '''
    overwrite_dict=dict()

    overwrite_dict["hmmsearch"]=(overwrite>=4)
    overwrite %= 4

    overwrite_dict["jackhmmer"]=(overwrite>=2)
    overwrite %= 2

    overwrite_dict["hhblits"]=(overwrite>=1)
    overwrite %= 1
    return overwrite_dict

def refilter_aln(prefix,tmpdir):
    ''' filter final MSA by -id 99 -cov 60 '''
    #### filter MSA ####
    final_prefix=os.path.join(tmpdir,"final")
    shutil.copyfile(prefix+".aln",final_prefix+".aln")
    cov=cov_cut[1]
    cmd=alnfilter_template.substitute(
        id_cut=id_cut[0],
        cov_cut=cov,
        prefix=final_prefix,
    )
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    shutil.copyfile(final_prefix+".%d.aln"%cov,prefix+".aln")

    #### calculate Nf ####
    cmd=calNf_template.substitute(
        infile   =final_prefix+".%d.aln"%cov,
        target_nf=target_nf[-1],
    )
    stdout,stderr=subprocess.Popen(cmd,
        shell=True,stdout=subprocess.PIPE).communicate()
    nf=float(stdout)
    sys.stdout.write("Re-filter final MSA to Nf >= %.1f\n"%nf)
    return nf

if __name__=="__main__":
    #### command line argument parsing ####
    db_dict=dict(
        hhblitsdb='',
        jackhmmerdb='',
        hmmsearchdb='',
    )
    tmpdir=''
    outdir='.'
    overwrite=0
    ncpu=1

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-hhblitsdb="):
            db_dict["hhblitsdb"]=os.path.abspath(arg[len("-hhblitsdb="):])
        elif arg.startswith("-jackhmmerdb="):
            db_dict["jackhmmerdb"]=os.path.abspath(arg[len("-jackhmmerdb="):])
        elif arg.startswith("-hmmsearchdb="):
            db_dict["hmmsearchdb"]=os.path.abspath(arg[len("-hmmsearchdb="):])
        elif arg.startswith("-tmpdir="):
            tmpdir=os.path.abspath(arg[len("-tmpdir="):])
        elif arg.startswith("-outdir="):
            outdir=os.path.abspath(arg[len("-outdir="):])
        elif arg.startswith("-ncpu="):
            ncpu=int(arg[len("-ncpu="):])
        elif arg.startswith("-overwrite="):
            overwrite=int(arg[len("-overwrite="):])
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! No such option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)!=1 or (not ''.join(db_dict.values())):
        sys.stderr.write(docstring)
        exit()
    
    #### check input format ####
    query_fasta=os.path.abspath(argv[0])
    if not os.path.isfile(query_fasta):
        sys.stderr.write("ERROR! No such query fasta %s\n"%query_fasta)
        exit()
    check_db(db_dict)
    sequence=read_one_sequence(query_fasta)
    tmpdir=make_tmpdir(tmpdir)
    prefix=os.path.splitext(query_fasta)[0]
    if outdir and outdir!='.':
        mkdir_if_not_exist(outdir)
        prefix=os.path.join(outdir,os.path.basename(prefix))

    #### start building MSA ####
    nf=build_MSA(prefix,sequence,tmpdir,db_dict,
        ncpu=ncpu, overwrite_dict=parse_overwrite_option(overwrite))

    #### filter final MSA if too large ####
    # this will not improve contact accuracy. it is solely for making the
    # MSA not too large so that it is manageable for contact prediction
    if nf>=target_nf[-1]:
        nf=refilter_aln(prefix,tmpdir)

    #### clean up ####
    #if os.path.isdir(tmpdir):
        #shutil.rmtree(tmpdir)
