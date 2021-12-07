import numpy as np 
from numpy import float16,float32
import os,sys,gzip
from io import BytesIO
import getpass
import subprocess
#######config#############################
ccmpred=os.path.join(os.path.dirname(os.path.abspath(__file__)),'bin/ccmpred22')
hhmake=os.path.join(os.path.dirname(os.path.abspath(__file__)),'bin/hhsuite2/bin/hhmake')
if hhmake=='':
    print('Please set hhmakePATH at features.py')
    assert False
ALPHA=22
username=getpass.getuser()
## ccmpred -r dummy.raw T0862-D1.aln T0862-D1.ccmpred |gzip > raw.gz
#######config#############################
def plotme(matrix,saveimg):
    import matplotlib
    matplotlib.use('agg')
    from matplotlib import pyplot as plt
    plt.imshow(matrix)
    plt.savefig(saveimg)
    


def apc(pre):
    pre=np.abs(pre)
    n=pre.shape[0]
    pcmat=np.zeros([n,n])
    pcmysum=np.zeros(n)
    pcmean=0
    for i in range(n):
        for j in range(i+1,n):
            pcmysum[i]+=pre[i,j]
            pcmysum[j]+=pre[i,j]
            pcmean+=pre[i,j]
    pcmean=pcmean*1.0/(n*(n-1)*0.5)
    for i in range(n):
        for j in range(i+1,n):
            pcmat[i,j]=pre[i,j]-pcmysum[i]*pcmysum[j]/((n-1)*(n-1)*pcmean)
            pcmat[j,i]=pcmat[i,j]
    return pcmat
def blockshaped(arr,seq,dim=ALPHA):
    p=arr.shape[0]
    re=np.zeros([dim*dim+4,p,p])
    for i in range(p):
        for j in range(p):
            re[:dim*dim,i,j]=arr[i,j].flatten()
            re[dim*dim,i,j]=np.linalg.norm(re[:-1,i,j])
            re[dim*dim+1,i,j]=np.linalg.norm(arr[i,j,:-1,:-1].flatten())
            re[dim*dim+2,i,j]=np.linalg.norm(arr[i,j,:-2,:-2].flatten())
            re[dim*dim+3,i,j]=arr[i,j,aadic[seq[i]],aadic[seq[j]]]
    return re
def read_potts(lines,l1,states=ALPHA):
    lines=lines.split(b"\n")
    num_lines=len(lines)
    startindex=0
    for i in range(num_lines):
        if b'Final fx' in lines[i]:
            startindex=i
            break
    lines=lines[startindex+2:]
    #print(lines[0])
    precision=np.zeros([l1,(l1),states,states])
    #first:   
    count=0
    for i in range(l1):
        for k in range(i+1,l1):
            count+=1
            for j in range(ALPHA):
                vec=np.genfromtxt(BytesIO(lines[l1+(states+1)*(count-1)+1+j]))
                precision[i,k,j]=vec
                precision[k,i,:,j]=np.transpose(vec)
    for i in range(l1):
        vec=np.genfromtxt(BytesIO(lines[i]))
        #vec=np.append(vec,0)
        precision[i,i]=np.diag(vec)
    return precision
def getpotts2(alnfile,outfile,istraining):
    seq=open(alnfile).readline().strip()
    length=len(seq)
    if not istraining:
        cmd=ccmpred+' -r tmp '+alnfile+' '+outfile+'.ccm'+' |gzip >'+outfile+'.gz'
    else:
        cmd=ccmpred+' -r tmp  '+alnfile+' '+outfile+'.ccm'+' |gzip >'+outfile+'.gz'

    os.system(cmd)
    with gzip.open(outfile+'.gz', "rb") as f:
        plm=read_potts(f,length)
    os.remove(outfile+'.gz')
    plm=blockshaped(plm,seq)
    #np.save(outfile+'.npy',plm.astype(float16))
    return plm
def getpotts(alnfile,outfile,istraining):
    seq=open(alnfile).readline().strip()
    length=len(seq)
    if not istraining:
        cmd=ccmpred+' -r tmp '+alnfile+' '+outfile+'.ccm'+' '
    else:
        cmd=ccmpred+' -r tmp  '+alnfile+' '+outfile+'.ccm'

    #os.system(cmd)
    MyOut=subprocess.Popen(cmd, shell=True , stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdout,stderr = MyOut.communicate()
    plm=read_potts(stdout,length)
    plm=blockshaped(plm,seq)
    return plm

aadic = {
    'A': 0,
    'B': 20,
    'C': 4,
    'D': 3,
    'E': 6,

    'F': 13,
    'G': 7,
    'H': 8,
    'I': 9,
    'J': 20,

    'K': 11,
    'L': 10,
    'M': 12,
    'N': 2,
    'O': 20,

    'P': 14,
    'Q': 5,
    'R': 1,
    'S': 15,
    'T': 16,
    'U': 20,

    'V': 19,
    'W': 17,
    'X': 20,
    'Y': 18,
    'Z': 20,
    '-': 21,
    '*': 21,
}
def read_msa(file_path,aadic=aadic):
    """
    :param file_path:
    :param aadic:
    :return: msa of int [B,L]
    """
    lines=open(file_path).readlines()
    lines=[line.strip() for line in lines]
    n=len(lines)
    d=len(lines[0]) 
    msa=np.zeros([n,d],dtype=int)
    for i in  range(n):
        aline=lines[i]
        for j in range(d):
            msa[i,j]=aadic[aline[j]]
    return msa
def cal_large_mi(msa,weight,ALPHA=ALPHA):
    #output:441*l*l
    
    pseudoc=1
    M=msa.shape[0]
    N=msa.shape[1]
    pab=np.zeros((ALPHA,ALPHA))
    pa=np.zeros((N,ALPHA))
    cov=np.zeros([N,N,ALPHA,ALPHA])
    for i in range(N):
        for aa in range(ALPHA):
            pa[i,aa] = pseudoc
        neff=0.0
        for k in range(M):
            pa[i,msa[k,i]]+=weight[k]
            neff+=weight[k]
        for aa in range(ALPHA):
            pa[i,aa] /=pseudoc * ALPHA * 1.0 + neff
    #print(pab)
    for i in range(N):
        for j in range(i,N):
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if i ==j :
                        if a==b :
                            pab[a,b]=pa[i,a]
                        else:
                            pab[a,b]=0.0
                    else:
                        pab[a,b] = pseudoc *1.0 /ALPHA
            if(i!=j):
                neff2=0;
                for k in range(M):
                    a=msa[k,i]
                    b=msa[k,j]
                    tmp=weight[k]
                    pab[a,b]+=tmp
                    neff2+=tmp
                for a in range(ALPHA):
                    for b in range(ALPHA):
                        pab[a,b] /= pseudoc*ALPHA*1.0 +neff2
            for a in range(ALPHA):
                for b in range(ALPHA):
                    if(i!=j or a==b):
                        if (pab[a][b] > 0.0):
                            cov[i,j,a,b]=pab[a][b] * np.log(pab[a][b]/(pa[i][a] * pa[j][b]))      
                            cov[j,i,b,a]=cov[i,j,a,b]
    return cov 
def getcov(testaln,savefile,istraining):
    msa=read_msa(testaln)
    seq=open(testaln).readline().strip()
    num_aligns=len(open(testaln).readlines())
    if not istraining:
        if not os.path.isfile(savefile+'.weight'):
            exe=os.path.join(os.path.dirname(os.path.abspath(__file__)),'bin/colors_weight')
            cmd=exe+' '+testaln+' '+str(0.8)+' >'+savefile+'.weight'
            os.system(cmd)
        if not os.stat(savefile+'.weight').st_size==0:
            weight=np.genfromtxt(savefile+'.weight').flatten()
        else:
            print('can not get weight, using 1 as default')
            weight=np.ones(num_aligns)
    else:
        weight=np.ones(num_aligns)

    cov=cal_large_mi(msa,weight)
    cov=blockshaped(cov,seq)

    return np.concatenate([cov,np.ones([1,len(seq),len(seq)])*(weight.sum())],0)
def getgapcov(testaln,istraining):
    msa=read_msa(testaln)
    msa[msa<21]=0
    msa[msa==21]=1
    seq=open(testaln).readline().strip()
    num_aligns=len(open(testaln).readlines())
    weight=np.ones(num_aligns)
    gapcov=cal_large_mi(msa,weight,ALPHA=2)
    p=gapcov.shape[0]
    re=np.zeros([4,p,p])
    for i in range(p):
        for j in range(p):
            re[:,i,j]=gapcov[i,j].flatten()
    return re


def get2dbias(plm,cov,seq,dim=ALPHA,aadic=aadic):
    # this includes bias and one hot
    l=plm.shape[-1]
    plmdiag=[   np.diag(plm[:dim*dim,i,i].reshape([dim,dim]))     for i in range(l)]
    plmseq=[ plmdiag[i][aadic[seq[i]]]   for   i in range(l)           ]
    plmdiag=np.array(plmdiag)
    plmseq=np.array(plmseq).reshape([l,1])

    midiag=[   np.diag(cov[:dim*dim,i,i].reshape([dim,dim]))     for i in range(l)]
    miseq=[ midiag[i][aadic[seq[i]]]   for   i in range(l)           ]
    midiag=np.array(midiag)
    miseq=np.array(miseq).reshape([l,1])
    
    seq=np.array([aadic[aseq] for aseq in list(seq)])
    seq=np.eye(dim)[seq]
    fea1d=np.concatenate([plmdiag,plmseq,midiag,miseq,seq],1)
    return fea1d    
def read_hmm(testa3m,hhm_file):
    
    cmd=hhmake+' -i '+testa3m+' -o '+hhm_file
    os.system(cmd)

    f = open(hhm_file)
    line=f.readline()
    while line[0]!='#':
        line=f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    seq = []
    extras = np.zeros([0,10])
    prob = np.zeros([0,20])
    line = f.readline()
    while line[0:2]!='//':
        lineinfo = line.split()
        seq.append(lineinfo[0])  
        probs_ = [2**(-float(lineinfo[i])/1000) if lineinfo[i]!='*' else 0. for i in range(2,22)]
        prob = np.concatenate((prob,np.matrix(probs_)),axis=0)
        
        line = f.readline()
        lineinfo = line.split()
        extras_ = [2**(-float(lineinfo[i])/1000) if lineinfo[i]!='*' else 0. for i in range(0,10)]
        extras = np.concatenate((extras,np.matrix(extras_)),axis=0)
        
        line = f.readline()
        assert len(line.strip())==0
        
        line = f.readline()
    #return (''.join(seq),prob,extras)
    return np.concatenate((prob,extras),axis=1)

def collect_allfeatures_old(testaln,outprefix,istraining=False):
    # if os.path.isfile(outprefix+'.1d.npy') and os.path.isfile(outprefix+'.2d.npy'):
    #     fea1d=np.load(outprefix+'.1d.npy').astype(float32)
    #     fea2d=np.load(outprefix+'.2d.npy').astype(float32)
    #     return fea1d,fea2d
    lines=open(testaln).readlines()
    lines=[aline.strip() for aline in lines]
    testa3m=outprefix+'.a3m'
    wa3m=open(testa3m,'w')
    for aline in lines:
        wa3m.write('>temp\n')
        wa3m.write(aline+'\n')
    wa3m.close()
    seq=open(testaln).readline().strip()
    plm=getpotts(testaln,outprefix,istraining)
    cov=getcov(testaln,outprefix,istraining)
    
    bias=get2dbias(plm,cov,seq)
    hmm=read_hmm(testa3m,outprefix+'.hmm')

    fea1d=np.concatenate([bias,hmm],1)
    fea2d=np.concatenate([plm,cov],0)
    
    return fea1d.transpose(),fea2d
    
def collect_allfeatures(testa3m,outprefix,istraining=False):
    # if os.path.isfile(outprefix+'.1d.npy') and os.path.isfile(outprefix+'.2d.npy'):
    #     fea1d=np.load(outprefix+'.1d.npy').astype(float32)
    #     fea2d=np.load(outprefix+'.2d.npy').astype(float32)
    #     return fea1d,fea2d
    testaln=outprefix+'.aln'
    cmd="egrep -v \"^>\" "+testa3m+" | sed 's/[a-z]//g' >"+testaln
    os.system(cmd)

    seq=open(testaln).readline().strip()
    plm=getpotts(testaln,outprefix,istraining)
    cov=getcov(testaln,outprefix,istraining)
    
    bias=get2dbias(plm,cov,seq)
    hmm=read_hmm(testa3m,outprefix+'.hmm')

    fea1d=np.concatenate([bias,hmm],1)
    fea2d=np.concatenate([plm,cov],0)
    #np.save(outprefix+'.1d.npy',fea1d.astype(float16))
    #np.save(outprefix+'.2d.npy',fea2d.astype(float16))
    return fea1d.transpose(),fea2d
