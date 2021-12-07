import os,sys
#from numba import jit

import numpy as np 
from io import BytesIO
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
    '-': 20,
    '*': 20,
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
def generate2d(L):
    m2=np.zeros([2,L,L])
    for i in range(L):
        for j in range(i,L):
            m2[0,i,j]=i
            m2[1,i,j]=j
            m2[0,j,i]=j
            m2[1,j,i]=i
    return m2*0.0
def generate1d(L):
    m1=np.arange(L)
    m1=np.array([m1])
    return m1*0.0
def blockshaped(arr,dim=21):
    p=arr.shape[0]
    re=np.zeros([dim*dim,p,p])
    for i in range(p):
        for j in range(p):
            re[:,i,j]=arr[i,j].flatten()
    return re

def read_ccmprior(ccm1,l1):
    lines=open(ccm1).readlines()
    precision=np.zeros([l1,(l1),21,21])
    #first:   
    count=0
    for i in range(l1):
        for k in range(i+1,l1):
            count+=1
            for j in range(21):
                vec=np.genfromtxt(BytesIO(lines[l1+22*(count-1)+1+j].encode()))
                precision[i,k,j]=vec
                precision[k,i,:,j]=np.transpose(vec)
    for i in range(l1):
        vec=np.genfromtxt(BytesIO(lines[i].encode()))
        vec=np.append(vec,0)
        precision[i,i]=np.diag(vec)
    return precision
def computeccm_bin(msafile,savefile):
    #exefile='bin/ccmpred'
    seq=open(msafile).readlines()[0].strip()
    length=len(seq)
    exefile=os.path.join(os.path.dirname(os.path.abspath(__file__)),'bin/ccmpred')
    #exefile="/oasis/projects/nsf/mia174/liyangum/restriplet2/bin/ccmpred"
    cmd=exefile+' -r '+savefile+'.str'+' '+msafile+' '+savefile+'.del'
    #p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT,close_fds=True)
    if not os.path.isfile(savefile+'.str'):
        os.system(cmd) 
    precision=read_ccmprior(savefile+'.str',length)
    typed_pre=np.zeros([1,length,length])
    for i in range(length):
        for j in range(i,length):
            typed_pre[0,i,j]=precision[i,j,aadic[seq[i]],aadic[seq[j]]]
            typed_pre[0,j,i]=precision[j,i,aadic[seq[j]],aadic[seq[i]]]
    ccm=np.genfromtxt(savefile+'.del')
    ccm=np.expand_dims(ccm,0)
    precision=np.concatenate([blockshaped(precision),typed_pre,ccm])
    os.remove(savefile+'.str')
    os.remove(savefile+'.del')
    return precision,seq 
def get1d(ccmfea,mi,seq):
    l=len(seq)
    plmdiag=[   np.diag(ccmfea[2:441+2,i,i].reshape([21,21]))     for i in range(l)]
    plmseq=[ plmdiag[i][aadic[seq[i]]]   for   i in range(l)           ]
    plmdiag=np.array(plmdiag)
    plmseq=np.array(plmseq).reshape([l,1])

    midiag=[   np.diag(mi[i,i])     for i in range(l)]
    miseq=[ midiag[i][aadic[seq[i]]]   for   i in range(l)           ]
    midiag=np.array(midiag)
    #print(midiag.shape,'midiagshape')
    miseq=np.array(miseq).reshape([l,1])
    
    seq=np.array([aadic[aseq] for aseq in list(seq)])
    seq=np.eye(21)[seq]
    fea1d=np.concatenate([plmdiag,plmseq,midiag,miseq,seq],1)
    return fea1d


def cal_large_mi(msa,weight):
    #output:441*l*l
    ALPHA=21
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
def addmis(ccmfea,msafile,savefile):
    #import libcontact
    seq=open(msafile).readlines()[0].strip()
    num_aligns=len(open(msafile).readlines())

    if not os.path.isfile(savefile+'.weight'):
        exe=os.path.join(os.path.dirname(os.path.abspath(__file__)),'bin/colors_weight')
        cmd=exe+' '+msafile+' '+str(0.8)+' >'+savefile+'.weight'
        os.system(cmd)
    if not os.stat(savefile+'.weight').st_size==0:
        weight=np.genfromtxt(savefile+'.weight').flatten()
    else:
        print('can not get weight, using 1 as default')
        weight=np.ones(num_aligns)
    mi=cal_large_mi(read_msa(msafile),weight)   
    l=mi.shape[0]
    contact=np.zeros([2,l,l])
    for i in range(l):
        for j in range(i+1,l):
            contact[0,i,j]=np.sum(mi[i,j])
            contact[0,j,i]=contact[0,i,j]
            contact[1,i,j]=mi[i,j,aadic[seq[i]],aadic[seq[j]] ]
            contact[1,j,i]=mi[j,i,aadic[seq[j]],aadic[seq[i]] ]
    
    ccmfea=np.concatenate([ccmfea,contact])
    
    #ccmfea=np.concatenate([generate2d(l),ccmfea])
    #print('qqq\n',ccmfea[-2:,:,:])
    np.save(savefile+'.2d',ccmfea)
    fea1d=get1d(ccmfea,mi,seq)
    #print(fea1d)
    #print(fea1d.shape,generate1d(l).shape)
    #fea1d=np.concatenate([generate1d(l).transpose(),fea1d],1)
    np.save(savefile+'.1d',fea1d)
    #os.remove(savefile+'.ccm.npy') 

def d2fea(msafile,savefile):
    plm=computeccm_bin(msafile,savefile+'.ccm')[0]
    addmis(plm,msafile,savefile)


def evaluate(mina,maxa,bina,Lfactor,sml):
    #mina,maxa,bina='2','16','55'
    min_dis,max_dis,bin_num=float(mina),float(maxa),int(bina)
    dis_region=np.zeros(bin_num)
    for i in range(bin_num):
        dis_region[i]=min_dis+(i+0.5)*(max_dis-min_dis)*1.0/bin_num
    out=np.load('test/T0862-D1_'+bina+'_'+mina+'_'+maxa+'_ca.npy')
    contact_dis=np.genfromtxt('test/T0862-D1/T0862-D1.contact_ca')
    noncontact_=out[-1,:,:]
    L=noncontact_.shape[-1]
    mask=np.zeros([L,L])
    if sml=='S':
        for i in range(L):
            for j in range(L):
                if j-i>=6 and j-i<12:
                    mask[i,j]=1
    if sml=='M':
        for i in range(L):
            for j in range(L):
                if j-i>=12 and j-i<24:
                    mask[i,j]=1
    if sml=='L':
        for i in range(L):
            for j in range(L):
                if j-i>=24:
                    mask[i,j]=1
    noncontact=(1-noncontact_)*mask
    all_error=0
    all_true=0
    all_ree=0
    log_all_error=0
    ij=0
    log_all_ree=0   
    log_meanpre,meanpre,meand=0,0,0
    for i in range(int(Lfactor*L)):
        indexes=np.unravel_index(noncontact.argmax(), noncontact.shape)
        noncontact[indexes[0],indexes[1]]=0
        preout=out[1:-1,indexes[0],indexes[1]]/np.sum(out[1:-1,indexes[0],indexes[1]])
        pred=np.sum(preout*dis_region)
        all_error+=np.abs(pred-contact_dis[indexes[0],indexes[1]])
        all_ree+=pred
        ij+=np.abs(indexes[0]-indexes[1])
        if contact_dis[indexes[0],indexes[1]]<max_dis:
            all_true+=1
        meand+=contact_dis[indexes[0],indexes[1]]
        meanpre+=pred

        pred=np.exp(np.sum(preout*np.log(dis_region)))
        log_all_error+=np.abs(pred-contact_dis[indexes[0],indexes[1]])
        log_all_ree+=np.abs(pred-contact_dis[indexes[0],indexes[1]])/(pred+contact_dis[indexes[0],indexes[1]])*2
        log_meanpre+=pred        
    
    print( Lfactor,'     ',sml," %.2f (%.2f %d %d)" %(all_error/int(Lfactor*L),all_ree/int(Lfactor*L),round(ij/int(Lfactor*L)),all_true))

if __name__ == '__main__':
    print("topxL range err_L(dij,|i-j|,Np) ")
    #print(re_out,label)
    for sml in ['S','M','L']:
        for Lfactor in [0.5,2.0,5.0]:
            
            evaluate(sys.argv[1],sys.argv[2],sys.argv[3],Lfactor,sml)

