#!/usr/bin/env python
docstring='''
npz2rr ResTriplet3_20.npz ResTriplet3.rr
    Convert trRosetta format distance prediction file to contact map
'''

import numpy as np
import sys

#min_dis,max_dis,num_bin=2.,20.,36

def npz2rr(infile,outfile):
    npz        = np.load(infile)
    cscore_list= []
    dist_mat   = np.concatenate((npz['dist' ][:,:,1:],npz['dist' ][:,:,:1]),axis=2)
    L          = len(dist_mat)
    for i in range(L):
        for j in range(i+1,L):
            cscore=dist_mat[i,j,:12].sum()
            cscore_list.append((cscore,i,j))
    cscore_list=sorted(cscore_list,reverse=True)
    txt=''
    for cscore,i,j in cscore_list:
        txt+=str(i+1)+' '+str(j+1)+' '+str(cscore)+'\n'
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return cscore_list

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    npz2rr(sys.argv[1],sys.argv[2])
