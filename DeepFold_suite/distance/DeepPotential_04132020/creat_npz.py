#!/nfs/amino-home/liyangum/miniconda3/bin/python 
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 16:30:45 2017

@author: lee
"""

import sys,os,re
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import torch
import torch.nn as nn
from torch.autograd import Variable
torch.manual_seed(6)
import random
random.seed(6)
import file_util as file_utils
print('cuda is ready? :',torch.cuda.is_available())
#import matplotlib
#matplotlib.use('agg')
#from matplotlib import pyplot as plt 
import features
import pickle
import importlib.util
develop=False
pre_list=       ['7','10','11','15',   '20',   '28']
ensemble_weight={'7':2.5,'10':2.0,'11':2.0,'15':2.0,'19':0,  '20':2,'22':2,'23':2,'24':2,'25':2,  '28':2,'29':2,'30':2,'31':2,'32':2}
plusfeatures=['10','11','15','19','20','28']
Netconfig={
    '7':[36+2,25,13,25,3,False],
    '10':[36+2,25,13,25,3,False],
    '11':[19,19,19,38,38,False],
    '15':[38,25,13,25,3,False],
    '19':[38,25,13,25,3,False],
    
    '20':[38,25,13,25,3,False],
    '22':[18,25,13,25,3,False],
    '23':[24,25,13,25,3,False],
    '24':[30,25,13,25,3,False],
    '25':[36,25,13,25,3,False],

    '28':[19,19,19,38,38,False],
    '29':[18,25,13,25,3,False],
    '30':[24,25,13,25,3,False],
    '31':[30,25,13,25,3,False],
    '32':[36,25,13,25,3,False]

}

def npytotxt(a_,outfile,sym=True):

    a=a_[0]
    if develop:
        np.save(outfile,a)
        return 0
    L=a.shape[-1]
    scores=''
    for i in range(0,L):
        for j in range(i+1,L):
            strs=[str(ap) for ap in a[:,i,j]]
            score=' '.join(strs)
            score=str(i+1)+' '+str(j+1)+' '+score+'\n'
            scores+=score
            if not sym:
                strs=[str(ap) for ap in a[:,j,i]]
                score=' '.join(strs)
                score=str(j+1)+' '+str(i+1)+' '+score+'\n'
                scores+=score
    wfile=open(outfile,'w')
    wfile.write(scores)
    wfile.close()

def return_add(re1,re2):
    if re1==():
        return re2
    re_sum=[]
    for an1,an2 in zip(re1,re2):
        re_sum.append(an1+an2)
    return tuple(re_sum)   

def return_div(re1,num):
    return [an1*1.0/num for an1 in re1]
      

def return_exp(re1):
    re_exp=[]
    for an1 in re1:
        re_exp.append( np.exp(an1.cpu().data. numpy()) )
    return tuple(re_exp)

def out_npz(cb,pomg,ptheta,pphi,savefile_):
    L=cb.shape[-1]
    pd=np.zeros([L,L,37])
    po=np.zeros([L,L,25])
    pt=np.zeros([L,L,25])
    pp=np.zeros([L,L,13])

    pd=np.concatenate([  np.swapaxes(cb[0],0,2)[:,:,[-1]] , np.swapaxes(cb[0],0,2)[:,:,1:-1]               ],axis=-1)
    po=np.concatenate([  np.swapaxes(pomg[0],0,2)[:,:,[-1]], np.swapaxes(pomg[0],0,2)[:,:,:-1]                ],axis=-1)
    pt=np.concatenate([   np.swapaxes(ptheta[0],0,2)[:,:,[-1]],  np.swapaxes(ptheta[0],0,2)[:,:,:-1]                ],axis=-1)
    pp=np.concatenate([  np.swapaxes(pphi[0],0,2)[:,:,[-1]] , np.swapaxes(pphi[0],0,2)[:,:,:-1]               ],axis=-1)
    contacts = {'pd':pd, 'po':po, 'pt':pt, 'pp':pp}
    np.savez_compressed(savefile_, dist=contacts['pd'], omega=contacts['po'], theta=contacts['pt'], phi=contacts['pp'])

def prepare_longdistance(msafile,savefile_):
    if msafile.endswith('.a3m'):
        fea_1d_,fea_2d_=features.collect_allfeatures(msafile,savefile_)
    else:
        fea_1d_,fea_2d_=features.collect_allfeatures_old(msafile,savefile_)
    re_dic={}
    for anpre in pre_list:
        re_dic[anpre]=()
    L=fea_1d_.shape[-1]
    dim2=fea_2d_.shape[0]
    dim1=fea_1d_.shape[0]
    weight=np.genfromtxt(savefile_+'.weight').sum()
    adddims=np.zeros([4,L])
    adddims[0,0]=1
    adddims[1,-1]=1
    #adddims[2,1:-1]=1
    adddims[3,:]=weight
    for acon in pre_list:
        print('predicting',acon)
        model_dirs=os.path.join(os.path.dirname(__file__),'models','exp_'+acon)
        models=os.listdir(model_dirs)
        models=[amodel for amodel in models if not amodel.endswith('.opt')]
        spec = importlib.util.spec_from_file_location("preNet", os.path.join(os.path.dirname(__file__),'scripts','exp_'+acon,'networks.py'))
        foo = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(foo)
        for model_file in models:
            if acon in plusfeatures:
                #fea_1d=np.concatenate([adddims,fea_1d],0)           
                model=foo.preNet(dim2,dim1+4,*Netconfig[acon])
            else:
                model=foo.preNet(dim2,dim1,*Netconfig[acon])

            model.cpu()
            model_dict=torch.load(os.path.join(model_dirs,model_file), map_location=lambda storage, loc: storage)
            model.load_state_dict(model_dict)
            model.eval()

            fea_1d,fea_2d=fea_1d_*1.0,fea_2d_*1.0
            if acon in plusfeatures:
                fea_1d=np.concatenate([adddims,fea_1d],0) 
            fea_1d,fea_2d=torch.FloatTensor(np.array(fea_1d)),torch.FloatTensor(np.array(fea_2d))
            with torch.no_grad():
                anre=return_exp  (model(fea_2d,fea_1d))
                re_dic[acon]= return_add(re_dic[acon],anre)
        re_dic[acon] = return_div(re_dic[acon],len(models))
    pcb_20=(re_dic['7'][0]*ensemble_weight['7'] + re_dic['10'][0]*ensemble_weight['10'] + re_dic['11'][-1]*ensemble_weight['11'] +\
            re_dic['15'][0]*ensemble_weight['15'] + re_dic['28'][-1]*ensemble_weight['28'] + re_dic['20'][0]*ensemble_weight['20']) / \
                (ensemble_weight['7'] + ensemble_weight['10'] + ensemble_weight['11'] + ensemble_weight['15']+ensemble_weight['19']+ensemble_weight['20']+ensemble_weight['28'])
    
    pomg_20=(re_dic['7'][3]*ensemble_weight['7'] + re_dic['10'][3]*ensemble_weight['10']  +\
            re_dic['15'][3]*ensemble_weight['15']  + re_dic['20'][4]*ensemble_weight['20']) / \
                (ensemble_weight['7'] + ensemble_weight['10']  + ensemble_weight['15']+ensemble_weight['19']+ensemble_weight['20'])
    
    pphi_20=(re_dic['7'][2]*ensemble_weight['7'] + re_dic['10'][2]*ensemble_weight['10']  +\
            re_dic['15'][2]*ensemble_weight['15'] + re_dic['20'][3]*ensemble_weight['20']) / \
                (ensemble_weight['7'] + ensemble_weight['10']  + ensemble_weight['15']+ensemble_weight['19']+ensemble_weight['20'])      

    ptheta_20=(re_dic['7'][1]*ensemble_weight['7'] + re_dic['10'][1]*ensemble_weight['10']  +\
            re_dic['15'][1]*ensemble_weight['15']  + re_dic['20'][2]*ensemble_weight['20']) / \
                (ensemble_weight['7'] + ensemble_weight['10']  + ensemble_weight['15']+ensemble_weight['19']+ensemble_weight['20'])      

    out_npz(pcb_20,pomg_20,ptheta_20,pphi_20,savefile_+'_20.npz')

    sys.stdout.write('\ndone.')




if __name__ == '__main__':
    alnfile=sys.argv[1]
    savefile=sys.argv[2]
    #plot_native('test/T0862-D1.zhang','test/T0862-D1.contact_ca',savefile)
    prepare_longdistance(alnfile,savefile)
    #testdis('test/T0862-D1/T0862-D1.pdb')


         
            
    
    
        
