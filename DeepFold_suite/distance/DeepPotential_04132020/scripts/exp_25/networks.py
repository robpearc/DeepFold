import os,sys
import torch
import torch.nn as nn
import torch.nn.functional as F



CUTOFF=0.8
ALPHA=21
USECUDA=torch.cuda.is_available()


def conv3x3(in_planes, out_planes,kernel_size,dialation):
    """3x3 convolution with padding"""
    return nn.Conv2d(in_planes, out_planes, kernel_size=kernel_size,dilation=dialation,
                     padding=dialation*(kernel_size-1)//2, bias=False)

def conv3(in_planes, out_planes,kernel_size,dialation):
    """3 convolution with padding"""
    return nn.Conv1d(in_planes, out_planes, kernel_size=kernel_size,dilation=dialation,
                     padding=dialation*(kernel_size-1)//2, bias=False)

class BasicBlock(nn.Module):
    def __init__(self, planes,kernel_sizes,dialations):
        super(BasicBlock,self).__init__()
        self.conv1 = conv3x3(planes, planes,kernel_sizes,dialations)
        self.bn1 = nn.InstanceNorm2d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv3x3(planes, planes,kernel_sizes,dialations)
        self.bn2 = nn.InstanceNorm2d(planes)
        self.droprate = 0.2    
    def forward(self,x):
        residual = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        if self.droprate > 0:
            out = F.dropout(out, p=self.droprate, training=self.training)
        out = self.conv2(out)
        out = self.bn2(out)
        out += residual
        out = self.relu(out)
        return out            
class BasicBlock1d(nn.Module):
    expansion = 1
    def __init__(self, planes,kernel_sizes,dialations):
        super(BasicBlock1d, self).__init__()
        self.conv1 = conv3(planes, planes,kernel_sizes,dialations)
        self.bn1 = nn.InstanceNorm1d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv3(planes, planes,kernel_sizes,dialations)
        self.bn2 = nn.InstanceNorm1d(planes)
        self.droprate = 0.2
    def forward(self, x):
        residual = x
        out = self.conv1(x)
        if self.droprate > 0:
            out = F.dropout(out, p=self.droprate, training=self.training)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.conv2(out)
        if self.droprate > 0:
            out = F.dropout(out, p=self.droprate, training=self.training)
        out = self.bn2(out)
        out += residual
        out = self.relu(out)
        return out

def d12d(x,inpane,L,batch_size):
    #x: B,F,L
    x1=torch.unsqueeze(x,3)
    x1=x1.repeat(1,1,1,L)
    x2=torch.unsqueeze(x,2)
    x2=x2.repeat(1,1,L,1)
    #x=x.expand(batch_size,inpane,L,L)
    x=torch.cat([x1,x2],1)
    return x

class preNet(nn.Module):
    def __init__(self,featurenum2d,featurenum1d,dis_dim,theta_dim,phi_dim,omg_dim,hbp_dim,do_attention):
        super(preNet,self).__init__()
        self.dim_2d=featurenum2d
        self.dim_1d=featurenum1d
        self.hidden_channel=64
        self.hidden_channel1d=32
        self.kernel_sizes=[3,5]*80
        self.dialations=[1,2,4,8,16]*40
        self.blocks1=[10]
        self.blocks2=[40]
        self.blocks1d=[10]
        self.blocks_in_att=[8,4]
        #self.blocks=[2,2,2,2]
        #self.blocks_in_att=[1,1]
        
        if True:
            self.conv_tran1=nn.Conv2d(featurenum2d,self.hidden_channel,kernel_size=1, bias=False)
            self.bn_tran1=nn.InstanceNorm2d(self.hidden_channel)

            self.conv_tran2=nn.Conv2d(self.hidden_channel+self.hidden_channel1d*2,self.hidden_channel,kernel_size=1, bias=False)
            self.bn_tran2=nn.InstanceNorm2d(self.hidden_channel)

            self.conv_tran1d=nn.Conv1d(featurenum1d,self.hidden_channel1d,kernel_size=1, bias=False)
            self.bn_tran1d=nn.InstanceNorm1d(self.hidden_channel1d)

            self.res1=self._make_layer(BasicBlock,self.blocks1[0],self.kernel_sizes[:],self.dialations[:])
            self.res2=self._make_layer(BasicBlock,self.blocks2[0],self.kernel_sizes[:],self.dialations[:])
            self.res1d=self._make_layer1d(BasicBlock1d,self.blocks1d[0],self.kernel_sizes[:],self.dialations[:])


        if True:
            self.conv_cb=nn.Conv2d(self.hidden_channel,dis_dim,3,padding=1)
            self.conv_omg=nn.Conv2d(self.hidden_channel,omg_dim,1,padding=0)
            self.conv_hbp=nn.Conv2d(self.hidden_channel,hbp_dim,3,padding=1)
            # un sys
            self.conv_theta=nn.Conv2d(self.hidden_channel,theta_dim,1,padding=0)
            self.conv_phi=nn.Conv2d(self.hidden_channel,phi_dim,1,padding=0)


                   
        self.relu = nn.ReLU(inplace=True)
        self.logsoft=nn.LogSoftmax(1)
    def _make_layer(self,block,num_blocks,kernel_sizes,dialations):
        layers=[block(self.hidden_channel,kernel_sizes[i],dialations[i]) for i in range(num_blocks)]
        return nn.Sequential(*layers)
    def _make_layer1d(self,block,num_blocks,kernel_sizes,dialations):
        layers=[block(self.hidden_channel1d,kernel_sizes[i],dialations[i]) for i in range(num_blocks)]
        return nn.Sequential(*layers)
    def forward(self,x,f1d,intra=True):
        L=f1d.shape[-1]
        if intra:# most cases
            f1d=self.conv_tran1d(f1d.unsqueeze(0))
            f1d=self.bn_tran1d(f1d)
            f1d=self.relu(f1d)
            f1d=self.res1d(f1d)

            x=self.conv_tran1(x.unsqueeze(0))
            x=self.bn_tran1(x)    
            x=self.relu(x)

            x=self.res1(x) 
            #print(f1d.shape,self.dim_1d,L,1)      
            x=torch.cat([d12d(f1d,self.hidden_channel1d,L,1),x],1)
            x=self.conv_tran2(x)
            x=self.bn_tran2(x)
            
            x=self.relu(x)
            x=self.res2(x)

            
            ptheta=self.logsoft(self.conv_theta(x)  )   
            pphi=self.logsoft(self.conv_phi(x) )         
            #x=(x.transpose(-1,-2)+x)*0.5
            #print('x shape',x.shape)
            
            pcb=self.logsoft(self.conv_cb(x))
            pomg=self.logsoft(self.conv_omg(x))

            phb=self.logsoft(self.conv_hbp(x))
            pcb=(pcb.transpose(-1,-2)+pcb)*0.5
            pomg=(pomg.transpose(-1,-2)+pomg)*0.5
            phb=(phb.transpose(-1,-2)+phb)*0.5
            return pcb,ptheta,pphi,pomg,phb




