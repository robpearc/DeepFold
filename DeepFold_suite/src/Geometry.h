/*******************************************************************************************
**  
**  Functions for converting between torsion/cartesian space and for applying gradients to 
**  the given degrees of freedom
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include "CommonParameters.h"
#include "Operations.h"

using namespace std;

class Geometry{
public:
    void loadsgpos2( const char *filename, int ndim );

    void randomdih( double *phi, double *psi );

    bool tor2coord( point3f *decstr, int seqnum );

    bool tor2coordsg( point3f *decstr, int seqnum );

    bool tor2coordsg_all( point3f *decstr, int seqnum );

    bool setinitdecoy( point3f *decstr, int numseq );
    bool setinitdecoyfromfile( point3f *decstr, int numseq, const char *phi_file, 
                               const char *psi_file );

    bool apply_torsions( point3f *decstr, double *vars, int numseq );

    virtual ~Geometry();

private:
    double ****sgposdat;
    bool flagsgpos2;

};

Geometry::~Geometry()
{
    release4DArr( 6, 20, 180, sgposdat );
}

void Geometry::loadsgpos2( const char *filename, int ndim )
{
    FILE *file;
    file = fopen(filename,"rt");
    if( !file ){
        printf( "Error! 'Geometry::loadsgpos2', No sgposdat file %s\n", filename );
        flagsgpos2 = false;
    }
    int i,j,k;
    char oneline[300];
    sgposdat = new4DArr( 6, 20, 180, 180 );
    for(i=0;i<20;i++)
    {
        for(j=0;j<ndim;j++)
        {
            for(k=0;k<ndim;k++)
            {
                fgets(oneline,300,file);
                sscanf(oneline,"%lf %lf %lf %lf %lf %lf",&sgposdat[0][i][j][k],&sgposdat[1][i][j][k],
				&sgposdat[2][i][j][k],&sgposdat[3][i][j][k],&sgposdat[4][i][j][k],
				&sgposdat[5][i][j][k]);

            }
        }
    }
    flagsgpos2 = true;
    return;
}

void Geometry::randomdih( double *phi, double *psi )
{
    double ran=Random(); 
    if(ran<=0.135)
    {
        *phi=-140.0+360.0;
        *psi=153.0;
    }
    else if(ran>0.135 and ran<=0.29)
    {    
        *phi=-72.0+360.0;
        *psi=145.0;//+180.0;
    }
    else if(ran>0.29 and ran<=0.363)
    {
        *phi=-122.0+360.0;
        *psi=117.0;//+180.0;
    }
    else if(ran>0.363 and ran<=0.485)
    {
        *phi=-82.0+360.0;
        *psi=-14.0+360.0;
    }
    else if(ran>0.485 and ran<=0.982)
    {
        *phi=-61.0+360.0;
	*psi=-41.0+360.0;
    }
    else
    {
        *phi=57.0;//+180.0;
        *psi=39.0;//+180.0;
    }
    return;
}

bool Geometry::tor2coord( point3f *decstr, int seqnum )
{
    int i;
    bool flagts;
    bool flagwhole=true;
    double len_n_ca=1.460;
    double len_ca_c=1.525;
    double len_c_n=1.338;
    double ang_n_ca_c=111.008;
    double ang_ca_c_n=116.617;
    double ang_c_n_ca=121.614;

    //if(decstr[0].len[1]<0) decstr[0].len[1]=//float(lennca);
    //if(decstr[0].len[2]<0) decstr[0].len[2]=//float(lencac);
    //if(decstr[0].ang[2]<0) decstr[0].ang[2]=//float(angncac);

    if(decstr[0].len_n_ca<0) decstr[0].len_n_ca=len_n_ca;//float(lennca);
    if(decstr[0].len_ca_c<0) decstr[0].len_ca_c=len_ca_c;//float(lencac);
    if(decstr[0].ang_n_ca_c<0) decstr[0].ang_n_ca_c=ang_n_ca_c;//float(angncac);
    decstr[0].ptn.x=0;
    decstr[0].ptn.y=0;
    decstr[0].ptn.z=0;
    decstr[0].x=decstr[0].len_n_ca;
    decstr[0].y=0;
    decstr[0].z=0;
    decstr[0].ptc.x=decstr[0].len_n_ca-decstr[0].len_ca_c*cos(decstr[0].ang_n_ca_c*raddeg);
    decstr[0].ptc.y=decstr[0].len_ca_c*sin(decstr[0].ang_n_ca_c*raddeg);
    decstr[0].ptc.z=0;
    //if(decstr[0].tor[0]>=0 && decstr[0].tor[2]>=0)
    /*
    if(decstr[0].psi.item<double>()>=0 && decstr[0].phi.item<double>()>=0)
    {
        //flagts=tor2pos22(zero_tensor,zero_tensor,zero_tensor,
        //                 len_n_ca,zero_tensor,zero_tensor,
        //                 len_n_ca+len_ca_c*at::sin(ang_n_ca_c*raddeg),-len_ca_c*at::cos(ang_n_ca_c*raddeg),zero_tensor,
        //                 decstr[0].psi*raddeg,len_c_n,ang_ca_c_n*raddeg,
	//		 &pn.x,&pn.y,&pn.z);
        //if(!flagts)
        //{
        //    flagwhole=false;
            //printf("wrong front coordinates n %d\n",0);
        //}
        //decstr[0].ptn.x=pn.x;
        //decstr[0].ptn.y=pn.y;
        //decstr[0].ptn.z=pn.z;
        if(decstr[0].omega.item<double>()<0) decstr[0].omega=tensor_180;
        flagts=tor2pos22(len_n_ca,zero_tensor,zero_tensor,
                         len_n_ca+len_ca_c*at::sin(ang_n_ca_c*raddeg),-len_ca_c*at::cos(ang_n_ca_c*raddeg),zero_tensor,
                         decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,
                         decstr[0].omega*raddeg,len_n_ca,ang_c_n_ca*raddeg,
                         &pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates ca %d\n",0);
        }
        decstr[0].x=pt.x;
        decstr[0].y=pt.y;
        decstr[0].z=pt.z;
        flagts=tor2pos22(len_n_ca+len_ca_c*at::sin(ang_n_ca_c*raddeg),-len_ca_c*at::cos(ang_n_ca_c*raddeg),zero_tensor,
                         decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,
                         decstr[0].x,decstr[0].y,decstr[0].z,
                         decstr[0].phi*raddeg,len_ca_c,ang_n_ca_c*raddeg,
			 &pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates c %d\n",0);
        }
        decstr[0].ptc.x=pc.x;
        decstr[0].ptc.y=pc.y;
        decstr[0].ptc.z=pc.z;
    }
    */
    for(i=1;i<seqnum;i++)
    {
        point3s pt,pn,pc;

        if(decstr[i].psi<0)
        {
            decstr[i].psi=120.0;
        }
        if(decstr[i].omega<0) decstr[i].omega=180.0;
        if(decstr[i].phi<0) decstr[i].phi=290.0;
        if(decstr[i].len_c_n<0) decstr[i].len_c_n=len_c_n;
        if(decstr[i].len_n_ca<0) decstr[i].len_n_ca=len_n_ca;
        if(decstr[i].len_ca_c<0) decstr[i].len_ca_c=len_ca_c;
        if(decstr[i].ang_ca_c_n<0) decstr[i].ang_ca_c_n=ang_ca_c_n;
        if(decstr[i].ang_c_n_ca<0) decstr[i].ang_c_n_ca=ang_c_n_ca;
        if(decstr[i].ang_n_ca_c<0) decstr[i].ang_n_ca_c=ang_n_ca_c;
        //0 1 2 n ca c
        //original phi i-1 psi i-1 omega i
        flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                         decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                         decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                         decstr[i].psi*raddeg,decstr[i].len_c_n,
                         decstr[i].ang_ca_c_n*raddeg,&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
            flagwhole=false;
            //printf("wrong front coordinates n %d %f %f %f\n",
            //       i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
        }
        decstr[i].ptn.x=pn.x;
        decstr[i].ptn.y=pn.y;
        decstr[i].ptn.z=pn.z;
        flagts=tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                         decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                         decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         decstr[i].omega*raddeg,decstr[i].len_n_ca,
                         decstr[i].ang_c_n_ca*raddeg,&pt.x,&pt.y,&pt.z);
        if(!flagts)
        {
            flagwhole=false;
            //printf("wrong front coordinates ca %d %f %f %f\n",
            //       i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
        }
        decstr[i].x=pt.x;
        decstr[i].y=pt.y;
        decstr[i].z=pt.z;
        flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                         decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         decstr[i].x,decstr[i].y,decstr[i].z,
                         decstr[i].phi*raddeg,decstr[i].len_ca_c,
                         decstr[i].ang_n_ca_c*raddeg,&pc.x,&pc.y,&pc.z);
        if(!flagts)
        {
            flagwhole=false;
            //printf("wrong front coordinates c %d %f %f %f\n",
            //       i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
        }
        decstr[i].ptc.x=pc.x;
        decstr[i].ptc.y=pc.y;
        decstr[i].ptc.z=pc.z;
    }
    return flagwhole;
}

bool Geometry::tor2coordsg( point3f *decstr, int seqnum )
{
    int i;
    int tind;

    bool flagts;
    bool flagwhole=true;

    for(i=1;i<seqnum;i++)
    {
        point3s pn1,pn2;
        //atom o
        double ang1=179.6715*raddeg;
        double len1=1.229;
        double len2=2.0961;
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                         decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                         ang1,len1,len2,&pn1.x,&pn1.y,&pn1.z);
        decstr[i-1].pto.x=pn1.x;
        decstr[i-1].pto.y=pn1.y;
        decstr[i-1].pto.z=pn1.z;
        //atom h
        double ang2=179.8174*raddeg;
        double len3=0.987;
        double len4=2.0814;
        flagts=tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                         decstr[i].x,decstr[i].y,decstr[i].z,
                         decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         ang2,len3,len4,&pn2.x,&pn2.y,&pn2.z);//0.9919f,2.0574f
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates d %d\n",i);
        }
        decstr[i].pth.x=pn2.x;
        decstr[i].pth.y=pn2.y;
        decstr[i].pth.z=pn2.z;
    }
    point3s pn3,pn4;
    //atom o
    double ang1=0.0;
    double len1=1.2439;
    double len2=2.0855;
    i=seqnum;
    flagts=tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                     decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                     decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                     ang1,len1,len2,&pn3.x,&pn3.y,&pn3.z);
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates o %d\n",i);
    }
    decstr[i-1].pto.x=pn3.x;
    decstr[i-1].pto.y=pn3.y;
    decstr[i-1].pto.z=pn3.z;
    //atom h 
    double ang2=60.0*raddeg;
    double len3=0.987;
    double len4=2.0306;
    i=0;
    flagts=tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                     decstr[i].x,decstr[i].y,decstr[i].z,
                     decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                     ang2,len3,len4,&pn4.x,&pn4.y,&pn4.z);//0.9972f
    if(!flagts)
    {
        flagwhole=false;
        printf("wrong front coordinates %d\n",i);
    }
    decstr[i].pth.x=pn4.x;
    decstr[i].pth.y=pn4.y;
    decstr[i].pth.z=pn4.z;
    //atom cb new
    for(i=0;i<seqnum;i++)
    {
        point3s pn5;
        tind=aminoid(decstr[i].aaa);
        if(tind>19) tind=19;
        double ang3=cbsta[tind][2]*raddeg;
        double len5=cbsta[tind][0];
        double len6=cbsta[tind][1];
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                         decstr[i].x,decstr[i].y,decstr[i].z,
                         ang3,len5,len6,&pn5.x,&pn5.y,&pn5.z);
        if(!flagts)
        {
            flagwhole=false;
            printf("wrong front coordinates cb2 %d\n",i);
        }
        decstr[i].ptb.x=pn5.x;
        decstr[i].ptb.y=pn5.y;
        decstr[i].ptb.z=pn5.z;
        if(decstr[i].aaa=='G')
        {
            decstr[i].ptb.x=decstr[i].x;
            decstr[i].ptb.y=decstr[i].y;
            decstr[i].ptb.z=decstr[i].z;
        }
    }
    return flagwhole;
}

bool Geometry::tor2coordsg_all( point3f *decstr, int seqnum )
{
    bool flagwhole=tor2coordsg( decstr, seqnum );
    int i;
    int tind;
    bool flagts;

    //atom sg
    int ti,tj;
    int cutnum=72;
    double delta=5.0;
    for ( i=0; i<seqnum; i++ )
    {
        point3s pn1,pn2;
        if( decstr[i].aaa=='G' )
        {   
            decstr[i].ptsg.x=decstr[i].x;
            decstr[i].ptsg.y=decstr[i].y;
            decstr[i].ptsg.z=decstr[i].z;
            continue;
        }
        tind=aminoid(decstr[i].aaa);
        if( tind>19 ) tind=19;
        if( i<seqnum-1 )
        {
            ti=int(decstr[i].phi/delta);
            tj=int(decstr[i+1].psi/delta);
            if ( ti>=0 && ti<cutnum && tj>=0 && tj<cutnum )
            {
                double ang1=sgposdat[2][tind][ti][tj];
                double len1=sgposdat[0][tind][ti][tj];
                double len2=sgposdat[1][tind][ti][tj];
                flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                                 decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                                 decstr[i].x,decstr[i].y,decstr[i].z,
                                 ang1,len1,len2,&pn1.x,&pn1.y,&pn1.z);
                if ( !flagts )
                {
                    flagwhole=false;
                    printf("wrong front coordinates full sg %d\n",i);
                }
                decstr[i].ptsg.x=pn1.x;
                decstr[i].ptsg.y=pn1.y;
                decstr[i].ptsg.z=pn1.z;
                continue;
            }
        }
        double ang1=sglatavg[tind][2];
        double len1=sglatavg[tind][0];
        double len2=sglatavg[tind][1];
        flagts=tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                         decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                         decstr[i].x,decstr[i].y,decstr[i].z,
                         ang1*raddeg,len1,len2,&pn2.x,&pn2.y,&pn2.z);
        if ( !flagts )
        {
            flagwhole=false;
            printf("wrong front coordinates sg %d\n",i);
        }
        decstr[i].ptsg.x=pn2.x;
        decstr[i].ptsg.y=pn2.y;
        decstr[i].ptsg.z=pn2.z;
    }
    //atom ct
    for ( i=0; i<seqnum; i++ )
    {
        point3s pn3;
        tind=aminoid(decstr[i].aaa);
        pn3.x=0;pn3.y=0;pn3.z=0;
        pn3.x+=decstr[i].ptn.x;
        pn3.y+=decstr[i].ptn.y;
        pn3.z+=decstr[i].ptn.z;
        pn3.x+=decstr[i].x;
        pn3.y+=decstr[i].y;
        pn3.z+=decstr[i].z;
        pn3.x+=decstr[i].ptc.x;
        pn3.y+=decstr[i].ptc.y;
        pn3.z+=decstr[i].ptc.z;
        pn3.x+=decstr[i].pto.x;
        pn3.y+=decstr[i].pto.y;
        pn3.z+=decstr[i].pto.z;
        pn3.x+=decstr[i].ptsg.x*sgatomnum[tind];
        pn3.y+=decstr[i].ptsg.y*sgatomnum[tind];
        pn3.z+=decstr[i].ptsg.z*sgatomnum[tind];
        decstr[i].ptg.x=pn3.x/double(sgatomnum[tind]+4.0);
        decstr[i].ptg.y=pn3.y/double(sgatomnum[tind]+4.0);
        decstr[i].ptg.z=pn3.z/double(sgatomnum[tind]+4.0);
    }
    return flagwhole;
}

bool Geometry::setinitdecoyfromfile( point3f *decstr, int numseq, const char *phi_file, const char *psi_file )
{
    for ( int i=0; i<numseq; i++ )
    {
        decstr[i].len_n_ca=1.460;//len_n_ca;//float(lennca);
        decstr[i].len_ca_c=1.525;//len_ca_c;//float(lencac);
        decstr[i].len_c_n=1.338;//len_c_n;
        decstr[i].ang_n_ca_c=111.008;//ang_n_ca_c;//float(angncac);
        decstr[i].ang_ca_c_n=116.617;//ang_ca_c_n;//float(angncac);
        decstr[i].ang_c_n_ca=121.614;//ang_c_n_ca;//float(angncac);
	decstr[i].omega=180.0;//tensor_180;
    }

    FILE *file;
    if ( ( file=fopen(phi_file,"rt") )!=NULL )
    {
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,file))
        {
            int res;
            double val=0.0;
            sscanf(line,"%d %lf",&res,&val);
	    //if(res==1) continue;
            if(val<0) val+=360.0;
            decstr[res-1].phi=val;
        }
    }
    fclose(file);
    decstr[0].psi=0.0;
    if ( ( file=fopen(psi_file,"rt") )!=NULL )
    {
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,file))
        {
            int res;
            double val=0.0;
            sscanf(line,"%d %lf",&res,&val);
	    if(res==numseq) continue;
            if(val<0) val+=360.0;
            decstr[res].psi=val;
        }
    }
    fclose(file);

    tor2coord( decstr, numseq );
    return tor2coordsg_all( decstr, numseq );
}



bool Geometry::setinitdecoy( point3f *decstr, int numseq )
{
    for ( int i=0; i<numseq; i++ )
    {
        decstr[i].len_n_ca=1.460;//len_n_ca;//float(lennca);
        decstr[i].len_ca_c=1.525;//len_ca_c;//float(lencac);
        decstr[i].len_c_n=1.338;//len_c_n;
        decstr[i].ang_n_ca_c=111.008;//ang_n_ca_c;//float(angncac);
        decstr[i].ang_ca_c_n=116.617;//ang_ca_c_n;//float(angncac);
        decstr[i].ang_c_n_ca=121.614;//ang_c_n_ca;//float(angncac);
        decstr[i].omega=180.0;//tensor_180;
    }

    for ( int i=0; i<numseq-1; i++ )
    {
        double phi,psi;
        randomdih(&phi,&psi);
        decstr[i].phi=phi;
        decstr[i+1].psi=psi;
    }
    double phi,psi;
    randomdih(&phi,&psi);
    decstr[numseq-1].phi=phi;

    tor2coord( decstr, numseq );
    return tor2coordsg_all( decstr, numseq );
}


bool Geometry::apply_torsions( point3f *decstr, double *vars, int numseq )
{
    for ( int i=0; i<numseq; i++ )
    {
        decstr[i].len_n_ca=1.460;//len_n_ca;//float(lennca);
        decstr[i].len_ca_c=1.525;//len_ca_c;//float(lencac);
        decstr[i].len_c_n=1.338;//len_c_n;
        decstr[i].ang_n_ca_c=111.008;//ang_n_ca_c;//float(angncac);
        decstr[i].ang_ca_c_n=116.617;//ang_ca_c_n;//float(angncac);
        decstr[i].ang_c_n_ca=121.614;//ang_c_n_ca;//float(angncac);
        decstr[i].omega=180.0;//tensor_180;
    }

    //dummy torsions
    decstr[0].psi = 0.0;
    decstr[0].phi = 0.0;
    for ( int i=1; i<numseq; i++ )
    {
        if ( vars[(i-1)*2]>=360.0 )
        {
            while ( vars[(i-1)*2]>=360.0 )
            {
                vars[(i-1)*2]-=360.0;
            }
        }
        else if ( vars[(i-1)*2]<0.0 )
        {
            while ( vars[(i-1)*2]<0.0 )
            {
                vars[(i-1)*2]+=360.0;
            }
        }

        if ( vars[(i-1)*2+1]>=360.0 )
        {
            while ( vars[(i-1)*2+1]>=360.0 )
            {
                vars[(i-1)*2+1]-=360.0;
            }
        }
        else if ( vars[(i-1)*2+1]<0.0 )
        {
            while ( vars[(i-1)*2+1]<0.0 )
            {
                vars[(i-1)*2+1]+=360.0;
            }
        }

        decstr[i].psi=vars[(i-1)*2];
        decstr[i].phi=vars[(i-1)*2+1];
    }
    tor2coord( decstr, numseq );
    return tor2coordsg_all( decstr, numseq );
}

#endif
