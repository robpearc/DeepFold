/*******************************************************************************************
 *   ______               ______    _     _ 
 *   |  _  \              |  ___|  | |   | |
 *   | | | |___  ___ _ __ | |_ ___ | | __| |
 *   | | | / _ \/ _ \ '_ \|  _/ _ \| |/ _` |
 *   | |/ /  __/  __/ |_) | || (_) | | (_| |
 *   |___/ \___|\___| .__/\_| \___/|_|\__,_|
 *                  | |                     
 *                  |_|                      
 *
 *  This program was written by Robin Pearce at the Zhang Lab
 *  Department of Computational Medicine and Bioinformatics 
 *  University of Michigan 
 *           
 *  Please report bugs and questions to robpearc@umich.edu
 *
*******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "CubicSpline.h"
#include "ParseSeq.h"
#include "Geometry.h"
#include "Energy.h"
#include "LBFGS_optimizer.h" 
#include "CommonParameters.h"
#include "Operations.h"

using namespace std;

Energy energyFunction;


void show_interface(){
    printf(
      "#########################################################################\n"
      "                                                                         \n"
      "          ______               ______    _     _                         \n"
      "          |  _  \\              |  ___|  | |   | |                        \n"
      "          | | | |___  ___ _ __ | |_ ___ | | __| |                        \n"
      "          | | | / _ \\/ _ \\ '_ \\|  _/ _ \\| |/ _` |                        \n"
      "          | |/ /  __/  __/ |_) | || (_) | | (_| |                        \n"
      "          |___/ \\___|\\___| .__/\\_| \\___/|_|\\__,_|                        \n"
      "                         | |                                             \n"
      "                         |_|                                             \n"
      "                                                                         \n"
      "   A program for ab initio protein structure prediction guided by        \n"
      "                 potentials from deep learning                           \n"
      "\n\n"
      "  Written by Robin Pearce at the Yang Zhang Lab\n"
      "  Dept. of Computational Medicine & Bioinformatics\n"
      "  University of Michigan\n"
      "  For questions email robpearc@umich.edu\n"
      "#########################################################################\n");
}


int main(int argc, char** argv){
    
    string datadir = argv[1]; 
    string libdir  = argv[2];

    show_interface();

    //------- read sequence and secondary structure ---------->
    ParseSeq ps,ps2;
    Geometry geo;
    string sgposfile = libdir+"/newsgdistriaanew72.txt";
    geo.loadsgpos2(sgposfile.c_str(), 72);
    
    string seqfile = datadir+"/seq.txt";
    string ssfile  = datadir+"/seq.dat.ss";
    if (!ps.loadseq(seqfile.c_str())) return 1;
    int numseq=ps.seqnum;
    ps.loadss2(ssfile.c_str());
    ps2.loadss2(ssfile.c_str());
    ps2.modifyss2();
    ps.genesse();
    ps2.genesse2();
    int longestH = ps2.getlongestH();

    string restraint_path = datadir+"/DeepPotential_20.npz";
    string weight_file1 = libdir+"/weight_stg1.txt";
    string weight_file2 = libdir+"/weight_stg2.txt";
    string weight_file3 = libdir+"/weight_stg3.txt";
    string ca_dist_file = datadir+"/DeepPotential_pca_20.txt";
    string cb_cont_file = datadir+"/DeepPotential_cb_contact.txt";
    string ca_cont_file = datadir+"/DeepPotential_ca_contact.txt";
    string phi_file=datadir+"/phi.txt";
    string psi_file=datadir+"/psi.txt";

    // Set up energy function
    energyFunction.load_files( longestH, numseq, datadir, libdir, weight_file1, restraint_path, 
                               ca_dist_file, cb_cont_file, ca_cont_file );

    // Generate initial decoy
    point3f *mcdecstr=new point3f[numseq];
    double *vars = new double[2*(numseq-1)];
    for ( int i=0; i<numseq; i++ )
    {
        mcdecstr[i].ss2=ps.ss2[i].ss;
        mcdecstr[i].aaa=ps.seqdata[i];
        mcdecstr[i].iaa=(unsigned char)(aminoid(ps.seqdata[i]));
        if(mcdecstr[i].iaa>19) mcdecstr[i].iaa=5;
    }
    geo.setinitdecoyfromfile( mcdecstr, numseq, phi_file.c_str(), psi_file.c_str() );
    //geo.setinitdecoy( mcdecstr, numseq );
    for ( int i=1; i<numseq; i++ )
    {
        vars[((i-1)*2)]=mcdecstr[i].psi;
        vars[((i-1)*2)+1]=mcdecstr[i].phi;
    }

    // Parameters for LBFGS
    int max_iterations = 2000;
    int scale_factor = 2.0;
    int problem_dim = 2 * ( numseq - 1 );
    int MAX_CYCLE = 10;
    string converge_test = "absolute";

    // Run LBFGS folding, probably only long proteins will use the full 10 cycles
    for ( int i=0; i<MAX_CYCLE; i++ ) 
    {
        Minimizer fold;
        cout << "Running LBFGS cycle " << i+1 << " of " << MAX_CYCLE << endl; 
        fold.run( vars, mcdecstr, problem_dim, max_iterations, converge_test );
        if( i==1 )
        {
            energyFunction.read_energy_weights( weight_file2.c_str() );
        }
        if( i==4 )
        {
            energyFunction.read_energy_weights( weight_file3.c_str() );
        }  
    }
    
    // Output final model
    {
        FILE *fp;
        string outfile_name = datadir+"/model.pdb";
        fp=fopen(outfile_name.c_str(),"wt");
        int indatom=1;
        for ( int i=0; i<numseq; i++ )
        {
            const char *resn;
            for ( int j=0; j<26;j++)
            {
                if (mcdecstr[i].aaa==aad1[j])
                {
                    resn=aad3[j];
                    break;
                }
            }

            fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                indatom++, " N  ", resn, i+1,
                mcdecstr[i].ptn.x,     mcdecstr[i].ptn.y,     mcdecstr[i].ptn.z);
            fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                indatom++, " CA ", resn, i+1,
                mcdecstr[i].x,     mcdecstr[i].y,     mcdecstr[i].z);
            if ( strcmp(resn,"GLY") !=0 )
            {
                fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
                    indatom++, " CB ", resn, i+1,
                    mcdecstr[i].ptb.x,     mcdecstr[i].ptb.y,     mcdecstr[i].ptb.z);
            }
            fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                indatom++, " C  ", resn, i+1,
                mcdecstr[i].ptc.x,     mcdecstr[i].ptc.y,     mcdecstr[i].ptc.z);
            fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",			 
                indatom++, " O  ", resn, i+1,
                mcdecstr[i].pto.x,     mcdecstr[i].pto.y,     mcdecstr[i].pto.z);

        }
        fclose(fp);
    }
  
    // Free memory
    delete[]mcdecstr;
    delete[]vars;
    mcdecstr=NULL;
    vars=NULL;

    return 0;
}


