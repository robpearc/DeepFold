/*******************************************************************************************
**  
**  Functions for calculating energy terms and gradients.
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef ENERGY_H
#define ENERGY_H

#include "Operations.h"
#include "CubicSpline.h"
#include "CommonParameters.h"
#include "Derivatives.h"

class Energy  
{
public:
    Energy();

    bool load_dfire( const char *filename );
    bool load_general_contact( const char *filename );
    void load_files( int longestH, int seqnum, string datadir, string libdir, string weight_file,
                         string restraint_path, string ca_dist_path, string cb_cont_path, string ca_cont_path );

    void read_ca_distance( string distance_path );
    void torsion_energy(point3f *decstr, double *enelist);
    void torsion_energy(point3f *decstr, double *enelist,double *dE_dvars );

    void read_pred_phi_psi(const char *phi_file,const char *psi_file);

    /* prepare pn compact data structure to host the xyz coordinate. initialize
     * Epair_mat, Epair_tmp, CA_SC_mat, CA_SC_tmp, CB_CT_mat, CB_CT_tmp */
    void genpn();

    void pairwise_orientation( point3f *decstr );

    //double energytorsion2(point3f *decstr,int numseq);//log
    void calcallenergy(point3f *decstr, double *enelist);

    /* Calculates energy and gradients when grad_flag is true. Gradients are stored in dE_dvars */
    double calcrmsdenergy( point3f *decstr, double *vars, double *dE_dvars, bool grad_flag );

    /* Only calculates energy without gradients */
    double calcrmsdenergy( point3f *decstr, double *vars );
    
    /* Calculates the backbone hydrogen bonding energy*/
    void energy_hbond( point3f *decstr, double *enelist );

    //void scalpredsolve(char *queseqdata);
    //bool loadsolseq(char *filename);

    bool read_energy_weights( const char *weightfile );
    void read_cb_cont( string cont_path );
    void read_ca_cont( string cont_path );
    void read_cb_distance( string distance_path );
    void read_omega( string omega_path );
    void read_phi( string phi_path );
    void read_theta( string theta_path );

    virtual ~Energy();

private:
    Derivatives deriv;
    Geometry geo;
  
    int SHORT_RANGE_MIN, SHORT_RANGE_MAX;
    int MED_RANGE_MIN, MED_RANGE_MAX;
    int LONG_RANGE_MIN;
    int numseq;
    bool flag_grad;
  
    int nbin = 37;
    int nbin_omega = 25;
    int nbin_theta = nbin_omega;
    int nbin_phi = 13;


    point3d **atom_derivatives_f1;
    point3d **atom_derivatives_f2;

    double **Pcontact;
    double **PcontactCA;
    double ***dist_energy;
    double  **cb_dist_prob;
    double   *cb_dist_bin_value; 
    double ***ca_dist_energy;
    double  **ca_dist_prob;
    double   *ca_dist_bin_value; 
    double ***omega_energy;
    double  **omega_prob;
    double   *omega_bin_value; 
    double ***theta_energy;
    double  **theta_prob;
    double   *theta_bin_value; 
    double ***phi_energy;
    double  **phi_prob;
    double   *phi_bin_value; 

    bool flagTheta, flagPhi, flagOmega, flagPdist;

    CubicSpline **spline_dist;
    CubicSpline **spline_dist_ca;
    CubicSpline **spline_omega;
    CubicSpline **spline_theta;
    CubicSpline **spline_phi;
    CubicSpline **spline_hb_aa;
    CubicSpline **spline_hb_bb;
    CubicSpline **spline_hb_cc;
    CubicSpline **spline_poladat;
    CubicSpline **spline_rw; 

    double *pred_phi,*pred_psi;
    bool flag_ca_dist,flag_ca_cont,flag_cb_cont;
    double dist_epsilon;

    /*
    int solnum;
    double *solseq;
    */

    double weights[MAX_ENERGY_TERM_NUM];
    double dist_cut_cb, dist_cut_ca, theta_cut, phi_cut, omg_cut, hb_cut;
    double cont_cb_cut, cont_ca_cut;

    // Coordinate arrays
    double **pn;
    double **Epair_mat;
    double **Epair_tmp;
    double **CA_SC_mat;
    double **CA_SC_tmp;
    double **CB_CT_mat;
    double **CB_CT_tmp;
    double **PHI_mat;
    double **PHI_tmp;
    double **THETA_mat;
    double **THETA_tmp;
    double **OMEGA_mat;
    double **OMEGA_tmp;

    double enelist[20];//used for output
    double ***kbp;
    double *poladat[20][120];
    double *poladat2[20][100];
    double longestdist;
    /////////////////////////restraints
     
    int protype;//-1  1 helix 2 beta
    int dptype;//-1 1 init.dat
    bool flagcont;//use contact restraints or not;
    double wtcont[5];
    double wtinitdp;//init.dat as dist prof
    double wdpE;//weight factor for distProf from fragment
    double d8,d10,da,db,dc,dd;
};

Energy::Energy()  // initialization of parameters/options
{
    SHORT_RANGE_MIN     =  1; 
    SHORT_RANGE_MAX     = 11;
    MED_RANGE_MIN       = 12; 
    MED_RANGE_MAX       = 23;
    LONG_RANGE_MIN      = 24;
    dist_epsilon = 1e-8;
    for ( int i=0; i<MAX_ENERGY_TERM_NUM; i++ )
    {
        weights[i] = 0.0;
    }   


    flag_grad=true;

    int i,j;
    for(i=0;i<20;i++) for(j=0;j<100;j++) poladat[i][j]=NULL;
    for(i=0;i<20;i++) for(j=0;j<100;j++) poladat2[i][j]=NULL;
    //solseq=NULL;

    pred_phi=NULL;
    pred_psi=NULL;
    spline_dist=NULL;
    spline_dist_ca=NULL;
    spline_omega=NULL;
    spline_theta=NULL;
    spline_phi=NULL;
    spline_hb_aa=NULL;
    spline_hb_bb=NULL;
    spline_hb_cc=NULL;
    spline_poladat=NULL;
    spline_rw=NULL;
    kbp=NULL;
   
    /* Chengxin's distance matrix */
    pn=NULL; // xyz
    Epair_mat=NULL;
    Epair_tmp=NULL;
    CA_SC_mat=NULL; // CA-CA distance and SC-SC distance
    CA_SC_tmp=NULL;
    CB_CT_mat=NULL; // CB-CB distance and center-center distance
    CB_CT_tmp=NULL;

    PHI_mat=NULL;   // (i,j) is <CAi-CBi-CBj>
    PHI_tmp=NULL;   // (i,j) is <CAi-CBi-CBj>
    THETA_mat=NULL; // (i,j) is <Ni-CAi-CBi-CBj>
    THETA_tmp=NULL; // (i,j) is <Ni-CAi-CBi-CBj>
    OMEGA_mat=NULL; // (i,j) and (j,i) are both <CAi-CBi-CBj-CAj>
    OMEGA_tmp=NULL; // (i,j) and (j,i) are both <CAi-CBi-CBj-CAj>
    /* chengxin distance matrix end */

    longestdist=0;
    protype=-1;
    dptype=-1;
    flagcont=false;
    wtcont[0]=0;wtcont[1]=0;wtcont[2]=0;wtcont[3]=0;wtcont[4]=0;
    wtinitdp=5.0;
}

Energy::~Energy()
{
}

void Energy::load_files( int longestH, int seqnum, string datadir, string libdir, string weight_file, 
                         string restraint_path, string ca_dist_path, string cb_cont_path, string ca_cont_path )
{
    read_energy_weights( weight_file.c_str() );
    
    // Side-chain center information
    string sgposfile = libdir+"/newsgdistriaanew72.txt";
    geo.loadsgpos2(sgposfile.c_str(), 72);

    // Radius of gyration parameter
    longestdist=2.30431+1.42339*longestH;
    numseq=seqnum;
    genpn();

    string phi_file=datadir+"/phi.txt";
    string psi_file=datadir+"/psi.txt";
    read_pred_phi_psi(phi_file.c_str(),psi_file.c_str());

    // Contact parameters
    double dwell=12;
    d8=8.0;
    if     (numseq<100){ dwell=6; }
    else if(numseq<120){ dwell=8; }
    else if(numseq<200){ dwell=10; }
    d10=d8+dwell;
    da=(d8+d10)/2; // middle of d8 and d10
    db=dwell;      // width of first well
    dc=(d10+80)/2; // middle of d10 and 80
    dd=80-d10;     // width of second well

    // Allocate memory for splines
    spline_dist = new2DArrT< CubicSpline >( numseq, numseq );
    spline_dist_ca = new2DArrT< CubicSpline >( numseq, numseq );
    spline_omega = new2DArrT< CubicSpline >( numseq, numseq );
    spline_phi = new2DArrT< CubicSpline >( numseq, numseq );
    spline_theta = new2DArrT< CubicSpline >( numseq, numseq );
    spline_hb_aa = new2DArrT< CubicSpline >( numseq, numseq );
    spline_hb_bb = new2DArrT< CubicSpline >( numseq, numseq );
    spline_hb_cc = new2DArrT< CubicSpline >( numseq, numseq );
    spline_poladat = new2DArrT< CubicSpline >( 20, 100 );
    spline_rw = new2DArrT< CubicSpline >( 160, 160 );

    // Allocate memory for atom derivatives
    atom_derivatives_f1 = new point3d*[numseq];
    atom_derivatives_f2 = new point3d*[numseq];
    for(int i=0;i<numseq;i++)
    {
        atom_derivatives_f1[i]=new point3d[8];
        atom_derivatives_f2[i]=new point3d[8];
    }

    // Deep learning restraints
    bool flagcon = checkfile(restraint_path);
    if( flagcon )
    {
        string distance_path = restraint_path + "_dist.txt";
        string theta_path = restraint_path + "_theta.txt";
        string omega_path = restraint_path + "_omega.txt";
        string phi_path = restraint_path + "_phi.txt";
        
        read_cb_distance( distance_path );
        read_omega( omega_path );
        read_theta( theta_path );
        read_phi( phi_path );
    }
    
    //----------------- CB Distance Energy ------------------->
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( cb_dist_bin_value, dist_energy[i][j], nbin-1 );
            spline_dist[i][j] = SPL;
        }
    }

    //----------------- CA Distance Energy ------------------->
    flag_ca_dist = checkfile(ca_dist_path);
    if ( flag_ca_dist )
    {
        read_ca_distance( ca_dist_path );
        for( int i=0; i<numseq; i++ ){
            for( int j=i+1; j<numseq; j++ ){
                CubicSpline SPL;
                SPL.set_points( ca_dist_bin_value, ca_dist_energy[i][j], nbin-1 );
                spline_dist_ca[i][j] = SPL;
            }
        }
    }
    else
    {
        cout << "Error: No CA distance file! Setting weight to 0, " 
             << "which may decrease the modeling performance." << endl;
        weights[3]=weights[4]=weights[5]=0.0;
    }

    //----------------- CB Contact Energy ------------------->    
    flag_cb_cont = checkfile(cb_cont_path);
    if ( flag_cb_cont )
    {
        read_cb_cont( cb_cont_path );
    }
    else
    {
        cout << "Error: No CB contact file! Setting weight to 0, " 
             << "which may decrease the modeling performance."<< endl;
        weights[6]=weights[7]=weights[8]=0.0;
    }

    //----------------- CA Contact Energy ------------------->    
    flag_ca_cont = checkfile(ca_cont_path);
    if ( flag_ca_cont )
    {
        read_ca_cont( ca_cont_path );
    }
    else
    {
        cout << "Error: No CA contact file! Setting weight to 0, "
             << "which may decrease the modeling performance."<< endl;
        weights[9]=weights[10]=weights[11]=0.0;
    }

    //----------------- OMEGA Energy ------------------->
    // Fit spline, note currently we do not institute periodic boundary conditions
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( omega_bin_value, omega_energy[i][j], nbin_omega-1 );
            spline_omega[i][j] = SPL;
        }
    }

    //----------------- PHI Energy ------------------->
    // Fit spline, note currently we do not institute periodic boundary conditions
    for( int i=0; i<numseq; i++ ){
        for( int j=0; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( phi_bin_value, phi_energy[i][j], nbin_phi-1+6 );
            spline_phi[i][j] = SPL;
        }
    }

    //----------------- THETA Energy ------------------->
    // Fit spline, note currently we do not institute periodic boundary conditions
    for( int i=0; i<numseq; i++ ){
        for( int j=0; j<numseq; j++ ){
            CubicSpline SPL;
            SPL.set_points( theta_bin_value, theta_energy[i][j], nbin_theta-1 );
            spline_theta[i][j] = SPL;
        }
    }

    //------------- General contact energy --------------->
    string filename=libdir+"/sgpolarity5.txt";
    load_general_contact(filename.c_str());
    int npoladat_bin=31;
    double *dist_bin_poladat = new double[npoladat_bin];//new1DArr( nkbp_bin );
    dist_bin_poladat[0]=0.25;
    for( int i=1; i<npoladat_bin; i++ )
    {
        dist_bin_poladat[i] = dist_bin_poladat[i-1]+0.5;
    }

    for(int i=0;i<20;i++)
    {
        for(int j=0;j<100;j++)
        {
            CubicSpline SPL;
            SPL.set_points( dist_bin_poladat, poladat2[i][j], npoladat_bin );
            spline_poladat[i][j] = SPL;
        }
    }
    delete[]dist_bin_poladat;
    dist_bin_poladat=NULL;

    //--------------- Load Dfire-like Potential ---------->
    filename=libdir+"/data.dat";
    load_dfire(filename.c_str());

    int nkbp_bin = 30;
    double *dist_bin_kbp = new double[nkbp_bin];
    dist_bin_kbp[0]=0.25;
    for( int i=1; i<nkbp_bin; i++ )
    {
        dist_bin_kbp[i] = dist_bin_kbp[i-1]+0.5;
    }

    // Fit spline
    for( int j=1; j<159; j++ )
    {
        for( int i=j; i<159; i++ )
	{
            CubicSpline SPL;
            SPL.set_points( dist_bin_kbp, kbp[i][j], nkbp_bin );
            spline_rw[i][j] = SPL;
            spline_rw[j][i] = SPL;
        }
    }
    delete[]dist_bin_kbp;
    dist_bin_kbp=NULL;

    //---- Load Predicted Solvation ----->
    /*
    filename="sol.txt";
    loadsolseq(filename);
    */
}

void Energy::read_pred_phi_psi( const char *phi_file, const char *psi_file )
{
    FILE *file;
    if ( ( file=fopen(phi_file,"rt") )!=NULL )
    {
        pred_phi=new double[numseq];
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,file))
        {
            int res;
            double val=0.0;
            sscanf(line,"%d %lf",&res,&val);
            if(val<0) val+=360.0;
            pred_phi[res-1]=val*raddeg;
        }
    }
    fclose(file);


    if ( ( file=fopen(psi_file,"rt") ) !=NULL )
    {
        pred_psi=new double[numseq];
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,file))
        {
            int res;
            double val=0.0;
            sscanf(line,"%d %lf",&res,&val);
            if(val<0) val+=360.0;
            pred_psi[res-1]=val*raddeg;
        }
    }
    fclose(file);
}

bool Energy::read_energy_weights( const char *weightfile )
{
    for(int i=0;i<MAX_ENERGY_TERM_NUM;i++)
    {
        weights[i]=0.0;
    }
    FILE* pf=fopen(weightfile,"r");
    if(pf!=NULL)
    {
        char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
        while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pf))
        {
            char term[MAX_LENGTH_ONE_LINE_IN_FILE+1];
            double val=0.0;
            sscanf(line,"%s %lf",term,&val);
                 if(!strcmp(term,"Deep_Short_Range_Dist_CB"))  weights[ 0]=val;
            else if(!strcmp(term,"Deep_Med_Range_Dist_CB"))    weights[ 1]=val;
            else if(!strcmp(term,"Deep_Long_Range_Dist_CB"))   weights[ 2]=val;
            else if(!strcmp(term,"Deep_Short_Range_Dist_CA"))  weights[ 3]=val;
            else if(!strcmp(term,"Deep_Med_Range_Dist_CA"))    weights[ 4]=val;
            else if(!strcmp(term,"Deep_Long_Range_Dist_CA"))   weights[ 5]=val;
            else if(!strcmp(term,"Deep_Short_Range_Cont_CB"))  weights[ 6]=val;
            else if(!strcmp(term,"Deep_Med_Range_Cont_CB"))    weights[ 7]=val;
            else if(!strcmp(term,"Deep_Long_Range_Cont_CB"))   weights[ 8]=val;
            else if(!strcmp(term,"Deep_Short_Range_Cont_CA"))  weights[ 9]=val;
            else if(!strcmp(term,"Deep_Med_Range_Cont_CA"))    weights[10]=val;
            else if(!strcmp(term,"Deep_Long_Range_Cont_CA"))   weights[11]=val;
            else if(!strcmp(term,"Deep_Short_Range_Omg"))      weights[12]=val;
            else if(!strcmp(term,"Deep_Med_Range_Omg"))        weights[13]=val;
            else if(!strcmp(term,"Deep_Long_Range_Omg"))       weights[14]=val;
            else if(!strcmp(term,"Deep_Short_Range_Theta"))    weights[15]=val;
            else if(!strcmp(term,"Deep_Med_Range_Theta"))      weights[16]=val;
            else if(!strcmp(term,"Deep_Long_Range_Theta"))     weights[17]=val;
            else if(!strcmp(term,"Deep_Short_Range_Phi"))      weights[18]=val;
            else if(!strcmp(term,"Deep_Med_Range_Phi"))        weights[19]=val;
            else if(!strcmp(term,"Deep_Long_Range_Phi"))       weights[20]=val;
            else if(!strcmp(term,"Deep_HB_AA"))                weights[21]=val;
            else if(!strcmp(term,"Deep_HB_BB"))                weights[22]=val;
            else if(!strcmp(term,"Deep_HB_CC"))                weights[23]=val;
            else if(!strcmp(term,"General_Energy_VDW"))        weights[24]=val;
            else if(!strcmp(term,"General_Energy_HB"))         weights[25]=val;
            else if(!strcmp(term,"General_Energy_Dfire"))      weights[26]=val;
            else if(!strcmp(term,"General_Energy_SG_Contact")) weights[27]=val;
            else if(!strcmp(term,"General_Energy_Rad_Gyr"))    weights[28]=val;
            else if(!strcmp(term,"Tor_Constr"))                weights[29]=val;
            else if(!strcmp(term,"PCUT_DIST_CB"))              dist_cut_cb=val;
            else if(!strcmp(term,"PCUT_DIST_CA"))              dist_cut_ca=val;
            else if(!strcmp(term,"PCUT_CONT_CB"))              cont_cb_cut=val;
            else if(!strcmp(term,"PCUT_CONT_CA"))              cont_ca_cut=val;
            else if(!strcmp(term,"PCUT_OMG"))                      omg_cut=val;
            else if(!strcmp(term,"PCUT_PHI"))                      phi_cut=val;
            else if(!strcmp(term,"PCUT_THETA"))                  theta_cut=val;
            else if(!strcmp(term,"PCUT_HB"))                        hb_cut=val;
        }
        fclose(pf);
    }
    else{
	cout << "The energy weights file was not found!" << endl;
	return false;
    }

    for(int i=0;i<MAX_ENERGY_TERM_NUM;i++)
    {
        weights[i]*=1.0;
    }

    return true;
}

void Energy::read_cb_distance( string distance_path )
{

    double ***Pdist_cb = new3DArr( numseq, numseq, nbin );
    dist_energy = new3DArr( numseq, numseq, nbin-1 );
    cb_dist_bin_value = new1DArr( nbin - 1 );
    cb_dist_prob = new2DArr( numseq, numseq );

    bool flag = checkfile( distance_path );
    if( !flag )
    {
        //transform_npz2txt( );
        flag = checkfile( distance_path );
        if( !flag )
        {
            cout<<"There is no predicted distance file: "<<distance_path<<endl;
            return;
        }
    }
        
    ifstream in ( distance_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin + 3; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            double prob = atof( content[i].c_str() );
            Pdist_cb[a-1][b-1][i-2] = prob;
            if(i==2){
                double prob2 = atof( content[i].c_str() );
                cb_dist_prob[a-1][b-1]=prob;
            }
        }
    }

    double dist_bin[37][2];
    dist_bin[0][0] = 20.0, dist_bin[0][1] = 999.0;
    dist_bin[1][0] = 2.5, dist_bin[1][1] = 2.5;
    cb_dist_bin_value[0] = -0.0001;
    cb_dist_bin_value[1] = 2.75;
    for(int i=3; i<nbin; i++ ){
        cb_dist_bin_value[i-1] = cb_dist_bin_value[i-2]+0.5;//0.5*( dist_bin[i][0] + dist_bin[i][1] );
    }

    // Calculate Distance Energy
    for( int i=0; i<numseq-1; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            double Pref = 1E-8;
            for( int k=1; k<nbin; k++ )
            {
                double score = -log( ( Pdist_cb[i][j][k] + Pref ) / ( Pref + Pdist_cb[i][j][nbin-1] ) );
                dist_energy[i][j][k-1] = score;
                dist_energy[j][i][k-1] = score;
            }
            double score0 = max( 10.0, dist_energy[i][j][1] + 4.0 );
            dist_energy[i][j][0] = score0;
            dist_energy[j][i][0] = score0;
       }
    }

    // get contact map from distance
/*    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            for ( int b=0; b<12; b++ ){
                Pcontact[i][j] += Pdist[i][j][b];
            }
        }
    }
*/
    //release3DArr( numseq, numseq, Pdist_cb );
    cout<<"The predicted distance file is available."<<endl;
    flagPdist = true;
    in.close();
    vector<string >().swap( content );
}

void Energy::read_cb_cont( string cont_path )
{
    Pcontact = new2DArr( numseq, numseq );
    
    ifstream in ( cont_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<3; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            double prob = atof( content[i].c_str() );
            Pcontact[a-1][b-1] = prob;
        }
    }
}

void Energy::read_ca_cont( string cont_path )
{
    PcontactCA = new2DArr( numseq, numseq );
    
    ifstream in ( cont_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<3; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            double prob = atof( content[i].c_str() );
            PcontactCA[a-1][b-1] = prob;
        }
    }
}

void Energy::read_ca_distance( string distance_path )
{
    double ***Pdist = new3DArr( numseq, numseq, nbin );
    ca_dist_energy = new3DArr( numseq, numseq, nbin-1 );
    //Pcontact_ca = new2DArr( numseq, numseq );
    ca_dist_bin_value = new1DArr( nbin - 1 );
    ca_dist_prob = new2DArr( numseq, numseq );

    bool flag = checkfile( distance_path );
    if( !flag )
    {
        //transform_npz2txt( );
        flag = checkfile( distance_path );
        if( !flag )
	{
            cout<<"There is no predicted distance file: "<<distance_path<<endl;
            return;
        }
    }

    ifstream in ( distance_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin + 3; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=3; i<content.size(); i++){
            double prob = atof( content[i].c_str() );
            Pdist[a-1][b-1][i-3] = prob;
            if(i==content.size()-1){
                double prob2 = atof( content[i].c_str() );
                ca_dist_prob[a-1][b-1]=prob;
            }
        }
    }

    double dist_bin[37][2];
    dist_bin[0][0] = 20.0, dist_bin[0][1] = 999.0;
    dist_bin[1][0] = 2.5, dist_bin[1][1] = 2.5;
    ca_dist_bin_value[0] = -0.0001;
    ca_dist_bin_value[1] = 2.75;
    for(int i=3; i<nbin; i++ ){
        ca_dist_bin_value[i-1] = ca_dist_bin_value[i-2]+0.5;//0.5*( dist_bin[i][0] + dist_bin[i][1] );
    }

    // Calculate Distance Energy
    for( int i=0; i<numseq-1; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            double Pref = 1E-8;
            for( int k=0; k<nbin-1; k++ )
            {
                double score = -log( ( Pdist[i][j][k] + Pref ) / ( Pref + Pdist[i][j][nbin-2] ) );
                ca_dist_energy[i][j][k] = score;
                ca_dist_energy[j][i][k] = score;
            }
            double score0 = max( 10.0, ca_dist_energy[i][j][1] + 4.0 );
            ca_dist_energy[i][j][0] = score0;
            ca_dist_energy[j][i][0] = score0;
       }
    }

    // get contact map from distance
/*    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            for ( int b=0; b<12; b++ ){
                Pcontact[i][j] += Pdist[i][j][b];
            }
        }
    }
*/
    //release3DArr( numseq, numseq, Pdist );
    cout<<"The predicted CA distance file is available."<<endl;
    flagPdist = true;
    in.close();
    vector<string >().swap( content );
}

void Energy::read_omega( string omega_path )
{
    double ***Pomega = new3DArr( numseq, numseq, nbin_omega );
    omega_energy = new3DArr( numseq, numseq, nbin_omega );
    omega_bin_value = new1DArr( nbin_omega-1 );
    omega_prob = new2DArr( numseq, numseq );

    bool flag = checkfile( omega_path );
    if( !flag ){
        //transform_npz2txt( );
        flag = checkfile( omega_path );
        if( !flag ){
            cout<<"There is no predicted omega angle: "<<omega_path<<endl;
            return;
        }
    }

    ifstream in ( omega_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin_omega + 2; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            Pomega[a-1][b-1][i-2] = atof( content[i].c_str() );
	    //cout << "Pomega " <<a <<" " <<b<<" "<<i-2 <<" " <<Pomega[a-1][b-1][i-2] <<endl;
            if(i==2) omega_prob[a-1][b-1] = atof( content[i].c_str() );
        }
    }

    double omega_bin[25][2];
    omega_bin[0][0]   = 180.0*raddeg, omega_bin[0][1] = 999.0;
    omega_bin[1][0]   = -180.0*raddeg, omega_bin[1][1] = -165.0*raddeg;
    omega_bin_value[0]= -180.0*raddeg+7.5*raddeg;
    for(int i=2; i<nbin_omega; i++ ){
        omega_bin[i][0] = omega_bin[i-1][0] + 15.0*raddeg;
        omega_bin[i][1] = omega_bin[i-1][1] + 15.0*raddeg;
        omega_bin_value[i-1] = 0.5*( omega_bin[i][0] + omega_bin[i][1] );
    }

    //calculate omega energy
    for( int i=0; i<numseq; i++ ){
        for( int j=i+1; j<numseq; j++ ){
            for( int k=0; k<nbin_omega; k++ ){
                double prob = Pomega[i][j][k+1];
                if( k == nbin_omega-1 ) prob = Pomega[i][j][1];
                double score = -1.0*log( ( prob + dist_epsilon ) );
                omega_energy[i][j][k] = score;
                omega_energy[j][i][k] = score;
            }
        }
    }

    release3DArr( numseq, numseq, Pomega );
    cout<<"The predicted omega angle file is available."<<endl;
    flagOmega = true;
    in.close();
    vector<string >().swap( content );
}

void Energy::read_theta( string Theta_path )
{
    double ***Ptheta = new3DArr( numseq, numseq, nbin_theta );
    theta_energy = new3DArr( numseq, numseq, nbin_theta );
    theta_prob = new2DArr( numseq, numseq );
    theta_bin_value = new1DArr( nbin_theta-1 );

    bool flag = checkfile( Theta_path );
    if( !flag )
    {
        flag = checkfile( Theta_path );
        if( !flag )
        {
            cout << "There is no predicted theta angle: " << Theta_path << endl;
            return;
        }
    }

    ifstream in ( Theta_path.c_str() );
    string line;
    vector<string > content;
    for(int i=0; i<nbin_theta + 2; i++){
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) ){
        stringstream word( line );
        for(int i=0; i<content.size(); i++){
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++){
            Ptheta[a-1][b-1][i-2] = atof( content[i].c_str() );
            if(i==2) theta_prob[a-1][b-1] = atof( content[i].c_str() );
        }
    }

    double theta_bin[25][2];
    theta_bin[0][0]    =  180.0*raddeg, theta_bin[0][1] = 999.0;
    theta_bin[1][0]    = -180.0*raddeg, theta_bin[1][1] = -165.0*raddeg;
    theta_bin_value[0] = -180.0*raddeg+7.5*raddeg;
    for(int i=2; i<nbin_theta; i++ ){
        theta_bin[i][0] = theta_bin[i-1][0] + 15.0*raddeg;
        theta_bin[i][1] = theta_bin[i-1][1] + 15.0*raddeg;
        theta_bin_value[i-1] = 0.5*( theta_bin[i][0] + theta_bin[i][1] );
    }

    //calculate theta energy
    for( int i=0; i<numseq; i++ ){
        for( int j=0; j<numseq; j++ ){
            for( int k=0; k<nbin_theta; k++ ){
                double prob = Ptheta[i][j][k+1];
                if( k == nbin_theta-1 ) prob = Ptheta[i][j][1];
                double score = -1.0*log( ( prob + dist_epsilon ) );
                theta_energy[i][j][k] = score;
            }
        }
    }

    cout<<"The predicted theta angle file is available."<<endl;
    flagTheta = true;
    in.close();
    vector<string >().swap( content );
}


void Energy::read_phi( string Phi_path )
{
    double ***Pphi = new3DArr( numseq, numseq, nbin_phi );
    phi_energy = new3DArr( numseq, numseq, nbin_phi -1 +6);
    phi_bin_value = new1DArr( nbin_phi-1+6 );
    phi_prob = new2DArr( numseq, numseq );

    bool flag = checkfile( Phi_path );
    if ( !flag ){
        ///transform_npz2txt( );
        flag = checkfile( Phi_path );
        if( !flag )
        {
            cout<<"There is no predicted phi angle: "<<Phi_path<<endl;
                        return;
        }
    }

    ifstream in ( Phi_path.c_str() );
    string line;
    vector<string > content;
    for ( int i=0; i<nbin_phi + 2; i++ )
    {
        string x = "";
        content.push_back(x);
    }
    while ( getline (in, line) )
    {
        stringstream word( line );
        for( int i=0; i<content.size(); i++ )
        {
            word >> content[i];
        }
        int a = atoi(content[0].c_str());
        int b = atoi(content[1].c_str());
        if( a > numseq || b > numseq ) continue;
        for(int i=2; i<content.size(); i++)
        {
            Pphi[a-1][b-1][i-2] = atof( content[i].c_str() );
            if(i==2)
            {
                phi_prob[a-1][b-1] = atof( content[i].c_str() );
            }
        }
    }

    double phi_bin[13][2];
    phi_bin[0][0] = PI, phi_bin[0][1] = 999.0;
    phi_bin[1][0] = 0.0, phi_bin[1][1] = 15.0*raddeg;
    phi_bin_value[0]=7.5*raddeg-(15.0*3)*raddeg;
    phi_bin_value[1]=7.5*raddeg-(15.0*2)*raddeg;
    phi_bin_value[2]=7.5*raddeg-(15.0)*raddeg;
    phi_bin_value[3]=7.5*raddeg;

    for(int i=2; i<nbin_phi; i++ )
    {
        phi_bin[i][0] = phi_bin[i-1][0] + (15.0)*raddeg;
        phi_bin[i][1] = phi_bin[i-1][1] + (15.0)*raddeg;
        phi_bin_value[i-1+3] = 0.5*( phi_bin[i][0] + phi_bin[i][1] );
    }

    phi_bin_value[nbin_phi-1+3]=phi_bin_value[nbin_phi-2+3]+15.0*raddeg;
    phi_bin_value[nbin_phi-1+4]=phi_bin_value[nbin_phi-2+3]+15.0*2*raddeg;
    phi_bin_value[nbin_phi-1+5]=phi_bin_value[nbin_phi-2+3]+15.0*3*raddeg;
        
    //calculate phi energy
    for( int i=0; i<numseq; i++ )
    {
        for( int j=0; j<numseq; j++ )
        {
            if( i==j ) continue;
            for( int k=0; k<nbin_phi-1; k++ )
            {
                double prob = Pphi[i][j][k+1];
                double score = -1.0*log( ( prob + dist_epsilon ) );
                phi_energy[i][j][k+3] = score;
            }
            phi_energy[i][j][0]              = phi_energy[i][j][2+3];
            phi_energy[i][j][1]              = phi_energy[i][j][1+3];
            phi_energy[i][j][2]              = phi_energy[i][j][0+3];
            phi_energy[i][j][nbin_phi-1+3]   = phi_energy[i][j][nbin_phi-1-1+3];
            phi_energy[i][j][nbin_phi-1+1+3] = phi_energy[i][j][nbin_phi-1-2+3];
            phi_energy[i][j][nbin_phi-1+2+3] = phi_energy[i][j][nbin_phi-1-3+3];
        }
    }

    cout<<"The predicted phi angle file is available."<<endl;
    flagPhi = true;
    in.close();
    vector<string >().swap( content );
}

bool Energy::load_dfire( const char *filename )
{
    FILE *file;
    if( ( file=fopen(filename,"rt") ) ==NULL )
    {
	cout << "Error loading dfire file " << filename << " from library!" << endl;
        return false;
    }

    int *** countptr;
    countptr = new3DIntArr( 160, 160, 30 );
    kbp = new3DArr(160, 160, 29);

    char bf[2000];
    int  i,j, k, num;
    for( k=0; k<30; k++ )
    {
        fgets(bf,2000,file);
        for( j=1; j<159; j++ )
        {
            fgets(bf,2000,file);
            for( i=1; i<159; i++ )
            {
                sscanf(bf+(i-1)*10,"%10d",&num);
                countptr[i][j][k]=num;
            }
        }
        fgets(bf,2000,file);
    }
    
    double RR = 1.9858775;
    double TT = 300.00;
    double ALPHA = 1.61; 
    for( k=0; k<30; k++ )
    {
        for( j=1; j<159; j++ )
        {
            for( i=j; i<159; i++ )
            {
                if(countptr[i][j][k]==0 && countptr[j][i][k]==0)
                    kbp[i][j][k]=0.1;//origin james is 0.1, the same in dfire
                else if(i==j)
                    kbp[i][j][k] = ( 0.016784 ) * ( -0.001 ) * RR * TT * log( ( countptr[i][j][k] )
                                 /  pow ( ( ( k + 0.5 ) / 29.5 ), ALPHA ) / ( countptr[i][j][29] ) );
                else
                    kbp[i][j][k] = ( 0.016784 ) * ( -0.001 ) * RR * TT 
                                 * log( ( countptr[i][j][k] + countptr[j][i][k] ) 
                                 / pow( ( ( k + 0.5 ) / 29.5 ) , ALPHA ) 
                                 / ( countptr[i][j][29] + countptr[j][i][29] ) );
                kbp[j][i][k]=kbp[i][j][k];
            }
        }
    }
    fclose(file);
    release3DArr( 160, 160, countptr );
    return true;
}

bool Energy::load_general_contact( const char *filename )
{
    FILE *file;
    if ( ( file=fopen(filename,"rt") ) ==NULL )
    {
        printf("Unable to open polaritydata %s\n",filename);
        return false;
    }

    int i,j,k;
    char oneline[200];
    for(i=0;i<20;i++)
    {
        for(j=0;j<100;j++)
        {
            poladat[i][j]=new double[61];
            poladat2[i][j]=new double[31];
        }
    }
    for(i=0;i<20;i++)
    {
        for(j=0;j<100;j++)
        {
            for(k=0;k<61;k++)
            {
                fgets(oneline,200,file);
                sscanf(oneline,"%lf",&poladat[i][j][k]);
                poladat[i][j][k]+=3.65;
            }
            for(k=0;k<31;k++) {
                    poladat[i][j][k]-=poladat[i][j][30];
                    poladat2[i][j][k]=poladat[i][j][k];
            }
        }
        //printf("%d %d %d %f\n",i,19,40,poladat[i][19][9]);
    }
    fclose(file);
    return true;
}
/*
bool Energy::loadsolseq( char *filename )
{
    FILE *file;
    if((file=fopen(filename,"rt"))==NULL)
    {
        printf("Unable to open solseq %s\n",filename);
        return false;
    }
    int i,j;
    double tval;
    char oneline[200];
    fgets(oneline,200,file);
    sscanf(oneline,"%d",&solnum);
    if(solseq) delete[]solseq;
    solseq=new double[solnum];
    for(i=0;i<solnum;i++)
    {
        fgets(oneline,200,file);
        sscanf(oneline,"%d %lf %lf",&j,&tval,&solseq[i]);
        solseq[i]+=1.0;
        solseq[i]/=2.0;
    }
    fclose(file);
    return true;
}

void Energy::scalpredsolve(char *queseqdata)
{
    if(!solseq) return;
    int i,indaa;
    
    for(i=0;i<solnum;i++)
    {
        indaa=aminoid(queseqdata[i]);
        if(indaa>=20) indaa=5;
        solseq[i]*=scaa[indaa]*0.91;
        if(solseq[i]>1.0) solseq[i]=1.0;
    }
}
*/

/* Chengxin's structure arrays */
void Energy::genpn()
{
    int seqnum=numseq;
    pn=new double*[8*seqnum];
    int i,j,ii;
    for (i=0;i<seqnum;i++)
    {
        for (ii=0;ii<8;ii++)
        {
            pn[i*8+ii]=new double[3];
            pn[i*8+ii][0]=0;
            pn[i*8+ii][1]=0;
            pn[i*8+ii][2]=0;
        }
    }

    Epair_mat=new double*[seqnum];
    Epair_tmp=new double*[seqnum];
    CA_SC_mat=new double*[seqnum];
    CA_SC_tmp=new double*[seqnum];
    CB_CT_mat=new double*[seqnum];
    CB_CT_tmp=new double*[seqnum];

    PHI_mat  =new double*[seqnum];
    PHI_tmp  =new double*[seqnum];
    THETA_mat=new double*[seqnum];
    THETA_tmp=new double*[seqnum];
    OMEGA_mat=new double*[seqnum];
    OMEGA_tmp=new double*[seqnum];
    for (i=0;i<seqnum;i++)
    {
        Epair_mat[i]=new double[seqnum];
        Epair_tmp[i]=new double[seqnum];
        CA_SC_mat[i]=new double[seqnum];
        CA_SC_tmp[i]=new double[seqnum];
        CB_CT_mat[i]=new double[seqnum];
        CB_CT_tmp[i]=new double[seqnum];

        PHI_mat[i]  =new double[seqnum];
        PHI_tmp[i]  =new double[seqnum];
        THETA_mat[i]=new double[seqnum];
        THETA_tmp[i]=new double[seqnum];
        OMEGA_mat[i]=new double[seqnum];
        OMEGA_tmp[i]=new double[seqnum];
        for (j=0;j<seqnum;j++)
        {
            Epair_mat[i][j]=0;
            Epair_tmp[i][j]=0;
            CA_SC_mat[i][j]=0;
            CA_SC_tmp[i][j]=0;
            CB_CT_mat[i][j]=0;
            CB_CT_tmp[i][j]=0;
            PHI_mat[i][j]  =360;
            PHI_tmp[i][j]  =360;
            THETA_mat[i][j]=360;
            THETA_tmp[i][j]=360;
            OMEGA_mat[i][j]=360;
            OMEGA_tmp[i][j]=360;
        }
    }
}

/* calculate orientations of decoy structure */
void Energy::pairwise_orientation( point3f *decstr )
{
    int i,j,r1,r2;
    double THETA,OMEGA;
    point3d pd,pd2;
    int idx;
    for ( i=0; i<numseq; i++ )
    {
        r1=i*8;
        for ( j=i+1; j<numseq; j++ )
        {
            if ( decstr[i].aaa=='G' || decstr[j].aaa=='G' ) continue;
            r2=j*8;
            idx=(j-1)+(2*numseq-3-i)*i/2;

            if ( theta_prob[i][j]<theta_cut )
            {
                // THETAij = <Ni-CAi-CBi-CBj>  , [-pi,pi]
                THETA=phi( pn[r1+1][0], pn[r1+1][1], pn[r1+1][2],
                           pn[r1  ][0], pn[r1  ][1], pn[r1  ][2],
                           pn[r1+4][0], pn[r1+4][1], pn[r1+4][2],
                           pn[r2+4][0], pn[r2+4][1], pn[r2+4][2] );
                if ( THETA>180 ) THETA-=360;
                THETA_tmp[i][j]=THETA;
            }

            if ( theta_prob[j][i]<theta_cut )
            {
                // THETAji = <Nj-CAj-CBj-CBi>  , [-pi,pi]
                THETA=phi( pn[r2+1][0], pn[r2+1][1], pn[r2+1][2],
                           pn[r2  ][0], pn[r2  ][1], pn[r2  ][2],
                           pn[r2+4][0], pn[r2+4][1], pn[r2+4][2],
                           pn[r1+4][0], pn[r1+4][1], pn[r1+4][2] );
                if ( THETA>180 ) THETA-=360;
                THETA_tmp[j][i]=THETA;
            }

            if ( phi_prob[i][j]<phi_cut )
            {
                // PHIij   = <CAi-CBi-CBj>     , [0,pi]
                pd.x  = pn[r1  ][0]-pn[r1+4][0];
                pd.y  = pn[r1  ][1]-pn[r1+4][1];
                pd.z  = pn[r1  ][2]-pn[r1+4][2];
                pd2.x = pn[r2+4][0]-pn[r1+4][0];
                pd2.y = pn[r2+4][1]-pn[r1+4][1];
                pd2.z = pn[r2+4][2]-pn[r1+4][2];
                PHI_tmp[i][j] = angv( pd, pd2 );
            }

            if ( phi_prob[j][i]<phi_cut )
            {
                // PHIji   = <CAj-CBj-CBi>     , [0,pi]
                pd.x  = pn[r2  ][0]-pn[r2+4][0];
                pd.y  = pn[r2  ][1]-pn[r2+4][1];
                pd.z  = pn[r2  ][2]-pn[r2+4][2];
                pd2.x = pn[r1+4][0]-pn[r2+4][0];
                pd2.y = pn[r1+4][1]-pn[r2+4][1];
                pd2.z = pn[r1+4][2]-pn[r2+4][2];
                PHI_tmp[j][i] = angv( pd, pd2 );
            }

            if ( omega_prob[i][j]<omg_cut )
            {
                // OMEGAij = <CAi-CBi-CBj-CAj> , [-pi,pi]
                OMEGA=phi( pn[r1  ][0], pn[r1  ][1], pn[r1  ][2],
                           pn[r1+4][0], pn[r1+4][1], pn[r1+4][2],
                           pn[r2+4][0], pn[r2+4][1], pn[r2+4][2],
                           pn[r2  ][0], pn[r2  ][1], pn[r2  ][2] );
                if ( OMEGA>180 ) OMEGA-=360;
                OMEGA_tmp[j][i] = OMEGA_tmp[i][j] = OMEGA;
            }
        }
    }
    return;
}


/* calculate most of the pairwise energy terms in QUARK 
 * enelist   - energy terms */
void Energy::calcallenergy( point3f *decstr, double *enelist )
{
    int i,j,ii,jj; // i,j are residue indexes; ii,jj are atom indexes
    int tbin,ind1,ind2,indm,ind5[5],ind6[5];
    int atomid[]={1,17,0,26,2};
    int sgmap[]={2,1,3,4,99,99,99,0};
    //point3d pin[8],pjn[8];
    point3d tp;
    point3d atom1,atom2,atom3,atom4;
    double tdist,tdist2,tdist3;
    double tdist_SC,tdist2_SC;
    double tdist_CA,tdist2_CA;
    double tdist_CB,tdist2_CB;
    double distcut=9.2;
    double sqdistcut=9.2*9.2;
    double scalf=0.00707;
    double cutca=4.0;
    double sqdist,tdists;
  
    double weight=0;
    double mu;
    double wtdp=0.60;
    double wtsg=0.10;
    if(protype==2) //non-alpha proteins
    {
        wtdp=3.00;
        wtsg=0.08;
    }
    else if(protype==1) wtdp=0.80;//alpha proteins
    wtdp*=wdpE; //reweight the weight of distProf for QE/QN

    int sdist=17;
    //double *soldat= new double[numseq];
    //for(i=0;i<numseq;i++) soldat[i]=0;
    point3d tcen;
    tcen.x=0;tcen.y=0;tcen.z=0;
    enelist[0]=0;//I-TASSER restraints
    enelist[1]=0;//Excluded volume, ww1
    enelist[4]=0;//E_RWplus_backbone_only, RW part tuned by wRW
    enelist[5]=0;//pairwise side-chain atomic potential, 'sgpolarity5.txt', ww5
    enelist[6]=0;//sequence-specific solvation potential, 'sol.txt', ww6
    enelist[7]=0;//radias gyration, RG part by ww7
    enelist[8]=0;//distant profile from either fragments (QE/QN) or init.dat (QP/QT), wtdp
    enelist[12]=0;//Edist ResTriplet2 distance map
    enelist[14]=0;//CA bond-length, ww14
    enelist[15]=0;//Econtact
    enelist[16]=0;//Econtact
    enelist[17]=0;//Econtact
    double tcont[5]={0,0,0,0,0};//tasser contact

    int r1,r2;
    int idx; // for accessing distance map
    for ( i=0; i<numseq; i++ )
    {
        r1=i*8;

        pn[r1  ][0]=decstr[i].x;     //CA
        pn[r1  ][1]=decstr[i].y;
        pn[r1  ][2]=decstr[i].z;
        pn[r1+1][0]=decstr[i].ptn.x; //N
        pn[r1+1][1]=decstr[i].ptn.y;
        pn[r1+1][2]=decstr[i].ptn.z;
        pn[r1+2][0]=decstr[i].ptc.x; //C
        pn[r1+2][1]=decstr[i].ptc.y;
        pn[r1+2][2]=decstr[i].ptc.z;
        pn[r1+3][0]=decstr[i].pto.x; //O
        pn[r1+3][1]=decstr[i].pto.y;
        pn[r1+3][2]=decstr[i].pto.z;
        pn[r1+4][0]=decstr[i].ptb.x; //CB
        pn[r1+4][1]=decstr[i].ptb.y;
        pn[r1+4][2]=decstr[i].ptb.z;
        pn[r1+5][0]=decstr[i].pth.x; //H, used in hbond energy
        pn[r1+5][1]=decstr[i].pth.y;
        pn[r1+5][2]=decstr[i].pth.z;
        pn[r1+6][0]=decstr[i].ptg.x; //CT, residue center for solvation
        pn[r1+6][1]=decstr[i].ptg.y;
        pn[r1+6][2]=decstr[i].ptg.z;
        pn[r1+7][0]=decstr[i].ptsg.x;//SG
        pn[r1+7][1]=decstr[i].ptsg.y;
        pn[r1+7][2]=decstr[i].ptsg.z;
    }


    /* fill up CA_SC_tmp and CB_CT_tmp */
    for(i=0;i<numseq;i++)
    {
        r1=i*8;

        for(j=i+1;j<numseq;j++)
        {
            idx=(j-1)+(2*numseq-3-i)*i/2;
            r2=j*8;

            /* CA-CA distance calculation */
            tp.x=pn[r1][0]-pn[r2][0];
            tp.y=pn[r1][1]-pn[r2][1];
            tp.z=pn[r1][2]-pn[r2][2];
            tdist2_CA=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
            CA_SC_tmp[i][j]=tdist2_CA;

            /* SC-SC distance calculation */
            tp.x=pn[r1+7][0]-pn[r2+7][0];
            tp.y=pn[r1+7][1]-pn[r2+7][1];
            tp.z=pn[r1+7][2]-pn[r2+7][2];
            tdist2_SC=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
            CA_SC_tmp[j][i]=tdist2_SC;

            CB_CT_tmp[i][j]=CB_CT_tmp[j][i]=0;
            /* CB-CB distance calculation 
             * CA-CB bond length is in the range of [1.51371,1.54507]
             * CA-SG bond length is in the range of [0,4.75663]
             * 1.54507*2+15=18.09014 
             * (1.54507+4.75663)*2+15=27.6034     */
	    
            //if (tdist2_CA<18.1  || PredictCB[i][j]>0.0 || (dist_map_type &&
            //    ((ResTriplet2CB && ResTriplet2CB[idx])||
            //     (DMPfoldCB     && DMPfoldCB[idx]    )||
            //     (trRosettaCB   && trRosettaCB[idx]  )||
            //     (ResTriplet3CB && ResTriplet3CB[idx]))))
            //{
                tp.x=pn[r1+4][0]-pn[r2+4][0];
                tp.y=pn[r1+4][1]-pn[r2+4][1];
                tp.z=pn[r1+4][2]-pn[r2+4][2];
                CB_CT_tmp[i][j]=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
            //}

            /* CT-CT distance square. Note that we did not take sqrt */
            //if (tdist2_CA<=20.0 || tdist2_SC<=15.0)
            //{
                tp.x=pn[r1+6][0]-pn[r2+6][0];
                tp.y=pn[r1+6][1]-pn[r2+6][1];
                tp.z=pn[r1+6][2]-pn[r2+6][2];
                sqdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z;
                if(sqdist<1.0) sqdist=1.0;
                CB_CT_tmp[j][i]=sqdist;
            //}

        }
    }

    /* pairwise orientation terms
                               Cj
                              /
                     THETAji /
                   CBj----CAj             ---  covalent bond
                  .__/       \            ...  non-covalent interaction
                 .   PHIji    \
              D .              Nj
         OMEGA .                  
              .                   
             .                    
          CBi \                   Symmetric properties:         
           | _/ PHIij                 Dij     = |CBi-CBj|
    THETAij|                          OMEGAij = <CAi-CBi-CBj-CAj> , [-pi,pi]
           |                      Asymmetric properties:
          CAi                         THETAij = <Ni-CAi-CBi-CBj>  , [-pi,pi]
         /   \                        THETAji = <Nj-CAj-CBj-CBi>  , [-pi,pi]
        /     \                       PHIij   = <CAi-CBi-CBj>     , [0,pi]
      Ni       Ci                     PHIji   = <CAj-CBj-CBi>     , [0,pi] */
    {
        pairwise_orientation( decstr );
       
        for ( i=0; i<numseq; i++ )
        {
            int r1=i*8;
            for (j=i+1;j<numseq;j++)
            {
                int r2=j*8;
                if ( decstr[i].aaa=='G' || decstr[j].aaa=='G' ) continue;

                double weight=0.0;
                if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
                {
                    weight=weights[12];	   
                }
                else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
                {
                    weight=weights[13];	   
                }
                else if ( abs(i-j)>=LONG_RANGE_MIN )
                {
                    weight=weights[14];	   
                }
                else{
                    weight=0.0;
                }

                if (weight>0.0 && omega_prob[i][j]<omg_cut && OMEGA_tmp[i][j]!=180.0)
                {
                    enelist[16]+=weight*spline_omega[i][j].cubic_spline( ( OMEGA_tmp[i][j]*raddeg ) );
                    if ( flag_grad )
                    {
                        double dfunc = weight*spline_omega[i][j].cubic_dspline( ( OMEGA_tmp[i][j]*raddeg ) );
                        double phi=0.0;
                        atom1 = setv( pn[r1  ][0], pn[r1  ][1], pn[r1  ][2] );
                        atom2 = setv( pn[r1+4][0], pn[r1+4][1], pn[r1+4][2] );
                        atom3 = setv( pn[r2+4][0], pn[r2+4][1], pn[r2+4][2] );
                        atom4 = setv( pn[r2  ][0], pn[r2  ][1], pn[r2  ][2] );
                        
                        point3d scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        point3d f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom1_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[i][0] = addv( atom_derivatives_f1[i][0], scal_f1 );
                        atom_derivatives_f2[i][0] = addv( atom_derivatives_f2[i][0], scal_f2 ); 

                        scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom2_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[i][4] = addv( atom_derivatives_f1[i][4], scal_f1 );
                        atom_derivatives_f2[i][4] = addv( atom_derivatives_f2[i][4], scal_f2 ); 
                        
                        scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom2_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[j][4] = addv( atom_derivatives_f1[j][4], scal_f1 );
                        atom_derivatives_f2[j][4] = addv( atom_derivatives_f2[j][4], scal_f2 ); 

                        scal_f1 = setv( 0.0, 0.0, 0.0 ), scal_f2 = setv( 0.0, 0.0, 0.0 );
                        f1 = setv( 0.0, 0.0, 0.0 ), f2 = setv( 0.0, 0.0, 0.0 );
                        deriv.atom1_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1 = scal( f1, dfunc );
                        scal_f2 = scal( f2, dfunc );
                        atom_derivatives_f1[j][0] = addv( atom_derivatives_f1[j][0], scal_f1 );
                        atom_derivatives_f2[j][0] = addv( atom_derivatives_f2[j][0], scal_f2 ); 
                    } 
                }

                if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
                {
                    weight=weights[18];	   
                }
                else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
                {
                    weight=weights[19];	   
                }
                else if ( abs(i-j)>=LONG_RANGE_MIN )
                {
                    weight=weights[20];	   
                }
                else{
                    weight=0.0;
                }

                if ( weight>0.0 && phi_prob[i][j]<phi_cut )
                {
                    enelist[18]+=weight*spline_phi[i][j].cubic_spline( ( PHI_tmp[i][j] ) );

                    if ( flag_grad )
                    {
                        double dfunc=weight*spline_phi[i][j].cubic_dspline( ( PHI_tmp[i][j] ) ); 
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1_p1, f2_p1, f1_p2, f2_p2, f1_p3, f2_p3;
                        point3d scalf1_p1, scalf2_p1, scalf1_p2, scalf2_p2, scalf1_p3, scalf2_p3;
                        double theta=0.0;
                        atom1=setv(pn[r1  ][0],pn[r1  ][1],pn[r1  ][2]);  //Ca i
                        atom2=setv(pn[r1+4][0],pn[r1+4][1],pn[r1+4][2]);  //Cb i
                        atom3=setv(pn[r2+4][0],pn[r2+4][1],pn[r2+4][2]);  //Cb j

                        deriv.angular_deriv(atom1,atom2,atom3,f1_p1,f2_p1,f1_p2,f2_p2,f1_p3,f2_p3);
                        scalf1_p1=scal(f1_p1,dfunc);
                        scalf2_p1=scal(f2_p1,dfunc);
                        atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scalf1_p1);
                        atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scalf2_p1);
                        scalf1_p2=scal(f1_p2,dfunc);
                        scalf2_p2=scal(f2_p2,dfunc);
                        atom_derivatives_f1[i][4]=addv(atom_derivatives_f1[i][4],scalf1_p2);
                        atom_derivatives_f2[i][4]=addv(atom_derivatives_f2[i][4],scalf2_p2);
                        scalf1_p3=scal(f1_p3,dfunc);
                        scalf2_p3=scal(f2_p3,dfunc);
                        atom_derivatives_f1[j][4]=addv(atom_derivatives_f1[j][4],scalf1_p3);
                        atom_derivatives_f2[j][4]=addv(atom_derivatives_f2[j][4],scalf2_p3);
                    }
                }
                
                if ( weight>0.0 && phi_prob[j][i]<phi_cut )
                {
                    enelist[18] += weight*spline_phi[j][i].cubic_spline( ( PHI_tmp[j][i] ) );

                    if ( flag_grad )
                    {
                        double dfunc = weight * spline_phi[j][i].cubic_dspline( ( PHI_tmp[j][i] ) ); 
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1_p1, f2_p1, f1_p2, f2_p2, f1_p3, f2_p3;
                        point3d scalf1_p1, scalf2_p1, scalf1_p2, scalf2_p2, scalf1_p3, scalf2_p3;
                        double theta=0.0;
                        atom1 = setv( pn[r2  ][0], pn[r2  ][1], pn[r2  ][2] );  //Ca j
                        atom2 = setv( pn[r2+4][0], pn[r2+4][1], pn[r2+4][2] );  //Cb j
                        atom3 = setv( pn[r1+4][0], pn[r1+4][1], pn[r1+4][2] );  //Cb i
                        deriv.angular_deriv( atom1, atom2, atom3, f1_p1, f2_p1, 
                                             f1_p2, f2_p2, f1_p3, f2_p3 );
                        scalf1_p1 = scal( f1_p1, dfunc );
                        scalf2_p1 = scal( f2_p1, dfunc);
                        atom_derivatives_f1[j][0] = addv( atom_derivatives_f1[j][0], scalf1_p1 );
                        atom_derivatives_f2[j][0] = addv( atom_derivatives_f2[j][0], scalf2_p1 );
                        scalf1_p2 = scal( f1_p2, dfunc );
                        scalf2_p2 = scal( f2_p2, dfunc );
                        atom_derivatives_f1[j][4] = addv( atom_derivatives_f1[j][4], scalf1_p2 );
                        atom_derivatives_f2[j][4] = addv( atom_derivatives_f2[j][4], scalf2_p2 );
                        scalf1_p3 = scal( f1_p3, dfunc );
                        scalf2_p3 = scal( f2_p3, dfunc );
                        atom_derivatives_f1[i][4] = addv( atom_derivatives_f1[i][4], scalf1_p3 );
                        atom_derivatives_f2[i][4] = addv( atom_derivatives_f2[i][4], scalf2_p3 );              
                    }
                }

                if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
                {
                    weight=weights[15];	   
                }
                else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
                {
                    weight=weights[16];	   
                }
                else if ( abs(i-j)>=LONG_RANGE_MIN )
                {
                    weight=weights[17];	   
                }
                else{
                    weight=0.0;
                }

                if( weight>0.0 && theta_prob[i][j]<theta_cut && THETA_tmp[i][j]!=180.0 )
                {
                    enelist[17]+=weight*spline_theta[i][j].cubic_spline( ( THETA_tmp[i][j]*raddeg ) );

                    if(flag_grad)
                    {
                        double dfunc = weight*spline_theta[i][j].cubic_dspline( ( THETA_tmp[i][j]*raddeg ) );
                        double phi=0.0;
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1, f2;
                        atom1=setv(pn[r1+1][0],pn[r1+1][1],pn[r1+1][2]);  //N  i
                        atom2=setv(pn[r1  ][0],pn[r1  ][1],pn[r1  ][2]);  //Ca i
                        atom3=setv(pn[r1+4][0],pn[r1+4][1],pn[r1+4][2]);  //Cb i
                        atom4=setv(pn[r2+4][0],pn[r2+4][1],pn[r2+4][2]);  //Cb j

                        deriv.atom1_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][1]=addv(atom_derivatives_f1[i][1],scal_f1);
                        atom_derivatives_f2[i][1]=addv(atom_derivatives_f2[i][1],scal_f2); 
   
                        deriv.atom2_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scal_f1);
                        atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scal_f2); 

                        deriv.atom2_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][4]=addv(atom_derivatives_f1[i][4],scal_f1);
                        atom_derivatives_f2[i][4]=addv(atom_derivatives_f2[i][4],scal_f2); 

                        deriv.atom1_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][4]=addv(atom_derivatives_f1[j][4],scal_f1);
                        atom_derivatives_f2[j][4]=addv(atom_derivatives_f2[j][4],scal_f2);
                    } 
                }


                if( weight>0.0 && theta_prob[j][i]<theta_cut && THETA_tmp[j][i]!=180.0 )
                {
                    enelist[17]+=weight*spline_theta[j][i].cubic_spline( ( THETA_tmp[j][i]*raddeg ) );

                    if(flag_grad)
                    {
                        double dfunc = weight*spline_theta[j][i].cubic_dspline( ( THETA_tmp[j][i]*raddeg ) );
                        double phi=0.0;
                        point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                        point3d f1, f2;
                        atom1=setv(pn[r2+1][0],pn[r2+1][1],pn[r2+1][2]);  //N  j
                        atom2=setv(pn[r2  ][0],pn[r2  ][1],pn[r2  ][2]);  //Ca j
                        atom3=setv(pn[r2+4][0],pn[r2+4][1],pn[r2+4][2]);  //Cb j
                        atom4=setv(pn[r1+4][0],pn[r1+4][1],pn[r1+4][2]);  //Cb i

                        deriv.atom1_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][1]=addv(atom_derivatives_f1[j][1],scal_f1);
                        atom_derivatives_f2[j][1]=addv(atom_derivatives_f2[j][1],scal_f2); 
   
                        deriv.atom2_dih_deriv( atom1, atom2, atom3, atom4, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][0]=addv(atom_derivatives_f1[j][0],scal_f1);
                        atom_derivatives_f2[j][0]=addv(atom_derivatives_f2[j][0],scal_f2); 

                        deriv.atom2_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[j][4]=addv(atom_derivatives_f1[j][4],scal_f1);
                        atom_derivatives_f2[j][4]=addv(atom_derivatives_f2[j][4],scal_f2); 

                        deriv.atom1_dih_deriv( atom4, atom3, atom2, atom1, phi, f1, f2 );
                        scal_f1=scal(f1,dfunc);
                        scal_f2=scal(f2,dfunc);
                        atom_derivatives_f1[i][4]=addv(atom_derivatives_f1[i][4],scal_f1);
                        atom_derivatives_f2[i][4]=addv(atom_derivatives_f2[i][4],scal_f2); 
                    } 
                }
            }
        }
    }

    if ( weights[28]>0.0 )
    {
        for ( i=0; i<numseq; i++ )
        {
            r1=i*8;

            /* for radius of gyration */
            tcen.x+=pn[r1][0];
            tcen.y+=pn[r1][1];
            tcen.z+=pn[r1][2];
        }
    }
 
    /* calculate the rest of terms */
    //int CBpair=0;
    //int CApair=0;
    for(i=0;i<numseq;i++)
    {
        r1=i*8;
        ind1=decstr[i].iaa; // identity of 20 amino acids [1-20]
        ind5[0]=Tab5[ind1][atomid[0]];
        ind5[1]=Tab5[ind1][atomid[1]];
        ind5[2]=Tab5[ind1][atomid[2]];
        ind5[3]=Tab5[ind1][atomid[3]];
        ind5[4]=Tab5[ind1][atomid[4]];
        for(j=i+1;j<numseq;j++)
        {
            r2=j*8;
            indm=i*numseq+j;
            idx=(j-1)+(2*numseq-3-i)*i/2;

            ind2=decstr[j].iaa;
            ind6[0]=Tab5[ind2][atomid[0]];
            ind6[1]=Tab5[ind2][atomid[1]];
            ind6[2]=Tab5[ind2][atomid[2]];
            ind6[3]=Tab5[ind2][atomid[3]];
            ind6[4]=Tab5[ind2][atomid[4]];

            /**** CB-CB specific energy ****/
            tdist2_CB=CB_CT_tmp[i][j];
            tdist_CB=tdist2_CB*tdist2_CB;

            double weight_cont=0.0;
            if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
            {
                weight_cont=weights[6];
            }
            else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
            {
                weight_cont=weights[7];
            }
            else if ( abs(i-j)>=LONG_RANGE_MIN )
            {
                weight_cont=weights[8];
            }
            else
            {
                weight_cont=0.0;
            }

            if( Pcontact[i][j]>cont_cb_cut && weight_cont>0.0 )
            {
                double dfunc=0.0;
                double weight=Pcontact[i][j];    //weight for well depth
                if(tdist2_CB<=d8)
                {	//r<8A
                    enelist[15]+=-weight*weight_cont; // list data when weight=1
                }
                else if(tdist2_CB<d10) // 8<r<10
                {
                    enelist[15]+=-weight_cont*weight*(1-sin((tdist2_CB-da)/db*PI))/2;

                    if(flag_grad)
                    {
                        dfunc=weight_cont*weight*cos((tdist2_CB-da)/db*PI)/(2*db*PI);
                        point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                        point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                        point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                        atom1.x=decstr[i].ptb.x,atom1.y=decstr[i].ptb.y,atom1.z=decstr[i].ptb.z;
                        atom2.x=decstr[j].ptb.x,atom2.y=decstr[j].ptb.y,atom2.z=decstr[j].ptb.z;
                        deriv.dist_deriv( atom1, atom2, tdist2_CB, f1, f2 );
                        scal_f1_res1=scal(f1,dfunc);
                        scal_f2_res1=scal(f2,dfunc);
                        scal_f1_res2=scal(f1,-dfunc);
                        scal_f2_res2=scal(f2,-dfunc);
                        atom_derivatives_f1[i][4]=addv(atom_derivatives_f1[i][4],scal_f1_res1);
                        atom_derivatives_f2[i][4]=addv(atom_derivatives_f2[i][4],scal_f2_res1);
                        atom_derivatives_f1[j][4]=addv(atom_derivatives_f1[j][4],scal_f1_res2);
                        atom_derivatives_f2[j][4]=addv(atom_derivatives_f2[j][4],scal_f2_res2);
                    }
                }
                else if(tdist2_CB<80) //10<r<80
                {
                    enelist[15]+=weight_cont*weight*(1+sin((tdist2_CB-dc)/dd*PI))/2;

                    if(flag_grad)
                    {
                        dfunc=weight_cont*weight*(cos((tdist2_CB-dc)/dd*PI))/(2*dd*PI);
                        point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                        point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                        point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                        atom1.x=decstr[i].ptb.x,atom1.y=decstr[i].ptb.y,atom1.z=decstr[i].ptb.z;
                        atom2.x=decstr[j].ptb.x,atom2.y=decstr[j].ptb.y,atom2.z=decstr[j].ptb.z;
                        deriv.dist_deriv( atom1, atom2, tdist2_CB, f1, f2 );
                        scal_f1_res1=scal(f1,dfunc);
                        scal_f2_res1=scal(f2,dfunc);
                        scal_f1_res2=scal(f1,-dfunc);
                        scal_f2_res2=scal(f2,-dfunc);
                        atom_derivatives_f1[i][4]=addv(atom_derivatives_f1[i][4],scal_f1_res1);
                        atom_derivatives_f2[i][4]=addv(atom_derivatives_f2[i][4],scal_f2_res1);
                        atom_derivatives_f1[j][4]=addv(atom_derivatives_f1[j][4],scal_f1_res2);
                        atom_derivatives_f2[j][4]=addv(atom_derivatives_f2[j][4],scal_f2_res2);
                    }
                }
                else //tdist2_CB>=80
                {
                    enelist[15]+=weight_cont*weight;
                }
            }

            double weight=0.0;
            if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
            {
                weight=weights[0];
            }
            else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
            {
                weight=weights[1];
            }
            else if ( abs(i-j)>=LONG_RANGE_MIN )
            {
                weight=weights[2];
            }
            else
            {
                weight=0.0;
            }

            if ( tdist2_CB<19.75 && (decstr[i].aaa!='G' && decstr[j].aaa!='G') && weight>0.0 
                 && ( cb_dist_prob[i][j]<dist_cut_cb ) )
            {
                enelist[12] += weight * spline_dist[i][j].cubic_spline( tdist2_CB );
                if(flag_grad)
                {
                    double dfunc = weight * spline_dist[i][j].cubic_dspline( tdist2_CB );
                    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                    point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                    point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                    atom1.x=decstr[i].ptb.x,atom1.y=decstr[i].ptb.y,atom1.z=decstr[i].ptb.z;
                    atom2.x=decstr[j].ptb.x,atom2.y=decstr[j].ptb.y,atom2.z=decstr[j].ptb.z;
                    deriv.dist_deriv( atom1, atom2, tdist2_CB, f1, f2 );
                    scal_f1_res1=scal(f1,dfunc);
                    scal_f2_res1=scal(f2,dfunc);
                    scal_f1_res2=scal(f1,-dfunc);
                    scal_f2_res2=scal(f2,-dfunc);
                    atom_derivatives_f1[i][4]=addv(atom_derivatives_f1[i][4],scal_f1_res1);
                    atom_derivatives_f2[i][4]=addv(atom_derivatives_f2[i][4],scal_f2_res1);
                    atom_derivatives_f1[j][4]=addv(atom_derivatives_f1[j][4],scal_f1_res2);
                    atom_derivatives_f2[j][4]=addv(atom_derivatives_f2[j][4],scal_f2_res2);
                }
            }


            /**** CA-CA specific energy ****/
            tdist2_CA=CA_SC_tmp[i][j];
            tdist_CA=tdist2_CA*tdist2_CA;

            double weight_cont_ca=0.0;
            if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
            {
                weight_cont_ca=weights[9];
            }
            else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
            {
                weight_cont_ca=weights[10];
            }
            else if ( abs(i-j)>=LONG_RANGE_MIN )
            {
                weight_cont_ca=weights[11];
            }
            else
            {
                weight_cont_ca=0.0;
            }

            if( PcontactCA[i][j]>cont_ca_cut && weight_cont_ca>0.0 )
            {
                double dfunc=0.0;
                weight=PcontactCA[i][j];    //weight for well depth
                if(tdist2_CA<=d8)
                {	//r<8A
                    enelist[15]+=-weight*weight_cont_ca; // list data when weight=1
                }
                else if(tdist2_CA<d10) // 8<r<10
                {
                    enelist[15]+=-weight_cont_ca*weight*(1-sin((tdist2_CA-da)/db*PI))/2;

                    if(flag_grad)
                    {
                        dfunc=weight_cont_ca*weight*cos((tdist2_CA-da)/db*PI)/(2*db*PI);
                        point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                        point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                        point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                        atom1.x=decstr[i].x,atom1.y=decstr[i].y,atom1.z=decstr[i].z;
                        atom2.x=decstr[j].x,atom2.y=decstr[j].y,atom2.z=decstr[j].z;
                        deriv.dist_deriv( atom1, atom2, tdist2_CA, f1, f2 );
                        scal_f1_res1=scal(f1,dfunc);
                        scal_f2_res1=scal(f2,dfunc);
                        scal_f1_res2=scal(f1,-dfunc);
                        scal_f2_res2=scal(f2,-dfunc);
                        atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scal_f1_res1);
                        atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scal_f2_res1);
                        atom_derivatives_f1[j][0]=addv(atom_derivatives_f1[j][0],scal_f1_res2);
                        atom_derivatives_f2[j][0]=addv(atom_derivatives_f2[j][0],scal_f2_res2);
                    }
                }
                else if(tdist2_CA<80) //10<r<80
                {
                    enelist[15]+=weight_cont_ca*weight*(1+sin((tdist2_CA-dc)/dd*PI))/2;

                    if(flag_grad)
                    {
                        dfunc=weight_cont_ca*weight*(cos((tdist2_CA-dc)/dd*PI))/(2*dd*PI);
                        point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                        point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                        point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                        atom1.x=decstr[i].x,atom1.y=decstr[i].y,atom1.z=decstr[i].z;
                        atom2.x=decstr[j].x,atom2.y=decstr[j].y,atom2.z=decstr[j].z;
                        deriv.dist_deriv( atom1, atom2, tdist2_CA, f1, f2 );
                        scal_f1_res1=scal(f1,dfunc);
                        scal_f2_res1=scal(f2,dfunc);
                        scal_f1_res2=scal(f1,-dfunc);
                        scal_f2_res2=scal(f2,-dfunc);
                        atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scal_f1_res1);
                        atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scal_f2_res1);
                        atom_derivatives_f1[j][0]=addv(atom_derivatives_f1[j][0],scal_f1_res2);
                        atom_derivatives_f2[j][0]=addv(atom_derivatives_f2[j][0],scal_f2_res2);
                    }
                }
                else //tdist2_CA>=80
                {
                    enelist[15]+=weight_cont_ca*weight;
                }
            }

            weight=0.0;
            if ( abs(i-j)>=SHORT_RANGE_MIN && abs(i-j)<=SHORT_RANGE_MAX )
            {
                weight=weights[3];
            }
            else if ( abs(i-j)>=MED_RANGE_MIN && abs(i-j)<=MED_RANGE_MAX )
            {
                weight=weights[4];
            }
            else if ( abs(i-j)>=LONG_RANGE_MIN )
            {
                weight=weights[5];
            }
            else
            {
                weight=0.0;
            }

            if ( tdist2_CA<19.75 && weight>0.0 && ( ca_dist_prob[i][j]<dist_cut_ca ) )
            {
                enelist[8] += weight * spline_dist_ca[i][j].cubic_spline( tdist2_CA );
                if(flag_grad)
                {
                    double dfunc = weight * spline_dist_ca[i][j].cubic_dspline( tdist2_CA );
                    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                    point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                    point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                    atom1.x=decstr[i].x,atom1.y=decstr[i].y,atom1.z=decstr[i].z;
                    atom2.x=decstr[j].x,atom2.y=decstr[j].y,atom2.z=decstr[j].z;
                    deriv.dist_deriv( atom1, atom2, tdist2_CA, f1, f2 );
                    scal_f1_res1=scal(f1,dfunc);
                    scal_f2_res1=scal(f2,dfunc);
                    scal_f1_res2=scal(f1,-dfunc);
                    scal_f2_res2=scal(f2,-dfunc);
                    atom_derivatives_f1[i][0]=addv(atom_derivatives_f1[i][0],scal_f1_res1);
                    atom_derivatives_f2[i][0]=addv(atom_derivatives_f2[i][0],scal_f2_res1);
                    atom_derivatives_f1[j][0]=addv(atom_derivatives_f1[j][0],scal_f1_res2);
                    atom_derivatives_f2[j][0]=addv(atom_derivatives_f2[j][0],scal_f2_res2);
                }
            }


            //---- CA contacts ---->
            /**** SC-SC specific energy ****/
            tdist2_SC=CA_SC_tmp[j][i];
            tdist_SC=tdist2_SC*tdist2_SC;

            //---- I-TASSER restraints (only for QT/QTC) ---->
            /*
            if(flagcont and tasspair[indm].ispair[2]) //comb.dat <----
            {
                tdist3=tdist2_SC-sgdavg[ind1][ind2]-2.5*sgdstd[ind1][ind2];//r-[d(Ai,Aj)+2.5*d(Ai,Aj)]
                //E=E+conf*[r-(d+2.5*dev(i,j)] if(r>d(Ai,Aj)+2.5*dev(Ai,Aj), similar to ITasser-old
                if(tdist3>0) tcont[2]+=tasspair[indm].pval[2]*tdist3; //E=conf*[r-(d+2.5sigma)]
            }
            */

            /* If both CA-CA distance (tdist2_CA) and side-chain center
             * distance (tdist2_SC) are so far away that no atoms are within
             * 15, skip full atomic calculation. CA-SG bond length is in the
             * range of [0,4.75663); CA-O bond length is 2.3985. Thus, full
             * atom energy can be skipped if two residues are apart by:
             * tdist2_CA > 15 + 2*2.3985 = 19.797       and 
             * tdist2_SC > 15
             * Solvation can be skipped if
             * max{tdist2_CA,tdist2_SC} > sqrt(9.2^2+4.75633^2) = 10.3568 */
            if ( tdist2_CA>20.0 && tdist2_SC>15.0 ) continue;

            /**** residue center-center specific energy ****/
            //---- solvation ---->
            /*
            sqdist=CB_CT_tmp[j][i];
            if(sqdist<=sqdistcut)
            {
                tdists=scalf/sqdist;
                soldat[i]+=tdists*asarea2[ind2];
                soldat[j]+=tdists*asarea2[ind1];
            }
            */
            /**** other backbone and side chain atom energy terms ****/
            for (ii=0;ii<8;ii++)
            {
                if (ii==5 || ii==6) continue;
                for (jj=0;jj<8;jj++)
                {
                    if (jj==5 || jj==6) continue;

                    if (ii==0 && jj==0) // CA-CA
                    {
                        tdist=tdist_CA;
                        tdist2=tdist2_CA;
                    }
                    else if (ii==7 && jj==7) // SC-SC
                    {
                        tdist=tdist_SC;
                        tdist2=tdist2_SC;
                    }
                    else if (ii==4 && jj==4) //&& PredictCB[i][j]>0.0) // CB-CB
                    {
                        tdist=tdist_CB;
                        tdist2=tdist2_CB;
                    }
                    else
                    {
                        //tp=bf.minu(pin[ii],pjn[jj]);
                        tp.x=pn[r1+ii][0]-pn[r2+jj][0];
                        tp.y=pn[r1+ii][1]-pn[r2+jj][1];
                        tp.z=pn[r1+ii][2]-pn[r2+jj][2];
                        tdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z; //r*r
                        tdist2=sqrt(tdist); //r
                    }
                    if (tdist2>=15.0) continue;
      
                    //---- excluded volume ---->
                    tdist-=vdwds[ii][jj];
                    if(tdist<0)
                    {
                        enelist[1] -= ( 0.03 * weights[24] ) * tdist;
                        
                        if(flag_grad)
                        {
                            double dfunc = weights[24]*0.03*-2.0*tdist2; //todo get rid of magic number weights
                            point3d f1, f2;
                            point3d scal_f1_res1,scal_f2_res1,scal_f1_res2,scal_f2_res2;
                            atom1.x=pn[r1+ii][0],atom1.y=pn[r1+ii][1],atom1.z=pn[r1+ii][2];
                            atom2.x=pn[r2+jj][0],atom2.y=pn[r2+jj][1],atom2.z=pn[r2+jj][2];
                            deriv.dist_deriv( atom1, atom2, tdist2, f1, f2 );
                            scal_f1_res1=scal(f1,dfunc);
                            scal_f2_res1=scal(f2,dfunc);
                            scal_f1_res2=scal(f1,-dfunc);
                            scal_f2_res2=scal(f2,-dfunc);
                            atom_derivatives_f1[i][ii]=addv(atom_derivatives_f1[i][ii],scal_f1_res1);
                            atom_derivatives_f2[i][ii]=addv(atom_derivatives_f2[i][ii],scal_f2_res1);
                            atom_derivatives_f1[j][jj]=addv(atom_derivatives_f1[j][jj],scal_f1_res2);
                            atom_derivatives_f2[j][jj]=addv(atom_derivatives_f2[j][jj],scal_f2_res2);
                        }
                    }

                    //------- RW potential ----->
                    if( ii<5 && jj<5 && tdist2<15.0 && !( ii==4 && ind1==5 ) &&
                        !( jj==4 && ind2==5 ) && weights[26]>0.0 )
                    {   //backbone atoms, ii=4(CB), ind1=5(GLY)
                        enelist[4]+=weights[26]*spline_rw[ind5[ii]][ind6[jj]].cubic_spline( tdist2 );

                        if(flag_grad)
                        {
                            double dfunc = weights[26]*spline_rw[ind5[ii]][ind6[jj]].cubic_dspline( tdist2 );
                            point3d f1, f2;
                            point3d scal_f1_res1,scal_f2_res1,scal_f1_res2,scal_f2_res2;
                            atom1.x=pn[r1+ii][0],atom1.y=pn[r1+ii][1],atom1.z=pn[r1+ii][2];
                            atom2.x=pn[r2+jj][0],atom2.y=pn[r2+jj][1],atom2.z=pn[r2+jj][2];
                            deriv.dist_deriv( atom1, atom2, tdist2, f1, f2 );
                            scal_f1_res1=scal(f1,dfunc);
                            scal_f2_res1=scal(f2,dfunc);
                            scal_f1_res2=scal(f1,-dfunc);
                            scal_f2_res2=scal(f2,-dfunc);
                            atom_derivatives_f1[i][ii]=addv(atom_derivatives_f1[i][ii],scal_f1_res1);
                            atom_derivatives_f2[i][ii]=addv(atom_derivatives_f2[i][ii],scal_f2_res1);
                            atom_derivatives_f1[j][jj]=addv(atom_derivatives_f1[j][jj],scal_f1_res2);
                            atom_derivatives_f2[j][jj]=addv(atom_derivatives_f2[j][jj],scal_f2_res2);
                        }
                    }

                    //------ pairwise side-chain center atomic potential ----->
                    if(ii==7 && (jj<4 || jj==7) && tdist2<=15.0 && weights[27]>0.0 )
                    { //SG-backbone and SG-SG
                        enelist[5]+=weights[27]*wtsg*spline_poladat[ind1][20*sgmap[jj]+ind2].cubic_spline( tdist2 );

                        if(flag_grad)
                        {
                            double dfunc = weights[27]*wtsg*spline_poladat[ind1][20*sgmap[jj]+ind2].cubic_dspline( tdist2 );
                            point3d f1, f2;
                            point3d scal_f1_res1,scal_f2_res1,scal_f1_res2,scal_f2_res2;
                            atom1.x=pn[r1+ii][0],atom1.y=pn[r1+ii][1],atom1.z=pn[r1+ii][2];
                            atom2.x=pn[r2+jj][0],atom2.y=pn[r2+jj][1],atom2.z=pn[r2+jj][2];
                            deriv.dist_deriv( atom1, atom2, tdist2, f1, f2 );
                            scal_f1_res1=scal(f1,dfunc);
                            scal_f2_res1=scal(f2,dfunc);
                            scal_f1_res2=scal(f1,-dfunc);
                            scal_f2_res2=scal(f2,-dfunc);
                            atom_derivatives_f1[i][ii]=addv(atom_derivatives_f1[i][ii],scal_f1_res1);
                            atom_derivatives_f2[i][ii]=addv(atom_derivatives_f2[i][ii],scal_f2_res1);
                            atom_derivatives_f1[j][jj]=addv(atom_derivatives_f1[j][jj],scal_f1_res2);
                            atom_derivatives_f2[j][jj]=addv(atom_derivatives_f2[j][jj],scal_f2_res2);
                        }
                    }
                    if((ii<4 || ii==7) && jj==7 && tdist2<=15.0 && weights[27]>0.0 )
                    {
                        enelist[5]+=weights[27]*wtsg*spline_poladat[ind2][20*sgmap[ii]+ind1].cubic_spline( tdist2 );

                        if(flag_grad)
                        {
                            double dfunc = weights[27]*wtsg*spline_poladat[ind2][20*sgmap[ii]+ind1].cubic_dspline( tdist2 );
                            point3d f1, f2;
                            point3d scal_f1_res1,scal_f2_res1,scal_f1_res2,scal_f2_res2;
                            atom1.x=pn[r1+ii][0],atom1.y=pn[r1+ii][1],atom1.z=pn[r1+ii][2];
                            atom2.x=pn[r2+jj][0],atom2.y=pn[r2+jj][1],atom2.z=pn[r2+jj][2];
                            deriv.dist_deriv( atom1, atom2, tdist2, f1, f2 );
                            scal_f1_res1=scal(f1,dfunc);
                            scal_f2_res1=scal(f2,dfunc);
                            scal_f1_res2=scal(f1,-dfunc);
                            scal_f2_res2=scal(f2,-dfunc);
                            atom_derivatives_f1[i][ii]=addv(atom_derivatives_f1[i][ii],scal_f1_res1);
                            atom_derivatives_f2[i][ii]=addv(atom_derivatives_f2[i][ii],scal_f2_res1);
                            atom_derivatives_f1[j][jj]=addv(atom_derivatives_f1[j][jj],scal_f1_res2);
                            atom_derivatives_f2[j][jj]=addv(atom_derivatives_f2[j][jj],scal_f2_res2);
                        }
                    }
                } // end of jj
            } //end of ii
        } //j
        //soldat[i]+=scalf*asarea2[ind1]/16.00;
    } //i


    //------ re-weight some energy terms ---------->
    if(dptype==1) // distProf from init.dat for QP/QT
        wtdp=wtinitdp;//if using init.dat


    //------- Radius of gyration --------------->
    if ( weights[28]>0.0 )
    {
        double trad=0.0;
        tcen.x/=double(numseq);tcen.y/=double(numseq);tcen.z/=double(numseq);
        for(i=0;i<numseq;i++)
        {
            r1=i*8;
            tp.x=pn[r1][0]-tcen.x;
            tp.y=pn[r1][1]-tcen.y;
            tp.z=pn[r1][2]-tcen.z;
            tdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z;
            trad+=tdist;
        }
        trad=sqrt(trad/numseq);
        double minradius=exp(0.840)*pow(numseq,0.358)-0.5;
        double maxradius=minradius+8.0;
        double tmaxval=longestdist*0.5*sqrt(0.6);
        if ( tmaxval > maxradius ) maxradius=tmaxval;
        if ( trad > maxradius )
        {
            enelist[7] = weights[28] * ( trad - maxradius ) * ( trad - maxradius );
            if( flag_grad )
            {
                for( i=0 ; i<numseq ; i++ )
                {
                    r1=i*8;
                    point3d drg_dx,f1,f2,temp,atom1;
                    atom1.x = pn[r1][0], atom1.y = pn[r1][1], atom1.z = pn[r1][2];
                    drg_dx.x = ( trad - maxradius ) * 2.0 * ( pn[r1][0] - tcen.x ) / ( trad * numseq );
                    drg_dx.y = ( trad - maxradius ) * 2.0 * ( pn[r1][1] - tcen.y ) / ( trad * numseq );
                    drg_dx.z = ( trad - maxradius ) * 2.0 * ( pn[r1][2] - tcen.z ) / ( trad * numseq );
                    f2 = scal( drg_dx, weights[28] );
                    atom_derivatives_f2[i][0] = addv( atom_derivatives_f2[i][0], f2 );
                    temp = scal( drg_dx, -1.0 );
                    f1 = prod( atom1, addv( temp, atom1 ) );
                    f1 = scal( f1, weights[28] );
                    atom_derivatives_f1[i][0] = addv( atom_derivatives_f1[i][0], f1 );
                }
            }
        }
        else if(trad<minradius)
        { 
            enelist[7] = weights[28] * (minradius-trad) * (minradius-trad);

            if( flag_grad )
            {
                for( i=0 ; i<numseq ; i++ )
                {
                    r1=i*8;
                    point3d drg_dx,f1,f2,temp,atom1;
                    atom1.x = pn[r1][0], atom1.y = pn[r1][1], atom1.z = pn[r1][2];
                    drg_dx.x = ( minradius - trad ) * -2.0 * ( pn[r1][0] - tcen.x ) / ( trad * numseq );
                    drg_dx.y = ( minradius - trad ) * -2.0 * ( pn[r1][1] - tcen.y ) / ( trad * numseq );
                    drg_dx.z = ( minradius - trad ) * -2.0 * ( pn[r1][2] - tcen.z ) / ( trad * numseq );
                    f2 = scal( drg_dx, weights[28] );
                    atom_derivatives_f2[i][0] = addv( atom_derivatives_f2[i][0], f2 );
                    temp = scal( drg_dx, -1.0 );
                    f1 = prod( atom1, addv( temp, atom1 ) );
                    f1 = scal( f1, weights[28] );
                    atom_derivatives_f1[i][0] = addv( atom_derivatives_f1[i][0], f1 );
                }
            }
        }
    }

    //---------- Energy for I-TASSER restraints --------------->
    //if(flagcont) enelist[0]+=wtcont[0]*tcont[0]  // dist
    //                       + wtcont[1]*tcont[1]  // distL
    //                       + wtcont[2]*tcont[2]  // comb
    //                       + wtcont[3]*tcont[3]  // combCA
    //                       + wtcont[4]*tcont[4]; // comb8CA
  
    //------- solvation potential - not currently instated -------------->
    /*
    for(i=0;i<numseq;i++)
    {
        soldat[i]*=scalaa[decstr[i].iaa];
        if(soldat[i]>1.0) soldat[i]=1.0;
        soldat[i]=1.0-soldat[i];//area from 0 to 1
        tdists=solseq[i]-soldat[i];
        if(tdists<0.0){
            tdists*=-1.0;
	}
        enelist[6]+=tdists;
    }
    enelist[6]*=4.00;
    delete[]soldat;
    soldat=NULL;
    */
}

void Energy::energy_hbond( point3f *decstr, double *enelist )
{
    const double HBOND_DISTANCE_CUTOFF_MAX = 3.0;
    const double HBOND_OPTIMAL_DISTANCE    = 1.9;
    const double HBOND_WELL_DEPTH = 1.0;
   
    int r1,r2;
    double dist_HA_ij,angle_DHA_ij,angle_HAB_ij;
    point3d p1,p2;
    point3d diss;
    for(int i=0;i<numseq-1;i++)
    {
        r1=i*8;
        for(int j=i+1;j<numseq;j++)
        {
            r2=j*8;
            
            { 
                p1 = setv(pn[r1+5][0],pn[r1+5][1],pn[r1+5][2]);
                p2 = setv(pn[r2+3][0],pn[r2+3][1],pn[r2+3][2]);
                dist_HA_ij = norm( minu( p2, p1 ) );
                if( dist_HA_ij>HBOND_DISTANCE_CUTOFF_MAX ) continue;
	    
                p1.x=pn[r1+1][0]-pn[r1+5][0];
                p1.y=pn[r1+1][1]-pn[r1+5][1];
                p1.z=pn[r1+1][2]-pn[r1+5][2];
                p2.x=pn[r2+3][0]-pn[r1+5][0];
                p2.y=pn[r2+3][1]-pn[r1+5][1];
                p2.z=pn[r2+3][2]-pn[r1+5][2];
                angle_DHA_ij = angv(p1,p2)*degrad;
                if( angle_DHA_ij<90 ) continue;
                
                p1.x=pn[r1+5][0]-pn[r2+3][0];
                p1.y=pn[r1+5][1]-pn[r2+3][1];
                p1.z=pn[r1+5][2]-pn[r2+3][2];
                p2.x=pn[r2+2][0]-pn[r2+3][0];
                p2.y=pn[r2+2][1]-pn[r2+3][1];
                p2.z=pn[r2+2][2]-pn[r2+3][2];
                angle_HAB_ij = angv(p1,p2)*degrad;
                if( angle_HAB_ij<80 ) continue;
            
                angle_DHA_ij*=raddeg;
                angle_HAB_ij*=raddeg;

                double energyR=0.0;
                if(dist_HA_ij<HBOND_OPTIMAL_DISTANCE)
                {
                    energyR=-1.0*HBOND_WELL_DEPTH*cos((dist_HA_ij-HBOND_OPTIMAL_DISTANCE)*PI);
                }
                else
                {
                    energyR=-0.5*cos(PI/(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE)*(dist_HA_ij-HBOND_OPTIMAL_DISTANCE))-0.5;
                }
                if(energyR > 0.0) energyR = 0.0;
 
                double cos_angle_DHA_ij = cos(angle_DHA_ij);
                double cos_angle_HAB_ij = cos(angle_HAB_ij-(150*raddeg));
                double energyTheta = -1.0*cos_angle_DHA_ij*cos_angle_DHA_ij*cos_angle_DHA_ij*cos_angle_DHA_ij;
                double energyPhi = -1.0*cos_angle_HAB_ij*cos_angle_HAB_ij*cos_angle_HAB_ij*cos_angle_HAB_ij;
                double energy = energyR + energyTheta + energyPhi;
                if(energy>0.0) energy=0.0;
                enelist[2]+=weights[25]*energy;
	    
                if(flag_grad && energy<0.0)
                {
                    double dfunc;
                    point3d atom1,atom2,atom3,atom4;
                    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                    point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                    point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                
                    if(energyR<0.0)
                    {
                        if(dist_HA_ij<HBOND_OPTIMAL_DISTANCE)
                        {
                            dfunc=1.0*HBOND_WELL_DEPTH*PI*sin((dist_HA_ij-HBOND_OPTIMAL_DISTANCE)*PI);
                        }
                        else
                        {
                            dfunc=0.5*PI*sin(PI/(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE)
                                 *(dist_HA_ij-HBOND_OPTIMAL_DISTANCE))
                                 /(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE);
                        }
                        dfunc*=weights[25];

                        atom1 = setv(pn[r1+5][0],pn[r1+5][1],pn[r1+5][2]);
                        atom2 = setv(pn[r2+3][0],pn[r2+3][1],pn[r2+3][2]);
                        deriv.dist_deriv( atom1, atom2, dist_HA_ij, f1, f2 );
                        scal_f1_res1=scal(f1,dfunc);
                        scal_f2_res1=scal(f2,dfunc);
                        scal_f1_res2=scal(f1,-dfunc);
                        scal_f2_res2=scal(f2,-dfunc);
                        atom_derivatives_f1[i][5]=addv(atom_derivatives_f1[i][5],scal_f1_res1);
                        atom_derivatives_f2[i][5]=addv(atom_derivatives_f2[i][5],scal_f2_res1);
                        atom_derivatives_f1[j][3]=addv(atom_derivatives_f1[j][3],scal_f1_res2);
                        atom_derivatives_f2[j][3]=addv(atom_derivatives_f2[j][3],scal_f2_res2);
                    }

                    dfunc=weights[25]*4*cos_angle_DHA_ij*cos_angle_DHA_ij*cos_angle_DHA_ij*sin(angle_DHA_ij);
                    point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                    point3d f1_p1, f2_p1, f1_p2, f2_p2, f1_p3, f2_p3;
                    point3d scalf1_p1, scalf2_p1, scalf1_p2, scalf2_p2, scalf1_p3, scalf2_p3;
                    double theta=0.0;
                    atom1=setv(pn[r1+1][0],pn[r1+1][1],pn[r1+1][2]);  //Ca i
                    atom2=setv(pn[r1+5][0],pn[r1+5][1],pn[r1+5][2]);  //Cb i
                    atom3=setv(pn[r2+3][0],pn[r2+3][1],pn[r2+3][2]);  //Cb j
                    deriv.angular_deriv(atom1,atom2,atom3,f1_p1,f2_p1,f1_p2,f2_p2,f1_p3,f2_p3);
                    scalf1_p1=scal(f1_p1,dfunc);
                    scalf2_p1=scal(f2_p1,dfunc);
                    atom_derivatives_f1[i][1]=addv(atom_derivatives_f1[i][1],scalf1_p1);
                    atom_derivatives_f2[i][1]=addv(atom_derivatives_f2[i][1],scalf2_p1);
                    scalf1_p2=scal(f1_p2,dfunc);
                    scalf2_p2=scal(f2_p2,dfunc);
                    atom_derivatives_f1[i][5]=addv(atom_derivatives_f1[i][5],scalf1_p2);
                    atom_derivatives_f2[i][5]=addv(atom_derivatives_f2[i][5],scalf2_p2);
                    scalf1_p3=scal(f1_p3,dfunc);
                    scalf2_p3=scal(f2_p3,dfunc);
                    atom_derivatives_f1[j][3]=addv(atom_derivatives_f1[j][3],scalf1_p3);
                    atom_derivatives_f2[j][3]=addv(atom_derivatives_f2[j][3],scalf2_p3);

                    dfunc=weights[25]*4*cos_angle_HAB_ij*cos_angle_HAB_ij*cos_angle_HAB_ij
                         *sin(angle_HAB_ij-(150*raddeg));
                    atom1=setv(pn[r1+5][0],pn[r1+5][1],pn[r1+5][2]);  //Ca i
                    atom2=setv(pn[r2+3][0],pn[r2+3][1],pn[r2+3][2]);  //Cb i
                    atom3=setv(pn[r2+2][0],pn[r2+2][1],pn[r2+2][2]);  //Cb j
                    deriv.angular_deriv(atom1,atom2,atom3,f1_p1,f2_p1,f1_p2,f2_p2,f1_p3,f2_p3);
                    scalf1_p1=scal(f1_p1,dfunc);
                    scalf2_p1=scal(f2_p1,dfunc);
                    atom_derivatives_f1[i][5]=addv(atom_derivatives_f1[i][5],scalf1_p1);
                    atom_derivatives_f2[i][5]=addv(atom_derivatives_f2[i][5],scalf2_p1);
                    scalf1_p2=scal(f1_p2,dfunc);
                    scalf2_p2=scal(f2_p2,dfunc);
                    atom_derivatives_f1[j][3]=addv(atom_derivatives_f1[j][3],scalf1_p2);
                    atom_derivatives_f2[j][3]=addv(atom_derivatives_f2[j][3],scalf2_p2);
                    scalf1_p3=scal(f1_p3,dfunc);
                    scalf2_p3=scal(f2_p3,dfunc);
                    atom_derivatives_f1[j][2]=addv(atom_derivatives_f1[j][2],scalf1_p3);
                    atom_derivatives_f2[j][2]=addv(atom_derivatives_f2[j][2],scalf2_p3);
                }
            }

            {
                p1 = setv(pn[r2+5][0],pn[r2+5][1],pn[r2+5][2]);
                p2 = setv(pn[r1+3][0],pn[r1+3][1],pn[r1+3][2]);
                dist_HA_ij = norm( minu( p2, p1 ) );
                if( dist_HA_ij>HBOND_DISTANCE_CUTOFF_MAX ) continue;
	    
                p1.x=pn[r2+1][0]-pn[r2+5][0];
                p1.y=pn[r2+1][1]-pn[r2+5][1];
                p1.z=pn[r2+1][2]-pn[r2+5][2];
                p2.x=pn[r1+3][0]-pn[r2+5][0];
                p2.y=pn[r1+3][1]-pn[r2+5][1];
                p2.z=pn[r1+3][2]-pn[r2+5][2];
                angle_DHA_ij = angv(p1,p2)*degrad;
                if( angle_DHA_ij<90 ) continue;
                
                p1.x=pn[r2+5][0]-pn[r1+3][0];
                p1.y=pn[r2+5][1]-pn[r1+3][1];
                p1.z=pn[r2+5][2]-pn[r1+3][2];
                p2.x=pn[r1+2][0]-pn[r1+3][0];
                p2.y=pn[r1+2][1]-pn[r1+3][1];
                p2.z=pn[r1+2][2]-pn[r1+3][2];
                angle_HAB_ij = angv(p1,p2)*degrad;
                if( angle_HAB_ij<80 ) continue;
            
                angle_DHA_ij*=raddeg;
                angle_HAB_ij*=raddeg;

                double energyR=0.0;
                if(dist_HA_ij<HBOND_OPTIMAL_DISTANCE)
                {
                    energyR=-1.0*HBOND_WELL_DEPTH*cos((dist_HA_ij-HBOND_OPTIMAL_DISTANCE)*PI);
                }
                else
                {
                    energyR=-0.5*cos(PI/(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE)*(dist_HA_ij-HBOND_OPTIMAL_DISTANCE))-0.5;
                }
                if(energyR > 0.0) energyR = 0.0;
 
                double cos_angle_DHA_ij = cos(angle_DHA_ij);
                double cos_angle_HAB_ij = cos(angle_HAB_ij-(150*raddeg));
                double energyTheta = -1.0*cos_angle_DHA_ij*cos_angle_DHA_ij*cos_angle_DHA_ij*cos_angle_DHA_ij;
                double energyPhi = -1.0*cos_angle_HAB_ij*cos_angle_HAB_ij*cos_angle_HAB_ij*cos_angle_HAB_ij;
                double energy = energyR + energyTheta + energyPhi;
                if(energy>0.0) energy=0.0;
                enelist[2]+=weights[25]*energy;
	    
                if(flag_grad && energy<0.0)
                {
                    double dfunc;
                    point3d atom1,atom2,atom3,atom4;
                    point3d f1=setv(0.0,0.0,0.0),f2=setv(0.0,0.0,0.0);
                    point3d scal_f1_res1=setv(0.0,0.0,0.0),scal_f2_res1=setv(0.0,0.0,0.0);
                    point3d scal_f1_res2=setv(0.0,0.0,0.0),scal_f2_res2=setv(0.0,0.0,0.0);
                
                    if(energyR<0.0)
                    {
                        if(dist_HA_ij<HBOND_OPTIMAL_DISTANCE)
                        {
                            dfunc=1.0*HBOND_WELL_DEPTH*PI*sin((dist_HA_ij-HBOND_OPTIMAL_DISTANCE)*PI);
                        }
                        else
                        {
                            dfunc=0.5*PI*sin(PI/(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE)
                                 *(dist_HA_ij-HBOND_OPTIMAL_DISTANCE))
                                 /(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE);
                        }
                        dfunc*=weights[25];

                        atom1 = setv(pn[r2+5][0],pn[r2+5][1],pn[r2+5][2]);
                        atom2 = setv(pn[r1+3][0],pn[r1+3][1],pn[r1+3][2]);
                        deriv.dist_deriv( atom1, atom2, dist_HA_ij, f1, f2 );
                        scal_f1_res1=scal(f1,dfunc);
                        scal_f2_res1=scal(f2,dfunc);
                        scal_f1_res2=scal(f1,-dfunc);
                        scal_f2_res2=scal(f2,-dfunc);
                        atom_derivatives_f1[j][5]=addv(atom_derivatives_f1[j][5],scal_f1_res1);
                        atom_derivatives_f2[j][5]=addv(atom_derivatives_f2[j][5],scal_f2_res1);
                        atom_derivatives_f1[i][3]=addv(atom_derivatives_f1[i][3],scal_f1_res2);
                        atom_derivatives_f2[i][3]=addv(atom_derivatives_f2[i][3],scal_f2_res2);
                    }

                    dfunc=weights[25]*4*cos_angle_DHA_ij*cos_angle_DHA_ij*cos_angle_DHA_ij*sin(angle_DHA_ij);
                    point3d scal_f1=setv(0.0,0.0,0.0),scal_f2=setv(0.0,0.0,0.0);
                    point3d f1_p1, f2_p1, f1_p2, f2_p2, f1_p3, f2_p3;
                    point3d scalf1_p1, scalf2_p1, scalf1_p2, scalf2_p2, scalf1_p3, scalf2_p3;
                    atom1=setv(pn[r2+1][0],pn[r2+1][1],pn[r2+1][2]);  //Ca i
                    atom2=setv(pn[r2+5][0],pn[r2+5][1],pn[r2+5][2]);  //Cb i
                    atom3=setv(pn[r1+3][0],pn[r1+3][1],pn[r1+3][2]);  //Cb j
                    deriv.angular_deriv(atom1,atom2,atom3,f1_p1,f2_p1,f1_p2,f2_p2,f1_p3,f2_p3);
                    scalf1_p1=scal(f1_p1,dfunc);
                    scalf2_p1=scal(f2_p1,dfunc);
                    atom_derivatives_f1[j][1]=addv(atom_derivatives_f1[j][1],scalf1_p1);
                    atom_derivatives_f2[j][1]=addv(atom_derivatives_f2[j][1],scalf2_p1);
                    scalf1_p2=scal(f1_p2,dfunc);
                    scalf2_p2=scal(f2_p2,dfunc);
                    atom_derivatives_f1[j][5]=addv(atom_derivatives_f1[j][5],scalf1_p2);
                    atom_derivatives_f2[j][5]=addv(atom_derivatives_f2[j][5],scalf2_p2);
                    scalf1_p3=scal(f1_p3,dfunc);
                    scalf2_p3=scal(f2_p3,dfunc);
                    atom_derivatives_f1[i][3]=addv(atom_derivatives_f1[i][3],scalf1_p3);
                    atom_derivatives_f2[i][3]=addv(atom_derivatives_f2[i][3],scalf2_p3);

                    dfunc=weights[25]*4*cos_angle_HAB_ij*cos_angle_HAB_ij*cos_angle_HAB_ij
                         *sin(angle_HAB_ij-(150*raddeg));
                    atom1=setv(pn[r2+5][0],pn[r2+5][1],pn[r2+5][2]);  //Ca i
                    atom2=setv(pn[r1+3][0],pn[r1+3][1],pn[r1+3][2]);  //Cb i
                    atom3=setv(pn[r1+2][0],pn[r1+2][1],pn[r1+2][2]);  //Cb j
                    deriv.angular_deriv(atom1,atom2,atom3,f1_p1,f2_p1,f1_p2,f2_p2,f1_p3,f2_p3);
                    scalf1_p1=scal(f1_p1,dfunc);
                    scalf2_p1=scal(f2_p1,dfunc);
                    atom_derivatives_f1[j][5]=addv(atom_derivatives_f1[j][5],scalf1_p1);
                    atom_derivatives_f2[j][5]=addv(atom_derivatives_f2[j][5],scalf2_p1);
                    scalf1_p2=scal(f1_p2,dfunc);
                    scalf2_p2=scal(f2_p2,dfunc);
                    atom_derivatives_f1[i][3]=addv(atom_derivatives_f1[i][3],scalf1_p2);
                    atom_derivatives_f2[i][3]=addv(atom_derivatives_f2[i][3],scalf2_p2);
                    scalf1_p3=scal(f1_p3,dfunc);
                    scalf2_p3=scal(f2_p3,dfunc);
                    atom_derivatives_f1[i][2]=addv(atom_derivatives_f1[i][2],scalf1_p3);
                    atom_derivatives_f2[i][2]=addv(atom_derivatives_f2[i][2],scalf2_p3);
                }
            }
        }
    }
}

void Energy::torsion_energy(point3f *decstr, double *enelist ){
    
    if(weights[29]==0.0) return;
    
    for(int i=1;i<numseq-1;i++)
    {
        point3d atom1,atom2,atom3,atom4;
        double phi=decstr[i].phi*raddeg;       
        double psi=decstr[i+1].psi*raddeg;       
        enelist[14]+=weights[29]*(1-cos((phi-pred_phi[i])));
        enelist[14]+=weights[29]*(1-cos((psi-pred_psi[i])));
    }
}

void Energy::torsion_energy(point3f *decstr, double *enelist,double *dE_dvars ){
   
    if(weights[29]==0.0) return;
    
    for(int i=1;i<numseq-1;i++)
    {
        double phi=decstr[i].phi*raddeg;       
        double psi=decstr[i+1].psi*raddeg;       
        enelist[14]+=weights[29]*(1-cos((phi-pred_phi[i])));
        enelist[14]+=weights[29]*(1-cos((psi-pred_psi[i])));

        if(flag_grad)
        {
            dE_dvars[(i-1)*2] += weights[29]*sin((psi-pred_psi[i]));
            dE_dvars[(i-1)*2+1] += weights[29]*sin((phi-pred_phi[i]));
        }
    }
}

double Energy::calcrmsdenergy( point3f *decstr, double *vars, double *dE_dvars, bool calc_gradients )
{
    double trmsd=0.0; // total energy E_total = sum of enelist
    int i,j,ii;
    for(i=0;i<20;i++) enelist[i]=0;
    geo.apply_torsions( decstr, vars, numseq );
    flag_grad=calc_gradients;
    
    // Reset gradients with respect to each variable (torsions)
    if(flag_grad)
    {
        for(i=0;i<(numseq-1)*2;i++)
        {
            dE_dvars[i]=0.0;
        }

        // Initialize f1/f2 as zero
        for(i=0;i<numseq;i++)
        {
            for(ii=0;ii<8;ii++)
            {
                atom_derivatives_f1[i][ii].x=0.0;
                atom_derivatives_f1[i][ii].y=0.0;
                atom_derivatives_f1[i][ii].z=0.0;
                atom_derivatives_f2[i][ii].x=0.0;
                atom_derivatives_f2[i][ii].y=0.0;
                atom_derivatives_f2[i][ii].z=0.0;
            }
        }
    }

    // Calculate most atomic energy terms and their gradients
    calcallenergy(decstr, enelist);

    // Calculate torsion energy from predicted phi/psi
    torsion_energy(decstr, enelist,dE_dvars);

    // Calculate hydrogen bonding energy;
    if(weights[25]>0.0) energy_hbond(decstr,enelist);
    
    for(i=0;i<20;i++) trmsd+=enelist[i];
    
    // Sum all f1/f2 contributions from atoms moved by each given torsion angle. This can be done
    // recurrently but the time savings is probably minimal
    if ( flag_grad )
    {
        for ( i=1; i<numseq; i++ )
        {
            decstr[i].f1_psi.x=0.0,decstr[i].f1_psi.y=0.0,decstr[i].f1_psi.z=0.0;
            decstr[i].f2_psi.x=0.0,decstr[i].f2_psi.y=0.0,decstr[i].f2_psi.z=0.0;
            decstr[i].f1_phi.x=0.0,decstr[i].f1_phi.y=0.0,decstr[i].f1_phi.z=0.0;
            decstr[i].f2_phi.x=0.0,decstr[i].f2_phi.y=0.0,decstr[i].f2_phi.z=0.0;
            
            // Add intra-residue f1/f2 contributions
            for ( ii=0; ii<8; ii++ )
            {
                if(ii==1 || ii==5) continue;
                decstr[i].f1_phi=addv(decstr[i].f1_phi,atom_derivatives_f1[i][ii]);
                decstr[i].f2_phi=addv(decstr[i].f2_phi,atom_derivatives_f2[i][ii]);
            }

            // Add inter-residue f1/f2 contributions
            for ( j=i+1; j<numseq; j++ )
            {
                for ( ii=0; ii<8; ii++ )
                {
                    decstr[i].f1_phi=addv(decstr[i].f1_phi,atom_derivatives_f1[j][ii]);
                    decstr[i].f2_phi=addv(decstr[i].f2_phi,atom_derivatives_f2[j][ii]);
                }
            }

            // Set psi f1/f2 equal to phi since all atoms moves by phi are moved by psi
            decstr[i].f1_psi.x=decstr[i].f1_phi.x, decstr[i].f1_psi.y=decstr[i].f1_phi.y;
            decstr[i].f1_psi.z=decstr[i].f1_phi.z;
            decstr[i].f2_psi.x=decstr[i].f2_phi.x, decstr[i].f2_psi.y=decstr[i].f2_phi.y, 
            decstr[i].f2_psi.z=decstr[i].f2_phi.z;

            // Add in contributions from the backbone nitrogen and hydrogen to psi
            decstr[i].f1_psi=addv(decstr[i].f1_psi,atom_derivatives_f1[i][1]);
            decstr[i].f2_psi=addv(decstr[i].f2_psi,atom_derivatives_f2[i][1]);
            decstr[i].f1_psi=addv(decstr[i].f1_psi,atom_derivatives_f1[i][5]);
            decstr[i].f2_psi=addv(decstr[i].f2_psi,atom_derivatives_f2[i][5]);
        }
    }
    
    // Calculate gradients with respect to each torsion angle. 
    if( flag_grad )
    { 
        for( i=1; i<numseq; i++ )
        {
            double dphi=0.0,dpsi=0.0;
            point3d axis;
            point3d end_pos;
            end_pos.x = decstr[i-1].ptc.x;
            end_pos.y = decstr[i-1].ptc.y;
            end_pos.z = decstr[i-1].ptc.z;
            axis.x = decstr[i-1].ptc.x-decstr[i-1].x; 
            axis.y = decstr[i-1].ptc.y-decstr[i-1].y; 
            axis.z = decstr[i-1].ptc.z-decstr[i-1].z;
            axis = unit( axis );
            dpsi -= raddeg * ( dotv( axis, decstr[i].f1_psi ) + 
                               dotv( prod( axis, end_pos ), decstr[i].f2_psi ) );
            dE_dvars[(i-1)*2] = dpsi;

            end_pos.x = decstr[i].x;
            end_pos.y = decstr[i].y;
            end_pos.z = decstr[i].z;
            axis.x = decstr[i].x-decstr[i].ptn.x;
            axis.y = decstr[i].y-decstr[i].ptn.y;
            axis.z = decstr[i].z-decstr[i].ptn.z;
            axis = unit( axis );
            dphi -= raddeg * ( dotv( axis, decstr[i].f1_phi ) + 
                               dotv( prod( axis, end_pos ), decstr[i].f2_phi ) );
            dE_dvars[(i-1)*2+1] = dphi;
        }
    }
    
    return trmsd;
}

double Energy::calcrmsdenergy( point3f *decstr, double *vars )
{
    double trmsd=0.0; // total energy E_total = sum of enelist
    int i,j,ii;
    for(i=0;i<20;i++) enelist[i]=0;
    geo.apply_torsions( decstr, vars, numseq );
    flag_grad=false; // we don't need to calculate gradients here
    
    calcallenergy(decstr, enelist);
    torsion_energy(decstr, enelist);
    if(weights[25]>0.0) energy_hbond(decstr,enelist);
    
    for(i=0;i<20;i++) trmsd+=enelist[i];
  
    return trmsd;
}



#endif

