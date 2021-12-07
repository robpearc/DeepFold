/*******************************************************************************************
**  
**  Common parameters and data tables
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef COMMONPARAMETERS_H
#define COMMONPARAMETERS_H

#include <string>// TO use string variables
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

using namespace std;

#define PI 3.14159265358979323846
#define raddeg PI/180.0
#define degrad 180.0/PI
#define lennca 1.460f
#define delnca 0.004f*20
#define lencac 1.525f
#define delcac 0.004f*20
#define lencn 1.338f
#define delcn 0.005f*15
#define lennn 3.131f
#define lencc 3.231f
#define lencaca 3.813f
#define delcaca 0.019f*10
#define lencan1 2.441f
#define delcan1 0.036f*4
#define lencca1 2.446f
#define delcca1 0.036f*4
#define lennc 2.460f
#define delnc 0.012f*10
#define angncac 111.008f
#define angcacn 116.617f
#define angcnca 121.614f
#define MAX_LENGTH_ONE_LINE_IN_FILE  1024
#define MAX_ENERGY_TERM_NUM           100

//-------------------- Common Data Structures ----------------------->
typedef struct ssef
{
    char res;
    char ss;
    float a[3];//coil helix sheet
}ssef;//no change

typedef struct point3s
{
    float x,y,z;
}point3s;

typedef struct point3d{
        double x,y,z;
}point3d;

typedef struct point3f
{
    float x,y,z;
    float psi ,phi, omega;
    float len_c_n,len_n_ca,len_ca_c;
    float ang_ca_c_n,ang_c_n_ca,ang_n_ca_c;
    float angl,leng;
    point3s ptn,ptc,pto,ptb,pth,ptha,ptsg,ptg;
    point3d f1_psi,f2_psi,f1_phi,f2_phi;

    char ss2;//psipred in the fragment
    char ssm;//mycalculation
    char aaa; // one letter amino acid code (realseq ACD)
    unsigned char iaa;//realACD index
}point3f;

typedef struct sssegment
{
    int init;
    int term;
    char ss;
}sssegment;

//---------------------- Data Tables -------------------->

 //////////////0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35
static int Tab5[20][36] ={
    {  3,  2,  5,999,999,999,999,999,999,999,999,999,999,999,999,999,999,  1,999,999,999,999,999,999,999,999,  4,999,999,999,999,999,999,999,999,999},//0
    {  8,  7, 10,999,999,999,999,999,999,999,999,999,999,999,999,999,999,  6,999,999,999,999,999,999,999,999,  9,999,999,999,999,999,999,999, 11,999},//1
    { 14, 13, 16,999,999,999,999,999,999,999, 17,999,999,999,999,999,999, 12,999,999,999,999,999,999,999,999, 15, 18, 18,999,999,999,999,999,999,999},//2
    { 21, 20, 23, 25,999,999,999,999,999,999, 24,999,999,999,999,999,999, 19,999,999,999,999,999,999,999,999, 22,999,999, 26, 26,999,999,999,999,999},//3
    { 29, 28, 31,999, 33, 33,999, 34, 34,999, 32,999,999,999, 35,999,999, 27,999,999,999,999,999,999,999,999, 30,999,999,999,999,999,999,999,999,999},//4
    { 38, 37,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999, 36,999,999,999,999,999,999,999,999, 39,999,999,999,999,999,999,999,999,999},//5
    { 42, 41, 44,999,999, 47,999, 48,999,999, 45,999,999,999,999,999,999, 40, 46,999,999,999, 49,999,999,999, 43,999,999,999,999,999,999,999,999,999},//6
    { 52, 51, 54,999, 57,999,999,999,999,999,999, 55, 56,999,999,999,999, 50,999,999,999,999,999,999,999,999, 53,999,999,999,999,999,999,999,999,999},//7
    { 60, 59, 62, 64,999,999, 65,999,999,999, 63,999,999,999,999,999,999, 58,999,999,999,999,999,999,999, 66, 61,999,999,999,999,999,999,999,999,999},//8
    { 69, 68, 71,999, 73, 73,999,999,999,999, 72,999,999,999,999,999,999, 67,999,999,999,999,999,999,999,999, 70,999,999,999,999,999,999,999,999,999},//9
    { 76, 75, 78,999,999,999, 81,999,999,999, 79,999,999,999,999,999,999, 74,999,999,999,999,999,999,999,999, 77,999,999,999,999,999,999,999,999, 80},//10
    { 84, 83, 86,999,999,999,999,999,999,999, 87,999,999,999,999,999,999, 82,999, 89,999,999,999,999,999,999, 85, 88,999,999,999,999,999,999,999,999},//11
    { 92, 91, 94, 96,999,999,999,999,999,999, 95,999,999,999,999,999,999, 90,999,999,999,999,999,999,999,999, 93,999,999,999,999,999,999,999,999,999},//12
    { 99, 98,101,103,999,999,999,999,999,999,102,999,999,999,999,999,999, 97,999,999,999,999,105,999,999,999,100,999,999,104,999,999,999,999,999,999},//13
    {108,107,110,112,999,999,999,999,999,999,111,999,999,999,114,999,999,106,999,999,113,999,999,115,115,999,109,999,999,999,999,999,999,999,999,999},//14
    {118,117,120,999,999,999,999,999,999,999,999,999,999,999,999,999,999,116,999,999,999,999,999,999,999,999,119,999,999,999,999,121,999,999,999,999},//15
    {124,123,126,999,999,999,999,999,999,999,999,999,128,999,999,999,999,122,999,999,999,999,999,999,999,999,125,999,999,999,999,999,127,999,999,999},//16
    {131,130,133,999,999,999,999,999,999,999,999,134,134,999,999,999,999,129,999,999,999,999,999,999,999,999,132,999,999,999,999,999,999,999,999,999},//17
    {137,136,139,999,141,142,999,999,144,145,140,999,999,148,999,146,147,135,999,999,999,143,999,999,999,999,138,999,999,999,999,999,999,999,999,999},//18
    {151,150,153,999,155,155,999,156,156,999,154,999,999,999,157,999,999,149,999,999,999,999,999,999,999,999,152,999,999,999,999,999,999,158,999,999}};//19

//ca    n     c      o    cb    h     ha    sg
static	double vdwd[8][8]={
3.600,2.300,3.700,2.900,3.500,3.200,3.500,1.000,//ca
3.700,2.500,3.500,2.500,3.500,2.100,3.600,1.000,//n
2.300,1.200,2.700,2.700,2.300,1.800,2.200,1.000,//c
2.500,2.100,2.400,2.300,2.600,1.700,2.200,1.000,//o
3.500,2.300,3.500,2.800,3.300,2.300,3.400,1.000,//cb
3.500,2.200,3.000,1.800,3.500,2.200,3.500,1.000,//h
3.500,2.200,3.300,2.800,3.300,2.100,3.500,1.000,//ha
1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000};//sg ori 1.5

static	double vdwds[8][8]={
12.96000, 5.29000,13.69000, 8.41000,12.25000,10.24000,12.25000, 1.00000,
13.69000, 6.25000,12.25000, 6.25000,12.25000, 4.41000,12.96000, 1.00000,
 5.29000, 1.44000, 7.29000, 7.29000, 5.29000, 3.24000, 4.84000, 1.00000,
 6.25000, 4.41000, 5.76000, 5.29000, 6.76000, 2.89000, 4.84000, 1.00000,
12.25000, 5.29000,12.25000, 7.84000,10.89000, 5.29000,11.56000, 1.00000,
12.25000, 4.84000, 9.00000, 3.24000,12.25000, 4.84000,12.25000, 1.00000,
12.25000, 4.84000,10.89000, 7.84000,10.89000, 4.41000,12.25000, 1.00000,
 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000};

static	double nhochhmean[22]={
	5.2172,8.7074,230.1397,4.2198,6.5391
};

static	double nhochhsigma[22]={
	0.3676,0.4375,10.2234,0.4080,0.3763
};
	
static	double nhocmeanval[][4]={//h-o h-o=c n-h-o n-h-o=c  
	2.85,89.0,110.5,199.5,//i+3
	2.00,147.0,159.0,160.0, //i+4  oldoh 
//	2.83,89.0,110.0,201.5,//i+3 newoh
//	2.00,148.0,159.0,155.0, //i+4

	2.00,155.0,164.0,180.0,//0
	2.00,155.0,164.0,180.0,//1
	2.00,151.0,163.0,192.0,//2
	2.00,151.0,163.0,192.0,//3

};
	
static	double nhocstdval[][4]={//h-o h-o=c n-h-o n-h-o=c 
	0.315504,7.697185,8.980366,7.932107,
	0.530859,10.582243,11.249764,25.360054,

	0.299730,11.770196,11.292558,68.955920,
	0.299730,11.770196,11.292558,68.955920,
	0.255088,12.376087,11.020081,69.165282,
	0.255088,12.376087,11.020081,69.165282,
};
	
static double cbsta[][3]={
1.52369,1.92448,122.35124, 1.52962,1.91875,122.28332, 1.53149,1.92096,122.53073, 1.53132,1.92149,122.55859,
1.53219,1.91936,122.36077, 1.51371,1.90607,121.58025, 1.53172,1.92135,122.58755, 1.54507,1.92240,122.99168,
1.53180,1.92265,122.48313, 1.53108,1.92040,122.28572, 1.53078,1.91922,122.34940, 1.53148,1.92241,122.84907,
1.52996,1.94084,115.54662, 1.53057,1.92128,122.53531, 1.53085,1.92046,122.42439, 1.52991,1.91734,122.39611,
1.54070,1.91400,122.79225, 1.54500,1.92132,123.02119, 1.53172,1.92002,122.56818, 1.53251,1.91842,122.36359
}; //n c ca cb  len ang tor

static double sglatavg[][3]={
1.5256,1.924,122.6346, 2.067316,2.04,127.254,  2.4798,2.0697,132.897,  3.129,2.087,136.15,    3.4099,2.104,137.9797,
0,1.2617,180.0,        3.1555,2.11659,135.919, 2.32,2.09,130.84,       3.569,2.0937,138.065,  2.6224,2.14,139.848,
2.9603,2.0976,140.648, 2.485,2.10996,133.8356, 1.87548,2.1875,73.1986, 3.107,2.1085,138.898,  4.16488,2.0575,136.45,
1.901,1.87785,119.78,  1.939,1.9906,123.672,   1.95866,1.99,131.016,   3.886,1.9977,143.1736, 3.781,2.0965,138.927 
}; // for sg

static double sglatavg2[][3]={
0.00000, 1.00000,  180.00, 2.35390, 2.14468,  130.58, 2.41571, 2.13243,  135.58, 3.15619, 2.13210,  136.99, 2.93226, 2.16722,  141.46,
0.00000, 1.00000,    0.00, 2.80663, 2.20671,  139.13, 2.51831, 2.13046,  131.80, 3.63956, 2.13159,  138.74, 2.84465, 2.18538,  144.33,
3.11849, 2.14252,  140.96, 2.44354, 2.17746,  137.26, 2.23972, 2.20653,   57.35, 3.14531, 2.15497,  138.77, 4.01554, 2.10163,  135.50,
2.05520, 1.86421,  116.85, 2.07510, 2.02899,  124.14, 2.13220, 2.02264,  134.07, 3.12109, 2.06139,  143.86, 3.11752, 2.15751,  142.69,
}; //exclude cb

static int sgatomnum[]={
1,2,4,5,7,
0,6,4,5,4,
4,4,3,5,7,
2,3,3,10,8,
0,0,0,0,0,0
};

static char aad1[]= {
'A','C','D','E','F',
'G','H','I','K','L',
'M','N','P','Q','R',
'S','T','V','W','Y',
'J','B','Z','X','O','U'
};

static const char aad3[][4] = {
"ALA","CYS","ASP","GLU","PHE",
"GLY","HIS","ILE","LYS","LEU",
"MET","ASN","PRO","GLN","ARG",
"SER","THR","VAL","TRP","TYR",
"XLE","ASX","GLX","UNK","PYL","SEC"
};

static double asarea2[]= {
240.270, 267.707, 283.510, 315.039, 363.235, 
215.365, 333.600, 312.879, 347.779, 320.895, 
323.209, 288.210, 269.461, 317.071, 379.887,
249.458, 275.129, 286.960, 404.584, 376.261, 
300,300,300,300,300,300
};

static double scalaa[]={
	0.99,1.01,1.01,1.00,1.01,
	1.02,0.99,1.02,0.95,1.02,
	0.99,1.00,0.98,0.99,0.96,
	1.00,0.99,1.01,0.98,0.99,
};//peraa with 4.0

static const double scaa[]={//ef.solvediff()
0.85,0.78,1.02,1.11,1.00,
0.79,1.09,0.99,1.19,1.02,
1.02,1.03,1.02,1.10,1.19,
0.90,1.00,1.00,1.02,1.11,
};

static double avgnormsol[][41]={
//0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19     0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19   avg
0.178,0.100,0.264,0.333,0.153,0.173,0.254,0.145,0.359,0.153,0.228,0.259,0.237,0.302,0.314,0.221,0.218,0.157,0.169,0.185,0.787,0.443,1.166,1.470,0.676,0.764,1.121,0.643,1.585,0.678,1.009,1.146,1.045,1.334,1.387,0.976,0.963,0.692,0.746,0.819,0.226,// aassC natss
0.108,0.052,0.233,0.269,0.093,0.086,0.188,0.081,0.299,0.090,0.110,0.202,0.207,0.240,0.263,0.160,0.147,0.084,0.112,0.129,0.665,0.319,1.436,1.657,0.570,0.532,1.160,0.500,1.844,0.552,0.675,1.247,1.276,1.477,1.621,0.987,0.906,0.516,0.693,0.793,0.162,// aassH
0.047,0.038,0.140,0.197,0.062,0.042,0.128,0.054,0.236,0.057,0.076,0.134,0.150,0.178,0.210,0.097,0.116,0.056,0.092,0.096,0.482,0.390,1.429,2.015,0.631,0.426,1.310,0.558,2.419,0.584,0.782,1.372,1.537,1.818,2.148,0.991,1.189,0.573,0.945,0.981,0.098,// aassE
0.122,0.069,0.240,0.282,0.103,0.140,0.207,0.088,0.314,0.101,0.144,0.227,0.223,0.254,0.273,0.181,0.172,0.093,0.125,0.138,0.692,0.390,1.361,1.603,0.588,0.795,1.175,0.503,1.784,0.572,0.818,1.290,1.266,1.441,1.549,1.030,0.978,0.528,0.712,0.784,0.176,// aa

0.180,0.092,0.270,0.334,0.147,0.168,0.250,0.147,0.358,0.156,0.206,0.260,0.226,0.313,0.320,0.216,0.211,0.158,0.170,0.185,0.799,0.407,1.196,1.478,0.650,0.743,1.107,0.649,1.587,0.691,0.913,1.154,1.003,1.385,1.418,0.957,0.934,0.701,0.753,0.821,0.226,//aass0 predss
0.105,0.050,0.229,0.265,0.079,0.083,0.179,0.075,0.292,0.086,0.092,0.207,0.185,0.235,0.256,0.147,0.138,0.079,0.100,0.109,0.677,0.324,1.470,1.701,0.506,0.535,1.151,0.482,1.876,0.551,0.593,1.331,1.192,1.509,1.643,0.942,0.886,0.510,0.644,0.700,0.156,//aass1
0.054,0.042,0.132,0.184,0.056,0.052,0.115,0.054,0.208,0.057,0.063,0.117,0.123,0.155,0.181,0.093,0.111,0.057,0.082,0.082,0.605,0.472,1.469,2.048,0.626,0.576,1.280,0.605,2.314,0.634,0.698,1.300,1.367,1.732,2.017,1.038,1.240,0.630,0.915,0.917,0.090,//aass2
0.123,0.066,0.244,0.279,0.094,0.143,0.201,0.084,0.307,0.098,0.124,0.231,0.212,0.253,0.267,0.176,0.166,0.090,0.118,0.127,0.711,0.382,1.415,1.617,0.543,0.827,1.167,0.489,1.781,0.568,0.721,1.337,1.229,1.466,1.547,1.022,0.961,0.521,0.685,0.734,0.173,//aa
};

// for sidecent.comm(ITASSER)
static double aalf[][3] = {
0.000,	0.000,	0.000, 
0.249,	-1.118, 0.976, 
0.169,	-1.369,	1.103, 
-0.073,	-1.201,	1.476, 
0.274, 	-1.162,	1.480, 
0.090, 	-1.296,	1.346, 
0.100, 	-1.363,	1.769, 
-0.743,	-1.563,	0.438, 
-0.049,	-1.246,	2.308, 
-0.221,	-1.249,	1.769, 
-0.357,	-1.096,	1.849, 
-0.057,	-1.161,	2.128, 
0.027,	-1.616,	2.597, 
-0.013,	-1.554,	2.219, 
-0.086,	-1.439,	2.296, 
0.113, 	-1.932,	2.933, 
-0.221,	-1.138,	2.165, 
0.111, 	-0.984,	2.447, 
0.128, 	-1.035,	2.604, 
0.476, 	-1.156,	2.541 
};

static double abet[][3] = {
0.000, 	0.000, 	0.000, 
0.113, 	-0.736, 	1.294, 
0.227, 	-0.966, 	1.427, 
0.084, 	-0.738, 	1.712, 
0.093, 	-0.583, 	1.799, 
0.070, 	-0.854, 	1.633, 
-0.105,	-0.601, 	2.135, 
-0.980,	-1.183, 	0.976, 
0.094, 	-0.723, 	2.610, 
0.334, 	-0.664, 	1.992, 
0.097, 	-0.699, 	1.962, 
0.003, 	-0.393, 	2.400, 
-0.019,	-0.745, 	2.972, 
0.101, 	-0.793, 	2.684, 
0.041, 	-0.707, 	2.666, 
-0.020,	-0.998, 	3.394, 
-0.133,	-0.598, 	2.363, 
-0.363,	-0.632, 	2.507, 
-0.375,	-0.601, 	2.706, 
-0.058,	-0.427, 	2.894	
};

static double balf[][3] = {
0.000, 	0.000,		0.000, 
0.249, 	-1.118, 	0.976, 
0.192, 	-1.089, 	1.016, 
0.189, 	-1.073, 	1.025, 
0.271, 	-1.095, 	1.019, 
0.218, 	-1.089, 	1.029, 
0.276, 	-1.100, 	1.015, 
0.154, 	-1.091, 	0.971, 
0.235, 	-1.098, 	1.002, 
0.140, 	-1.049, 	1.050, 
0.047, 	-1.034, 	1.064, 
0.244, 	-1.107, 	0.996, 
0.215, 	-1.105, 	0.994, 
0.239, 	-1.110, 	0.990, 
0.216, 	-1.102, 	0.999, 
0.222, 	-1.102, 	0.996, 
0.152, 	-1.070, 	1.031, 
0.210, 	-1.085, 	1.012, 
0.207, 	-1.087, 	1.010, 
0.228, 	-1.095, 	1.006
};

static double bbet[][3] = {
0.000, 	0.000, 		0.000, 
0.113, 	-0.736, 	1.294, 
0.130, 	-0.765, 	1.276, 
0.098, 	-0.726, 	1.308, 
0.090, 	-0.659, 	1.373, 
0.145, 	-0.740, 	1.310, 
0.085, 	-0.654, 	1.376, 
0.096, 	-0.827, 	1.258, 
0.102, 	-0.702, 	1.327, 
0.093, 	-0.749, 	1.287, 
0.073, 	-0.738, 	1.292, 
0.107, 	-0.692, 	1.334, 
0.108, 	-0.718, 	1.311, 
0.109, 	-0.717, 	1.309, 
0.101, 	-0.710, 	1.316, 
0.100, 	-0.719, 	1.310, 
0.082, 	-0.729, 	1.305, 
0.097, 	-0.706, 	1.323, 
0.098, 	-0.709, 	1.322, 
0.108, 	-0.707, 	1.322
};

#endif