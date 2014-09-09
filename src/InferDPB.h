/** =============================================================================
 *  Description: preprocess for InforDPB.
 *  Author:  Songpeng Zu, Bioinformatics Division, TNLIST, Tsinghua University
 *  Email: zsp07@mails.tsinaghua.edu.cn
 *  Last Update: June 21, 2014
 ===============================================================================**/

#ifndef __INFERDPB_H__
#define __INFERDPB_H__

#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<cmath>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_matrix.h>
#include<time.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<cstring>
#include<vector>
#include<boost/thread/thread.hpp>
#include<boost/bind.hpp>
#include<unistd.h>
#include<ctype.h>

using namespace std;

/*Constant Definition*/
const int NFILE = 10; // Number of files can be input.
const int NAMELEN = 256; // Length of file names.
const int GETIN = 65536; // Maximum of buffer
const int STEP = 300;  // EM loop number.
//const char *SPLIT = "\t"; // symbol to split the row into slices.

/*Structure Definition*/

struct inputpara
{
    double fn;
    double fp;
    int threadNum;
    char filepath[NFILE][NAMELEN];
};

struct matrix_dim
{
    int row;
    int col;
};

struct parameter
{
	int DrugNum;
	int SubNum;
	int ProNum;
	int DomNum;
	int threadNum;
        double fn;
        double fp;
};

/*Function Declaration*/
void read_matrix_int(const char *str,gsl_matrix_int* Matrix);

void read_matrix_double(const char *str,gsl_matrix* Matrix);

void write_matrix_double(const char *str,gsl_matrix* Matrix);

void write_vector_double(const char *str,vector<double> &delta);

void vector_inition(vector<int> empty_v,vector<vector<int> > &vector_use,int row);

void vector_assign_int(vector<vector<int> > &vector_use,gsl_matrix_int* Matrix);

void E_step_Boost(gsl_matrix* O,vector<vector<int> > &Drug_Sub,vector<vector<int> > &Protein_Domain,\
                          gsl_matrix* lamda,parameter* stead,int thread);

void E_step(gsl_matrix* O,vector<vector<int> > &Drug_Sub,vector<vector<int> > &Protein_Domain,\
                         gsl_matrix* lamda,double fn,double fp,int DrugNum,int ProNum,int threadNum);

void M_step_Boost(vector<vector<int> > &Sub_Drug,vector<vector<int> > &Domain_Protein,\
                         gsl_matrix_int* Drug_Protein,gsl_matrix* lamda,gsl_matrix* O,parameter* stead,int thread);

void M_step(int SubNum,int DomNum,vector<vector<int> > &Sub_Drug,vector<vector<int> > &Domain_Protein,\
                   gsl_matrix_int* Drug_Protein,gsl_matrix* lamda,double fn,gsl_matrix* O,int threadNum);

void LoglikehoodEM(double &loglikehood,int DrugNum,int ProNum,vector<vector<int> > &Drug_Sub,\
                          vector<vector<int> > &Protein_Domain,\
                          gsl_matrix_int* Drug_Protein,gsl_matrix* lamda,double fn,double fp);

void EM(vector<vector<int> > &Drug_Sub,gsl_matrix_int* Drug_Protein,vector<vector<int> > &Sub_Drug,\
               vector<vector<int> > &Protein_Domain,vector<vector<int> > &Domain_Protein,\
               gsl_matrix* lamda,gsl_matrix* O,double fn,double fp,double &loglikehood,int &step,\
               vector<double> &delta,int &DrugNum,int &DomNum,int &ProNum,int &SubNum,gsl_matrix* ZD,int threadNum);


#endif
