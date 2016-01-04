/** =============================================================================
 *  Description: preprocess for InforDPB.
 *  Author:  Songpeng Zu, Bioinformatics Division, TNLIST, Tsinghua University
 *  Email: zsp07@mails.tsinaghua.edu.cn
 *  Last Update: June 21, 2014
 ===============================================================================**/

#include "PreInfer.h"


int main(int argc,char *argv[])
{
    inputpara para;
    matrix_dim d2s,d2p,p2d;

    bool rown = true;
    bool coln = true;
    
    const char *split = "\t";

    int flag = getoption(argc,argv,para);
    if(flag==1){
        cerr<< "function getoption runs wrong!" <<endl;
    }

    flag = fileglance(para.filepath[0],split,rown,coln,d2s);
    if(flag==1){
        cerr <<"Error in read <<para.filepath[0]"<<endl;
    }
    
    flag = fileglance(para.filepath[1],split,rown,coln,p2d);
    if(flag==1){
        cerr <<"Error in read <<para.filepath[1]"<<endl;
    }
    
    flag = fileglance(para.filepath[2],split,rown,coln,d2p);
    if(flag==1){
        cerr <<"Error in read <<para.filepath[2]"<<endl;
    }
    

    int DrugNum = d2s.row;
    int SubNum = d2s.col;
    int ProNum = p2d.row;
    int DomNum = p2d.col;
    cout<<" DrugNum "<<DrugNum<<" SubNum "<<SubNum<<" ProNum "<<ProNum<<" DomNum "<<DomNum<<endl;

    int threadNum = para.threadNum;
    double fn = para.fn;
    double fp = para.fp;
    int step = STEP;

    double loglikehood;

    vector<int> empty_v;

    vector<double> delta;

    vector<vector<int> > Drug_Sub;
    vector<vector<int> > Sub_Drug;
    vector<vector<int> > Protein_Domain;
    vector<vector<int> > Domain_Protein;

    vector_inition(empty_v,Drug_Sub,DrugNum);
    vector_inition(empty_v,Sub_Drug,SubNum);
    vector_inition(empty_v,Protein_Domain,ProNum);
    vector_inition(empty_v,Domain_Protein,DomNum);

    gsl_matrix_int* Drug_Protein = gsl_matrix_int_alloc(DrugNum,ProNum);
    gsl_matrix_int* Drug_Sub_tmp = gsl_matrix_int_alloc(DrugNum,SubNum);
    gsl_matrix_int* Sub_Drug_tmp = gsl_matrix_int_alloc(SubNum,DrugNum);
    gsl_matrix_int* Protein_Domain_tmp = gsl_matrix_int_alloc(ProNum,DomNum);
    gsl_matrix_int* Domain_Protein_tmp = gsl_matrix_int_alloc(DomNum,ProNum);


    read_matrix_int(para.filepath[0],Drug_Sub_tmp);
    gsl_matrix_int_transpose_memcpy(Sub_Drug_tmp,Drug_Sub_tmp);

    read_matrix_int(para.filepath[1],Protein_Domain_tmp);
    gsl_matrix_int_transpose_memcpy(Domain_Protein_tmp,Protein_Domain_tmp);
    
    read_matrix_int(para.filepath[2],Drug_Protein);

    vector_assign_int(Drug_Sub,Drug_Sub_tmp);
    vector_assign_int(Sub_Drug,Sub_Drug_tmp);
    vector_assign_int(Protein_Domain,Protein_Domain_tmp);
    vector_assign_int(Domain_Protein,Domain_Protein_tmp);

    gsl_matrix* O = gsl_matrix_alloc(DrugNum,ProNum);
    gsl_matrix* ZD = gsl_matrix_alloc(SubNum,DomNum);
    gsl_matrix* lamda = gsl_matrix_alloc(SubNum,DomNum);

    read_matrix_double(para.filepath[3],lamda);

    // =====================Run EM Algorithm===================//
    EM(Drug_Sub,Drug_Protein,Sub_Drug,Protein_Domain,Domain_Protein,lamda,O,fn,fp,\
       loglikehood,step,delta,DrugNum,DomNum,ProNum,SubNum,ZD,threadNum);

    // =====================Save Result========================//

    write_matrix_double(para.filepath[4],lamda);

    return 0;
}
