/** =============================================================================
 *  Description: preprocess for InforDPB.
 *  Author:  Songpeng Zu, Bioinformatics Division, TNLIST, Tsinghua University
 *  Email: zsp07@mails.tsinaghua.edu.cn
 *  Last Update: June 21, 2014
 ===============================================================================**/


#include "InferDPB.h"

// =============================================================================//
// =======================Function Defination===================================//
// =============================================================================//
void read_matrix_int(const char *array,gsl_matrix_int* Matrix)
{
    FILE *inputfile;
    //const char *array = str.c_str();
    inputfile = fopen(array,"r");
    if(inputfile == NULL)
    {
        cout<<"Error: There is a problem in opening "<<array<<endl;
        exit(1);
    }
    int status = gsl_matrix_int_fscanf(inputfile, Matrix);
    if(status == GSL_EFAILED)
    {
        cout<<"Error: The dimension in "<<array<<" may not equal the dimension inserted."<<endl;
        exit(1);
    }
    fclose(inputfile);
   //cout<<"I found the bug~"<<endl;
}

void read_matrix_double(const char *str,gsl_matrix* Matrix)
{
    FILE *inputfile;
    //const char *array = str.c_str();
    inputfile = fopen(str,"r");
    if(inputfile == NULL)
    {
        cout<<"Error: There is a problem in opening "<<str<<endl;
        exit(1);
    }
    int status = gsl_matrix_fscanf(inputfile, Matrix);
    if(status == GSL_EFAILED)
    {
        cout<<"Error: The dimension in "<<str<<" may not equal the dimension inserted."<<endl;
        exit(1);
    }
    fclose(inputfile);
    //cout<<"I found the bug~"<<endl;
}

void write_matrix_double(const char *str,gsl_matrix* Matrix)
ail{
    FILE *outputfile;
    //const char *array = str.c_str();
    outputfile = fopen(str,"w");
    if(outputfile == NULL)
    {
        cout<<"Error:There is a problem in writing "<<str<<endl;
        exit(1);
    }
    int status = gsl_matrix_fprintf(outputfile,Matrix,"%f");
    if(status == GSL_EFAILED)
    {
        cout<<"Error:There is a problem in writing "<<str<<" maybe the format or matrix is wrong."<<endl;
        exit(1);
    }
    fclose(outputfile);
}

void write_vector_double(const char *str,vector<double> &delta)
{
    FILE *outputfile;
    //const char *array = str.c_str();
    outputfile = fopen(str,"w");
    if(outputfile == NULL)
    {
        cout<<"Error:There is a problem in writing"<<str<<endl;
        exit(1);
    }
    int i;
    for(i=0;i<delta.size();i++)
    {
        fprintf(outputfile,"%f\n",delta[i]);
    }
    fclose(outputfile);
}

void vector_inition(vector<int> empty_v,vector<vector<int> > &vector_use,int row)
{
    int i;
    for(i=0;i<row;i++) {vector_use.push_back(empty_v);}
}

void vector_assign_int(vector<vector<int> > &vector_use,gsl_matrix_int* Matrix)
{
    int row = Matrix->size1;
    int column = Matrix->size2;
    int i,j;
    for(i=0;i<row;i++)
    {
        for(j=0;j<column;j++)
        {
            if(gsl_matrix_int_get(Matrix,i,j)>0)
            {
                vector_use[i].push_back(j);
            }
        }
    }
}

void E_step_Boost(gsl_matrix* O,vector<vector<int> > &Drug_Sub,vector<vector<int> > &Protein_Domain,
            gsl_matrix* lamda,parameter* stead,int thread)
{
    int i,j;
    int m,n;
    double tmp,YP;
    int Zm,Dn;
    for(i=thread; i<stead->DrugNum; i+=stead->threadNum)
    {
        for(j=0;j<stead->ProNum;j++)
        {
            tmp = 0;
            Zm = Drug_Sub[i].size();
            Dn = Protein_Domain[j].size();
            for(m=0;m<Zm;m++)
            {
                for(n=0;n<Dn;n++)
                {
                    tmp +=  log(1-gsl_matrix_get(lamda,Drug_Sub[i][m],Protein_Domain[j][n]));
                }
            }
	    tmp = exp(tmp);
            YP = (1-stead->fn)*(1-tmp)+stead->fp*tmp;
            gsl_matrix_set(O,i,j,YP);
        }
    }
}

void E_step(gsl_matrix* O,vector<vector<int> > &Drug_Sub,vector<vector<int> > &Protein_Domain,\
            gsl_matrix* lamda,double fn,double fp,int DrugNum,int ProNum,int threadNum)
{
    int i;
    boost::thread_group *x;
    x = new boost::thread_group;
    boost::thread *y;
    struct parameter stead;
    stead.fn = fn;
    stead.fp = fp;
    stead.DrugNum = DrugNum;
    stead.ProNum = ProNum; 
    stead.threadNum = threadNum;
    for (i=0; i<threadNum; i++)
    {
        y=new boost::thread(E_step_Boost,O,Drug_Sub,Protein_Domain,lamda,&stead,i);
        x->add_thread(y);
    }
    x->join_all();
    delete x;
}

void M_step_Boost(vector<vector<int> > &Sub_Drug,vector<vector<int> > &Domain_Protein,\
            gsl_matrix_int* Drug_Protein,gsl_matrix* lamda,gsl_matrix* O,parameter* stead,int thread)
{
    int m,n;
    int Y,P;
    int i,j;
    double tmp;
    for(m=thread;m<stead->SubNum;m+=stead->threadNum)
    {
        for(n=0;n<stead->DomNum;n++)
        {
           Y = Sub_Drug[m].size();
           P = Domain_Protein[n].size();
           tmp = 0;
           for(i=0;i<Y;i++)
           {
               for(j=0;j<P;j++)
               {
                    if (gsl_matrix_int_get(Drug_Protein,Sub_Drug[m][i],Domain_Protein[n][j])>0)
                    {
                        tmp = tmp + (1-stead->fn)/gsl_matrix_get(O,Sub_Drug[m][i],Domain_Protein[n][j]);
                    }
                    else
                    {
                        tmp = tmp + stead->fn/(1-gsl_matrix_get(O,Sub_Drug[m][i],Domain_Protein[n][j]));
                    }
               }
           }
           tmp = log(gsl_matrix_get(lamda,m,n)) + log(tmp/(Y*P));
           gsl_matrix_set(lamda,m,n,exp(tmp));
        }
    }
}
void M_step(int SubNum,int DomNum,vector<vector<int> > &Sub_Drug,vector<vector<int> > &Domain_Protein,\
            gsl_matrix_int* Drug_Protein,gsl_matrix* lamda,double fn,gsl_matrix* O,int threadNum)
{
    int m;
    boost::thread_group *x;
    x = new boost::thread_group;
    boost::thread *y;
    struct parameter stead;
    stead.fn = fn;
    stead.SubNum = SubNum;
    stead.DomNum = DomNum;
    stead.threadNum = threadNum;
    for(m=0;m<threadNum;m++)
    {
        y = new boost::thread(M_step_Boost,Sub_Drug,Domain_Protein,Drug_Protein,lamda,O,&stead,m);
        x->add_thread(y);
    }
    x->join_all();
    delete x;
}

void LoglikehoodEM(double &loglikehood,int DrugNum,int ProNum,vector<vector<int> > &Drug_Sub,\
                   vector<vector<int> > &Protein_Domain,\
                   gsl_matrix_int* Drug_Protein,gsl_matrix* lamda,double fn,double fp)
{
    int i,j;
    int m,n;
    int Y,P;
    double tmp,YP;
    loglikehood = 0;
    for(i=0;i<DrugNum;i++)
    {
        for(j=0;j<ProNum;j++)
        {
            tmp = 1;
            Y = Drug_Sub[i].size();
            P = Protein_Domain[j].size();
            for(m=0;m<Y;m++)
            {
                for(n=0;n<P;n++)
                {
                    tmp = tmp * (1-gsl_matrix_get(lamda,Drug_Sub[i][m],Protein_Domain[j][n]));
                }
            }
            YP = (1-fn)*(1-tmp)+fp*tmp;
            if (gsl_matrix_int_get(Drug_Protein,i,j)==1)
                {loglikehood = loglikehood + log(YP);}
            else
                {loglikehood = loglikehood + log(1-YP);}
        }
    }
}

void EM(vector<vector<int> > &Drug_Sub,gsl_matrix_int* Drug_Protein,vector<vector<int> > &Sub_Drug,\
        vector<vector<int> > &Protein_Domain,vector<vector<int> > &Domain_Protein,\
        gsl_matrix* lamda,gsl_matrix* O,double fn,double fp,double &loglikehood,int &step,\
        vector<double> &delta,int &DrugNum,int &DomNum,int &ProNum,int &SubNum,gsl_matrix* ZD,int threadNum)
{
    int st;
    delta.push_back(1);
    clock_t ss,tt;
    double time;
    double max1,max2;
    for(st=1;st<=step;st++)
    {
        gsl_matrix_memcpy(ZD,lamda);

        // E step first

        ss = clock();
        E_step(O,Drug_Sub,Protein_Domain,lamda,fn,fp,DrugNum,ProNum,threadNum);
        tt = clock();
        time = (double)(tt-ss)/CLOCKS_PER_SEC;
        cout<<"E Step time is"<<time<<endl;

        // M step then

        ss = clock();
        M_step(SubNum,DomNum,Sub_Drug,Domain_Protein,Drug_Protein,lamda,fn,O,threadNum);
        tt = clock();
        time = (double)(tt-ss)/CLOCKS_PER_SEC;
        cout<<"M Step time is"<<time<<endl;

        // Judge last

        gsl_matrix_sub(ZD,lamda);
        max1 = gsl_matrix_max(ZD);
        max2 = -1*gsl_matrix_min(ZD);
        if(max1 >= max2)
        {
           delta.push_back(max1);
        }
        else
        {
            delta.push_back(max2);
        }
        cout<<"The Max of DeltaMatrix is "<<delta[st]<<endl;
    }
    LoglikehoodEM(loglikehood,DrugNum,ProNum,Drug_Sub,Protein_Domain,Drug_Protein,lamda,fn,fp);
    cout<<"The Loglikehood is "<<loglikehood<<endl;
}

