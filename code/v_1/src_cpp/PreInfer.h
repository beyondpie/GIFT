/** =============================================================================
 *  Description: preprocess for InforDPB.
 *  Author:  Songpeng Zu, Bioinformatics Division, TNLIST, Tsinghua University
 *  Email: zsp07@mails.tsinaghua.edu.cn
 *  Last Update: June 21, 2014
 ===============================================================================**/

#ifndef PREINFERDPB_H
#define PREINFERDPB_H

#include "InferDPB.h"

/*Function Declaration*/

void Usage();
int getoption(int iargc,char **iargv, inputpara &inputparameter);
int fileglance(const char *filepath,const char *sep,bool rh,bool ch, matrix_dim &dim);

#endif
