
/** =============================================================================
 *  Description: preprocess for InforDPB.
 *  Author:  Songpeng Zu, Bioinformatics Division, TNLIST, Tsinghua University
 *  Email: zsp07@mails.tsinaghua.edu.cn
 *  Last Update: June 21, 2014
 ===============================================================================**/


#include "PreInfer.h"

/*Function Declaration*/

void Usage()
{
    string boundary(80,'=');
    cout \
        << "======================Infer Drug Protein Interaction Positions==================="<<endl
        << endl
        << "Usage:                                                                           "<<endl
        << "    <InferDPB> -f 0.7 -p 0.0001 -t 8 -s drug2sub -d pro2dom -i dpinter -o outfile"<<endl
        <<"-f: false negative                                                                "<<endl
        <<"-p: false positive                                                                "<<endl
        <<"-t: multi thread number                                                           "<<endl
        <<"-s: drug2sub matrix file                                                          "<<endl
        <<"-d: protein2domain matrix file                                                    "<<endl
        <<"-i: drug2protein interaction file                                                 "<<endl
        <<"-n: sub2domain inition file                                                       "<<endl
        <<"-o: outputfile                                                                    "<<endl
        <<"-h: help.                                                                         "<<endl
        <<'\n'<<boundary<<endl;
}


int getoption(int iargc,char **iargv, inputpara &para)
{
    int c;
    opterr = 1;
    while((c = getopt(iargc,iargv,"f:p:t:s:d:i:n:o:h") )!= -1){
        switch (c){
            case 'f':
                para.fn = atof(optarg);
                break;
            case 'p':
                para.fp = atof(optarg);
                break;
            case 't':
                para.threadNum = atoi(optarg);
                break;
            case 's':
                strcpy(para.filepath[0],optarg);
                break;
            case 'd':
                strcpy(para.filepath[1],optarg);
                break;
            case 'i':
                strcpy(para.filepath[2],optarg);
                break;
            case 'n':
                strcpy(para.filepath[3],optarg);
                break;
            case 'o':
                strcpy(para.filepath[4],optarg);
                break;
            case 'h':
                Usage();
                exit(0);
                break;
            case '?':
                if(optopt=='c')
                    fprintf(stderr,"Option -%c requires an argument.\n",optopt);
                else if(isprint(optopt))
                    fprintf(stderr,"Unknown option '-%c' .\n",optopt);
                else
                    fprintf(stderr,"Unknown option character '\\x%x' .\n",optopt);
                return 1;
            default:
                Usage();
                exit(1);
        }
    }
    if(optind < iargc){
        printf("Non-option argument:");
        for(int index=optind;index<iargc;index++){
            printf(" %s",iargv[index]);
        }
        printf("\nType -h to obtain usage.\n");
        exit(1);
    }
    
    return 0;
}

int fileglance(const char *filepath,const char *sep,bool rh,bool ch, matrix_dim &p)
{
	ifstream file(filepath);
	char buf[GETIN];
	if (!file){
		cerr << "error: unable to open the data file: No such file or directory" << endl;
		return 1;
	}
	file.getline(buf, GETIN);
	string instr = buf;
	string::size_type pos1, pos2;
	pos1 = 0;
	pos2 = instr.find(sep);
	int col = 1;
	while(string::npos != pos2){
		pos1 = pos2+1;
		pos2 = instr.find(sep, pos1);
		++col;
	}
	if (ch){
		p.col = col;
	}else{
		p.col = col - 1;
	}
	int row_num = 1;
	while(!file.eof()){
		file.getline(buf, GETIN);
		if (buf[0] != '\0'){
			++row_num;
		}
	}
	file.clear();
	file.close();
	if (rh){
		p.row = row_num;
	}else{
		p.row = row_num - 1;
	}
	return 0;
}


