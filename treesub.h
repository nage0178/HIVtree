#ifndef TREESUB_H
#define TREESUB_H
double lfun(double x[], int np);
double lfundG(double x[], int np);
int OutTreeN(FILE *fout, int spnames, int printopt);
int ReadTreeN (FILE *ftree, int *haslength, int copyname, int popline);
int ReadTreeSeqs(FILE *fout);
int DeRoot (void);
int printPatterns(FILE *fout);
void printSeqs(FILE *fout, unsigned char *z[], unsigned char *spnames[], int ns, int ls, int npatt, double fpatt[], int *pose, char keep[], int format);
int RemoveIndel(void);
int PatternWeight (void);
void EncodeSeqs (void);
void PrintTree(int timebranches);
void NodeToBranch (void);
int OutTreeB (FILE *fout);
int SetNodeScale(int inode);
int NodeScale(int inode, int pos0, int pos1);
int eigenTN93(int model, double kappa1, double kappa2, double pi[],
	      int *nR, double Root[], double Cijk[]);
int fx_r(double x[], int np);
void ReadLatentSeqs(void);

#include <regex.h>

#endif
