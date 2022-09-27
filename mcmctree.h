#ifndef MCMCTREE_H
#define MCMCTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <search.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <assert.h>



#define NS            500
#define NBRANCH      (NS*2-2)
#define NNODE        (NS*2-1)
#define MAXNSONS      3
#define NGENE         8000          /* used for gnodes[NGENE] */
#define LSPNAME       100
#define NCODE         64
#define NCATG         50
#define MaxNFossils   200


int GetOptions(char *ctlf);
int ProcessNodeAnnotation();
int ReadBlengthGH(char infile[]);
int GenerateBlengthGH(char infile[]);
int GetMem(void);
void FreeMem(void);
int UseLocus(int locus, int copycondP, int setmodel, int setSeqName);
int AcceptRejectLocus(int locus, int accept);
void switchconPin(void);
int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[]);
void getSinvDetS(double space[]);
int GetInitials(void);
int GenerateGtree(int locus);
int printGtree(int printBlength);
int ConditionalPNode(int inode, int igene, double x[]);
double lnpData(double lnpDi[]);
double lnpD_locus(int locus);
double lnpD_locus_Approx(int locus);
double lnptNCgiventC(void);
double lnptC(void);
double lnptCalibrationDensity(double t, int fossil, double p[7]);
int SetupPriorTimesFossilErrors(void);
double lnpriorTimesBDS_Approach1(void);
double lnpriorTimesBDS_Approach1Latent(void);
double lnpriorTimesTipDate_Approach2(void);
double lnpriorTimesTipDateEquation6Stadler2010(void);
double lnpriorTimes(void);
double lnpriorRates(void);
void copySptree(void);
void printStree(void);
double InfinitesitesClock(double *FixedDs);
double Infinitesites(FILE *fout);
int collectx(FILE* fout, double x[]);
int MCMC(FILE* fout);
int LabelOldCondP(int spnode);
int UpdateTimes(double *lnL, double steplength[], char accept[]);
int UpdateTimesClock23(double *lnL, double steplength[], char accept[]);
int UpdateParaRates(double *lnL, double steplength[], char accept[], double space[]);
int UpdateRates(double *lnL, double steplength[], char accept[]);
int UpdateParameters(double *lnL, double steplength[], char accept[]);
int mixing(double *lnL, double steplength, char *accept);
int mixingTipDate(double *lnL, double steplength, char *accept);
int mixingCladeStretch(double *lnL, double steplength, char *accept);
int UpdatePFossilErrors(double steplength, char *accept);
int getPfossilerr(double postEFossil[], double nround);
int DescriptiveStatisticsSimpleMCMCTREE(FILE *fout, char infile[]);
double lnpriorTimesBDS_Approach1_OnlyRoot(void);



#endif
