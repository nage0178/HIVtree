#ifndef TOOLS_H
#define TOOLS_H
#include <stdio.h>
int testXMat (double x[]);
double distance (double x[], double y[], int n);
int printsma (FILE*fout, char*spname[], unsigned char*z[], int ns, int l, int lline, int gap, int seqtype,
    int transformed, int simple, int pose[]);
int fillxc (double x[], double c, int n);
double sum (double x[], int n);
int abyx (double a, double x[], int n);
int matIout (FILE *fout, int x[], int n, int m);
int PopEmptyLines (FILE* fseq, int lline, char line[]);
int blankline (char *str);
void strcase (char *str, int direction);
int ScanFastaFile(FILE *f, int *ns, int *ls, int *aligned);
double QuantileChi2 (double prob, double v);
double QuantileNormal (double prob);
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
double logCDFNormal (double x);
double PDFt (double x, double loc, double scale, double df, double lnConst);
double CDFt (double x, double loc, double scale, double df, double lnbeta);
double PDFNormal (double x, double mu, double sigma2);
double CDFNormal (double x);
int splitline(char line[], int nfields, int fields[]);
double factorial (int n);
double rndgamma(double alpha);
double LnGamma(double alpha);
int matout2 (FILE *fout, double x[], int n, int m, int wid, int deci);
double CDFBeta(double x, double p, double q, double lnbeta);
double logPDFSkewN(double x, double loc, double scale, double shape);
double PDFSkewT (double x, double loc, double scale, double shape, double df);
double PDFSkewN (double x, double loc, double scale, double shape);
double Binomial(double n, int k, double *scale);
double reflect(double x, double a, double b);
double logPriorRatioGamma(double xnew, double xold, double a, double b);
int scanfile (FILE*fin, int *nrecords, int *nx, int *HasHeader, char line[], int ifields[]);
double Eff_IntegratedCorrelationTime (double x[], int n, double *mx, double *vx, double *rho1);
int HPDinterval(double x[], int n, double HPD[2], double alpha);
char* printtime(char timestr[]);
int zero (double x[], int n);
int ResetStepLengths(FILE *fout, double Pjump[], double finetune[], int nsteps);

#define mBactrian  0.95
#define sBactrian  sqrt(1 - mBactrian*mBactrian)
#define LnBeta(p,q) (LnGamma(p) + LnGamma(q) - LnGamma(p+q))
#define DGammaUseMedian 0
#define MAXNFIELDS 320000
#define FOR(i,n) for(i=0; i<n; i++)
#define min2(a,b) ((a)<(b)?(a):(b))
#define FPN(file) fputc('\n', file)
#define max2(a,b) ((a)>(b)?(a):(b))
#define square(a) ((a)*(a))
#define Pi  3.1415926535897932384626433832795
#define F0 stdout
#define pamlVerStr "paml version 4.9j, February 2020"
#define QuantileGamma(prob,alpha,beta) QuantileChi2(prob,2.0*(alpha))/(2.0*(beta))

typedef enum {PrBranch=1, PrNodeNum=2, PrLabel=4, PrNodeStr=8, PrAge=16, PrOmega=32} outTreeOptions;
typedef enum {BASEseq=0, CODONseq, AAseq, CODON2AAseq, BINARYseq, BASE5seq} seqTypes;

void starttimer(void);
int xtoy (double x[], double y[], int n);
int UseLocus (int locus, int copyconP, int setmodel, int setSeqName);
int DiscreteGamma (double freqK[], double rK[], double alpha, double beta, int K, int UseMedian);
int PMatK80 (double P[], double t, double kapa);
int PMatTN93 (double P[], double a1t, double a2t, double bt, double pi[]);
int matout (FILE *file, double x[], int n, int m);



FILE *gfopen(char *filename, char *mode);
void error2(char * message);
void SetSeed (int seed, int PrintSeed);
double rnduM0V1 (void);
double rndNormal(void);
double rndTriangle(void);
double rndLaplace (void);
double rndBactrian(void);
double rndBactrianTriangle(void);
double rndBactrianLaplace(void);
int appendfile(FILE*fout, char*filename);
double rndu (void);
int comparedouble (const void *a, const void *b);

#endif
