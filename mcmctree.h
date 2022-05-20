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



struct CommonInfo {
   unsigned char *z[NS];
   char *spname[NS];
   char seqf[2048], outf[2048], treef[2048], daafile[2048], mcmcf[2048], inBVf[2048], checkpointf[2048], latentf[2048];
   char oldconP[NNODE];       /* update conP for node? (0 yes; 1 no) */
   int seqtype, ns, ls, ngene, posG[2], lgene[3], *pose, npatt, readpattern;
   int np, ncode, ntime, nrate, nrgene, nalpha, npi, ncatG, print, verbose, checkpoint;
   int cleandata, ndata;
   int model, clock, fix_kappa, fix_alpha, fix_rgene, Mgene;
   int method, icode, codonf, aaDist, NSsites;
   int latentData, latentBound;
   double *fpatt, kappa, alpha, TipDate, TipDate_TimeUnit, tipDateOld;
   double latentBoundDate;
   double rgene[NGENE], piG[NGENE][NCODE];  /* not used */
   double(*plfun)(double x[], int np), freqK[NCATG], rK[NCATG], *conP, *fhK;
   double pi[NCODE];
   int curconP;                       /* curconP = 0 or 1 */
   size_t sconP;
   double *conPin[2], space[100000];  /* change space[] to dynamic memory? */
   int conPSiteClass, NnodeScale;
   char *nodeScale;          /* nScale[ns-1] for interior nodes */
   double *nodeScaleF;       /* nScaleF[npatt] for scale factors */
}  com;

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


struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
}  tree;
struct TREEN { /* ipop is node number in species tree */
   int father, nson, sons[2], ibranch, ipop;
   int latent;
   double branch, age, label, label2, *conP, latentAge;
   char fossil, *name, *annotation;
}  *nodes, **gnodes, nodes_t[2 * NS - 1];
/* nodes_t[] is working space.  nodes is a pointer and moves around.
gnodes[] holds the gene trees, subtrees constructed from the master species
tree.  Each locus has a full set of rates (rates) for all branches on the
master tree, held in stree.nodes[].rates.  Branch lengths in the gene
tree are calculated by using those rates and the divergence times.

gnodes[][].label in the gene tree is used to store branch lengths estimated
under no clock when mcmc.usedata=2 (normal approximation to likelihood).
*/


struct SPECIESTREE {
   int nspecies, nbranch, nnode, root, nfossil, duplication, nLatent;
   double RootAge[4];
   struct TREESPN {
      char name[LSPNAME + 1], fossil, usefossil;    /* fossil: 0, 1(L), 2(U), 3(B), 4(G) */
      int father, nson, sons[2], label;
      double age, latentAge, pfossil[7];                       /* parameters in fossil distribution */
      double *rates;
      int latent;                             /* log rates for loci */
   } nodes[2 * NS - 1];
}    stree;
/* all trees are binary & rooted, with ancestors unknown. */


struct DATA { /* locus-specific data and tree information */
   int ns[NGENE], ls[NGENE], npatt[NGENE], ngene, lgene[NGENE];
   int root[NGENE + 1], conP_offset[NGENE];
   int priortime, priorrate;
   char datatype[NGENE], *z[NGENE][NS], cleandata[NGENE];
   double *zmorph[NGENE][NS * 2 - 1], *Rmorph[NGENE];
   double zpopvar[NGENE], ldetRm[NGENE]; /* MdR */
   double *fpatt[NGENE], lnpT, lnpR, lnpDi[NGENE], pi[NGENE][NCODE];
   double kappa[NGENE], alpha[NGENE];
   double BDS[4];  /* parameters in the birth-death-sampling model */
   double kappagamma[2], alphagamma[2];
   double pfossilerror[3], /* (p_beta, q_beta, NminCorrect) */ Pfossilerr, *CcomFossilErr;
   int    rgeneprior;         /* 0: gamma-Dirichlet; 1: conditional iid */
   double rgene[NGENE + 1], sigma2[NGENE + 1], rgenepara[3], sigma2para[3];
   double *blMLE[NGENE], *Gradient[NGENE], *Hessian[NGENE];
   int    transform;
}  data;

struct MCMCPARAMETERS {
   int burnin, nsample, sampfreq, usedata, saveconP, print;
   int nsteplength, steplengthOffset[8];
   char  *accept;
   double *steplength, *Pjump;
}  mcmc; /* control parameters */

enum { JC69, K80, F81, F84, HKY85, T92, TN93, REV } MODELS;
enum { BASE, AA, CODON, MORPHC } DATATYPE;
enum { LOWER_F = 1, UPPER_F, BOUND_F, GAMMA_F, SKEWN_F, SKEWT_F, S2N_F, UNIF_F } FOSSIL_FLAGS;

#endif
