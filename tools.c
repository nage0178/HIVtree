/*
    Copyright (C) 2022 Anna Nagel
    Originally modified from PAML version 4.9 by Ziheng Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "mcmctree.h"
#include "tools.h"

char BASEs[] = "TCAGUYRMKSWHBVD-N?";
char *EquateBASE[] = { "T","C","A","G", "T", "TC","AG","CA","TG","CG","TA",
     "TCA","TCG","CAG","TAG", "TCAG","TCAG","TCAG" };
char AAs[] = "ARNDCQEGHILKMFPSTWYV-*?X";
char BINs[] = "TC";
char CODONs[256][4];
char nChara[256], CharaMap[256][64];

int noisy = 0, NFunCall = 0;
static time_t time_start;
double PjumpOptimum = 0.30; /* this is the optimum for the Bactrian move. */
unsigned int z_rndu = 666, w_rndu = 1237;

int splitline(char line[], int nfields, int fields[])
{
/* This finds out how many fields there are in the line, and marks the starting positions of the fields.
   if(nfield>0), only nfield fiends are read.  Otherwise read until end of line or until MAXNFIELDS.
   Fields are separated by spaces, and texts are allowed as well.
   returns the number of fields read.
*/
   int i, nfieldsread = 0, InSpace = 1;
   char *p = line;

   for (i = 0; (nfields==-1 || nfieldsread<nfields) && *p && *p != '\n'; i++, p++) {
      if (isspace(*p))
         InSpace = 1;
      else {
         if (InSpace) {
            InSpace = 0;
            fields[nfieldsread ++] = i;
            if (nfieldsread > MAXNFIELDS)
               puts("raise MAXNFIELDS?");
         }
      }
   }
   return(nfieldsread);
}


int testXMat(double x[])
{
   /* test whether X matrix is acceptable (0) or not (-1) */
   int it = 0, i, j;
   double t;
   for (i = 0, t = 0; i < 4; i++) FOR(j, 4) {
      if (x[i * 4 + j] < 0 || x[i * 4 + j]>1)  it = -1;
      t += x[i * 4 + j];
   }
   if (fabs(t - 1) > 1e-4) it = -1;
   return(it);
}


double distance(double x[], double y[], int n)
{
   int i; double t = 0;
   for (i = 0; i < n; i++) t += square(x[i] - y[i]);
   return(sqrt(t));
}

double sum(double x[], int n)
{
   int i; double t = 0;  for (i = 0; i < n; i++) t += x[i];    return(t);
}

int abyx(double a, double x[], int n)
{
   int i; for (i = 0; i < n; x[i] *= a, i++) {}  return(0);
}

int fillxc(double x[], double c, int n)
{
   int i; for (i = 0; i < n; i++) x[i] = c; return (0);
}

int printsma(FILE*fout, char*spname[], unsigned char*z[], int ns, int l, int lline, int gap, int seqtype,
   int transformed, int simple, int pose[])
{
   /* print multiple aligned sequences.
      use spname==NULL if no seq names available.
      pose[h] marks the position of the h_th site in z[], useful for
      printing out the original sequences after site patterns are collapsed.
      Sequences z[] are coded if(transformed) and not if otherwise.
   */
   int igroup, ngroup, lt, h, hp, i, b, b0 = -1, igap, lspname = 30, lseqlen = 7;
   char indel = '-', ambi = '?', equal = '.';
   char *pch = (seqtype <= 1 ? BASEs : (seqtype == 2 ? AAs : BINs));
   char codon[4] = "   ";

   if (l == 0) return(1);
   codon[0] = -1;  /* to avoid warning */
   if (gap == 0) gap = lline + 1;
   ngroup = (l - 1) / lline + 1;
   fprintf(fout, "\n");
   for (igroup = 0; igroup < ngroup; igroup++) {
      lt = min2(l, (igroup + 1)*lline);  /* seqlen mark at the end of block */
      igap = lline + (lline / gap) + lspname + 1 - lseqlen - 1; /* spaces */
      if (igroup + 1 == ngroup)
         igap = (l - igroup*lline) + (l - igroup*lline) / gap + lspname + 1 - lseqlen - 1;
      /* fprintf (fout,"%*s[%*d]\n", igap, "", lseqlen,lt); */
      for (i = 0; i < ns; i++) {
         if (spname) fprintf(fout, "%-*s  ", lspname, spname[i]);
         for (h = igroup*lline, lt = 0, igap = 0; lt < lline && h < l; h++, lt++) {
            hp = (pose ? pose[h] : h);
            if (seqtype == CODONseq && transformed) {
               fprintf(fout, " %s", CODONs[(int)z[i][hp]]);
               continue;
            }
            b0 = (int)z[0][hp];
            b = (int)z[i][hp];
            if (transformed) {
               b0 = pch[b0];
               b = pch[b];
            }
            if (i&&simple && b == b0 && b != indel && b != ambi)
               b = equal;
            fputc(b, fout);
            if (++igap == gap) {
               fputc(' ', fout); igap = 0;
            }
         }
         fprintf(fout, "\n");
      }
      fprintf(fout, "\n");
   }
   fprintf(fout, "\n");
   return(0);
}

int PopEmptyLines(FILE* fseq, int lline, char line[])
{
   /* pop out empty lines in the sequence data file.
      returns -1 if EOF.
   */
   char *eqdel = ".-?", *p;
   int i;

   for (i = 0; ; i++) {
      p = fgets(line, lline, fseq);
      if (p == NULL) return(-1);
      while (*p)
         if (*p == eqdel[0] || *p == eqdel[1] || *p == eqdel[2] || isalpha(*p))
            /*
                     if (*p==eqdel[0] || *p==eqdel[1] || *p==eqdel[2] || isalnum(*p))
            */
            return(0);
         else p++;
   }
}

int matIout(FILE *fout, int x[], int n, int m)
{
   int i, j;
   fprintf(fout, "\n");
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++) fprintf(fout, "  %4d", x[i*m + j]);
      fprintf(fout, "\n");
   }
   return (0);
}


int blankline(char *str)
{
   char *p = str;
   while (*p) if (isalnum(*p++)) return(0);
   return(1);
}

void strcase(char *str, int direction)
{
   /* direction = 0: to lower; 1: to upper */
   char *p = str;
   if (direction)  while (*p) { *p = (char)toupper(*p); p++; }
   else           while (*p) { *p = (char)tolower(*p); p++; }
}


int ScanFastaFile(FILE *fin, int *ns, int *ls, int *aligned)
{
   /* This scans a fasta alignment file to get com.ns & com.ls.
      Returns -1 if the sequences are not aligned and have different lengths.
   */
   int len = 0, ch, starter = '>', stop = '/';  /* both EOF and / mark the end of the file. */
   char name[200], *p;

   if (noisy) printf("\nprocessing fasta file");
   for (*aligned = 1, *ns = -1, *ls = 0; ; ) {
      ch = fgetc(fin);
      if (ch == starter || ch == EOF || ch == stop) {
         if (*ns >= 0) {  /* process end of the sequence */
            if (noisy) printf(" %7d sites", len);

            if (*ns > 1 && len != *ls) {
               *aligned = 0;
               printf("previous sequence %s has len %d, current seq has %d\n", name, *ls, len);
            }
            if (len > *ls) *ls = len;
         }
         (*ns)++;      /* next sequence */
         if (ch == EOF || ch == stop) break;
         /* fscanf(fin, "%s", name); */
         p = name;
         while ((ch = getc(fin)) != '\n' && ch != EOF) *p++ = ch;
         *p = '\0';
         if (noisy) printf("\nreading seq#%2d %-50s", *ns + 1, name);
         len = 0;
      }
      else if (isgraph(ch)) {
         if (*ns == -1)
            error2("seq file error: use '>' in fasta format.");
         len++;
      }
   }
   rewind(fin);
   return(0);
}


double QuantileNormal(double prob)
{
   /* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
      returns (-9999) if in error
      Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
      Applied Statistics 22: 96-97 (AS70)

      Newer methods:
        Wichura MJ (1988) Algorithm AS 241: the percentage points of the
          normal distribution.  37: 477-484.
        Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
          points of the normal distribution.  26: 118-121.
   */
   double a0 = -.322232431088, a1 = -1, a2 = -.342242088547, a3 = -.0204231210245;
   double a4 = -.453642210148e-4, b0 = .0993484626060, b1 = .588581570495;
   double b2 = .531103462366, b3 = .103537752850, b4 = .0038560700634;
   double y, z = 0, p = prob, p1;

   p1 = (p < 0.5 ? p : 1 - p);
   if (p1 < 1e-20) z = 999;
   else {
      y = sqrt(log(1 / (p1*p1)));
      z = y + ((((y*a4 + a3)*y + a2)*y + a1)*y + a0) / ((((y*b4 + b3)*y + b2)*y + b1)*y + b0);
   }
   return (p < 0.5 ? -z : z);
}

double IncompleteGamma(double x, double alpha, double ln_gamma_alpha)
{
   /* returns the incomplete gamma ratio I(x,alpha) where x is the upper
              limit of the integration and alpha is the shape parameter.
      returns (-1) if in error
      ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
      (1) series expansion,     if (alpha>x || x<=1)
      (2) continued fraction,   otherwise
      RATNEST FORTRAN by
      Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
      19: 285-287 (AS32)
   */
   int i;
   double p = alpha, g = ln_gamma_alpha;
   double accurate = 1e-10, overflow = 1e60;
   double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0, pn[6];

   if (x == 0) return (0);
   if (x < 0 || p <= 0) return (-1);

   factor = exp(p*log(x) - x - g);
   if (x > 1 && x >= p) goto l30;
   /* (1) series expansion */
   gin = 1;  term = 1;  rn = p;
l20:
   rn++;
   term *= x / rn;   gin += term;
   if (term > accurate) goto l20;
   gin *= factor / p;
   goto l50;
l30:
   /* (2) continued fraction */
   a = 1 - p;   b = a + x + 1;  term = 0;
   pn[0] = 1;  pn[1] = x;  pn[2] = x + 1;  pn[3] = x*b;
   gin = pn[2] / pn[3];
l32:
   a++;
   b += 2;
   term++;
   an = a*term;
   for (i = 0; i < 2; i++)
      pn[i + 4] = b*pn[i + 2] - an*pn[i];
   if (pn[5] == 0) goto l35;
   rn = pn[4] / pn[5];
   dif = fabs(gin - rn);
   if (dif > accurate) goto l34;
   if (dif <= accurate*rn) goto l42;
l34:
   gin = rn;
l35:
   for (i = 0; i < 4; i++) pn[i] = pn[i + 2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i = 0; i < 4; i++) pn[i] /= overflow;
   goto l32;
l42:
   gin = 1 - factor*gin;

l50:
   return (gin);
}

double factorial(int n)
{
   double fact = 1, i;
   if (n > 100) printf("factorial(%d) may be too large\n", n);
   for (i = 2; i <= n; i++) fact *= i;
   return (fact);
}

double CDFNormal(double x)
{
   /* Hill ID (1973) The normal integral. Applied Statistics, 22:424-427.  (Algorithm AS 66)
      Adams AG (1969) Algorithm 39. Areas under the normal curve. Computer J. 12: 197-198.
      adapted by Z. Yang, March 1994.
   */
   int invers = 0;
   double p, t = 1.28, y = x*x / 2;

   if (x < 0) { invers = 1;  x = -x; }
   if (x < t)
      p = .5 - x * (.398942280444 - .399903438504 * y
         / (y + 5.75885480458 - 29.8213557808
            / (y + 2.62433121679 + 48.6959930692
               / (y + 5.92885724438))));
   else {
      p = 0.398942280385 * exp(-y) /
         (x - 3.8052e-8 + 1.00000615302 /
         (x + 3.98064794e-4 + 1.98615381364 /
            (x - 0.151679116635 + 5.29330324926 /
            (x + 4.8385912808 - 15.1508972451 /
               (x + 0.742380924027 + 30.789933034 /
               (x + 3.99019417011))))));
   }
   return (invers ? p : 1 - p);
}


double logCDFNormal(double x)
{
   /* logarithm of CDF of N(0,1).

      The accuracy is good for the full range (-inf, 38) on my 32-bit machine.
      When x=38, log(F(x)) = -2.88542835e-316.  When x > 38, log(F(x)) can't be
      distinguished from 0.  F(5) = 1 - 1.89E-8, and when x>5, F(x) is hard to
      distinguish from 1.  Instead the smaller tail area F(-5) is used for the
      calculation, using the expansion log(1-z) = -z(1 + z/2 + z*z/3), where
      z = F(-5) is small.
      For 3 < x < 7, both approaches are close, but when x = 8, Mathematica and
      log(CDFNormal) give the incorrect answer -6.66133815E-16, while the correct
      value is log(F(8)) = log(1 - F(-8)) ~= -F(-8) = -6.22096057E-16.
      Note on 2019.1.5: In R, pnorm(8,0,1, log.p=TRUE) gives -6.220961e-16, which is correct.

      F(x) when x<-10 is reliably calculatd using the series expansion, even though
      log(CDFNormal) works until F(-38) = 2.88542835E-316.

      Regarding calculation of the logarithm of Pr(a < X < b), note that
      F(-9) - F(-10) = F(10) - F(9), but the left-hand side is better computationally.
   */
   double lnF, z = fabs(x), C, low = -10, high = 5;

   /* calculate the log of the smaller area */
   if (x >= low && x <= high)
      return log(CDFNormal(x));
   if (x > high && x < -low)
      lnF = log(CDFNormal(-z));
   else {
      C = 1 - 1 / (z*z) + 3 / (z*z*z*z) - 15 / (z*z*z*z*z*z) + 105 / (z*z*z*z*z*z*z*z);
      lnF = -z*z / 2 - log(sqrt(2 * Pi)*z) + log(C);
   }
   if (x > 0) {
      z = exp(lnF);
      lnF = -z*(1 + z / 2 + z*z / 3 + z*z*z / 4 + z*z*z*z / 5);
   }
   return(lnF);
}

double PDFt(double x, double loc, double scale, double df, double lnConst)
{
   /* CDF of t distribution with lococation, scale, and degree of freedom
   */
   double z = (x - loc) / scale, lnpdf = lnConst;

   if (lnpdf == 0) {
      lnpdf = LnGamma((df + 1) / 2) - LnGamma(df / 2) - 0.5*log(Pi*df);
   }
   lnpdf -= (df + 1) / 2 * log(1 + z*z / df);
   return exp(lnpdf) / scale;
}

double CDFt(double x, double loc, double scale, double df, double lnbeta)
{
   /* CDF of t distribution with location, scale, and degree of freedom
   */
   double z = (x - loc) / scale, cdf;
   double lnghalf = 0.57236494292470008707;  /* log{G(1/2)} = log{sqrt(Pi)} */

   if (lnbeta == 0) {
      lnbeta = LnGamma(df / 2) + lnghalf - LnGamma((df + 1) / 2);
   }
   cdf = CDFBeta(df / (df + z*z), df / 2, 0.5, lnbeta);

   if (z >= 0) cdf = 1 - cdf / 2;
   else     cdf /= 2;
   return(cdf);
}

double PDFNormal(double x, double mu, double sigma2)
{
   return 1 / sqrt(2 * Pi*sigma2)*exp(-.5 / sigma2*(x - mu)*(x - mu));
}


double QuantileChi2(double prob, double v)
{
   /* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
      returns -1 if in error.   0.000002<prob<0.999998
      RATNEST FORTRAN by
          Best DJ & Roberts DE (1975) The percentage points of the
          Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
      Converted into C by Ziheng Yang, Oct. 1993.
   */
   double e = .5e-6, aa = .6931471805, p = prob, g, smallv = 1e-6;
   double xx, c, ch, a = 0, q = 0, p1 = 0, p2 = 0, t = 0, x = 0, b = 0, s1, s2, s3, s4, s5, s6;

   if (p < smallv)     return(0);
   if (p > 1 - smallv) return(9999);
   if (v <= 0)         return (-1);

   g = LnGamma(v / 2);
   xx = v / 2;   c = xx - 1;
   if (v >= -1.24*log(p)) goto l1;

   ch = pow((p*xx*exp(g + xx*aa)), 1 / xx);
   if (ch - e < 0) return (ch);
   goto l4;
l1:
   if (v > .32) goto l3;
   ch = 0.4;   a = log(1 - p);
l2:
   q = ch;  p1 = 1 + ch*(4.67 + ch);  p2 = ch*(6.73 + ch*(6.66 + ch));
   t = -0.5 + (4.67 + 2 * ch) / p1 - (6.73 + ch*(13.32 + 3 * ch)) / p2;
   ch -= (1 - exp(a + g + .5*ch + c*aa)*p2 / p1) / t;
   if (fabs(q / ch - 1) - .01 <= 0) goto l4;
   else                       goto l2;

l3:
   x = QuantileNormal(p);
   p1 = 0.222222 / v;
   ch = v*pow((x*sqrt(p1) + 1 - p1), 3.0);
   if (ch > 2.2*v + 6)
      ch = -2 * (log(1 - p) - c*log(.5*ch) + g);
l4:
   q = ch;   p1 = .5*ch;
   if ((t = IncompleteGamma(p1, xx, g)) < 0)
      error2("\nIncompleteGamma");
   p2 = p - t;
   t = p2*exp(xx*aa + g + p1 - c*log(ch));
   b = t / ch;  a = 0.5*t - b*c;

   s1 = (210 + a*(140 + a*(105 + a*(84 + a*(70 + 60 * a))))) / 420;
   s2 = (420 + a*(735 + a*(966 + a*(1141 + 1278 * a)))) / 2520;
   s3 = (210 + a*(462 + a*(707 + 932 * a))) / 2520;
   s4 = (252 + a*(672 + 1182 * a) + c*(294 + a*(889 + 1740 * a))) / 5040;
   s5 = (84 + 264 * a + c*(175 + 606 * a)) / 2520;
   s6 = (120 + c*(346 + 127 * c)) / 5040;
   ch += t*(1 + 0.5*t*s1 - b*c*(s1 - b*(s2 - b*(s3 - b*(s4 - b*(s5 - b*s6))))));
   if (fabs(q / ch - 1) > e) goto l4;

   return (ch);
}

double rndgamma(double a)
{
   /* This returns a random variable from gamma(a, 1).
      Marsaglia and Tsang (2000) A Simple Method for generating gamma variables",
      ACM Transactions on Mathematical Software, 26 (3): 363-372.
      This is not entirely safe and is noted to produce zero when a is small (0.001).
    */
   double a0 = a, c, d, u, v, x, smallv = 1E-300;

   if (a < 1) a++;

   d = a - 1.0 / 3.0;
   c = (1.0 / 3.0) / sqrt(d);

   for (; ; ) {
      do {
         x = rndNormal();
         v = 1.0 + c * x;
      } while (v <= 0);

      v *= v * v;
      u = rndu();

      if (u < 1 - 0.0331 * x * x * x * x)
         break;
      if (log(u) < 0.5 * x * x + d * (1 - v + log(v)))
         break;
   }
   v *= d;

   if (a0 < 1)    /* this may cause underflow if a is small, like 0.01 */
      v *= pow(rndu(), 1 / a0);
   if (v == 0)   /* underflow */
      v = smallv;
   return v;
}

double LnGamma(double x)
{
   /* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
      Stirling's formula is used for the central polynomial part of the procedure.

      Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
      Communications of the Association for Computing Machinery, 9:684
   */
   double f = 0, fneg = 0, z, lng;
   int nx = (int)x;

   if ((double)nx == x && nx >= 0 && nx <= 11)
      lng = log(factorial(nx - 1));
   else {
      if (x <= 0) {
         printf("LnGamma(%.6f) not implemented", x);
         if ((int)x - x == 0) { puts("lnGamma undefined"); return(-1); }
         for (fneg = 1; x < 0; x++) fneg /= x;
         if (fneg < 0)
            error2("strange!! check lngamma");
         fneg = log(fneg);
      }
      if (x < 7) {
         f = 1;
         z = x - 1;
         while (++z < 7)
            f *= z;
         x = z;
         f = -log(f);
      }
      z = 1 / (x*x);
      lng = fneg + f + (x - 0.5)*log(x) - x + .918938533204673
         + (((-.000595238095238*z + .000793650793651)*z - .002777777777778)*z + .083333333333333) / x;
   }
   return  lng;
}

int matout2(FILE * fout, double x[], int n, int m, int wid, int deci)
{
   int i, j;
   fprintf(fout, "\n");
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++)
         fprintf(fout, " %*.*g", wid - 1, deci, x[i*m + j]);
      fprintf(fout, "\n");
   }
   return (0);
}

double CDFBeta(double x, double pin, double qin, double lnbeta)
{
   /* Returns distribution function of the standard form of the beta distribution,
      that is, the incomplete beta ratio I_x(p,q).

      This is also known as the incomplete beta function ratio I_x(p, q)

      lnbeta is log of the complete beta function; provide it if known,
      and otherwise use 0.

      This is called from QuantileBeta() in a root-finding loop.

       This routine is a translation into C of a Fortran subroutine
       by W. Fullerton of Los Alamos Scientific Laboratory.
       Bosten and Battiste (1974).
       Remark on Algorithm 179, CACM 17, p153, (1974).
   */
   double ans, c, finsum, p, ps, p1, q, term, xb, xi, y, smallv = 1e-15;
   int n, i, ib;
   static double eps = 0, alneps = 0, sml = 0, alnsml = 0;

   if (x < smallv)        return 0;
   else if (x > 1 - smallv) return 1;
   if (pin <= 0 || qin <= 0) {
      printf("p=%.4f q=%.4f: parameter outside range in CDFBeta", pin, qin);
      return (-1);
   }

   if (eps == 0) {/* initialize machine constants ONCE */
      eps = pow((double)FLT_RADIX, -(double)DBL_MANT_DIG);
      alneps = log(eps);
      sml = DBL_MIN;
      alnsml = log(sml);
   }
   y = x;  p = pin;  q = qin;

   /* swap tails if x is greater than the mean */
   if (p / (p + q) < x) {
      y = 1 - y;
      p = qin;
      q = pin;
   }

   if (lnbeta == 0) lnbeta = LnBeta(p, q);

   if ((p + q) * y / (p + 1) < eps) {  /* tail approximation */
      ans = 0;
      xb = p * log(max2(y, sml)) - log(p) - lnbeta;
      if (xb > alnsml && y != 0)
         ans = exp(xb);
      if (y != x || p != pin)
         ans = 1 - ans;
   }
   else {
      /* evaluate the infinite sum first.  term will equal */
      /* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */
      ps = q - floor(q);
      if (ps == 0)
         ps = 1;

      xb = LnGamma(ps) + LnGamma(p) - LnGamma(ps + p);
      xb = p * log(y) - xb - log(p);

      ans = 0;
      if (xb >= alnsml) {
         ans = exp(xb);
         term = ans * p;
         if (ps != 1) {
            n = (int)max2(alneps / log(y), 4.0);
            for (i = 1; i <= n; i++) {
               xi = i;
               term = term * (xi - ps) * y / xi;
               ans = ans + term / (p + xi);
            }
         }
      }

      /* evaluate the finite sum. */
      if (q > 1) {
         xb = p * log(y) + q * log(1 - y) - lnbeta - log(q);
         ib = (int)(xb / alnsml);  if (ib < 0) ib = 0;
         term = exp(xb - ib * alnsml);
         c = 1 / (1 - y);
         p1 = q * c / (p + q - 1);

         finsum = 0;
         n = (int)q;
         if (q == (double)n)
            n = n - 1;
         for (i = 1; i <= n; i++) {
            if (p1 <= 1 && term / eps <= finsum)
               break;
            xi = i;
            term = (q - xi + 1) * c * term / (p + q - xi);
            if (term > 1) {
               ib = ib - 1;
               term = term * sml;
            }
            if (ib == 0)
               finsum = finsum + term;
         }
         ans = ans + finsum;
      }
      if (y != x || p != pin)
         ans = 1 - ans;
      if (ans > 1) ans = 1;
      if (ans < 0) ans = 0;
   }
   return ans;
}

double logPDFSkewN(double x, double loc, double scale, double shape)
{
   double z = (x - loc) / scale, lnpdf = 2 / scale;

   lnpdf = 0.5*log(2 / (Pi*scale*scale)) - z*z / 2 + logCDFNormal(shape*z);
   return lnpdf;
}

double PDFSkewT(double x, double loc, double scale, double shape, double df)
{
   double z = (x - loc) / scale, pdf;
   double lnghalf = 0.57236494292470008707;    /* log{G(1/2)} = log{sqrt(Pi)} */
   double lngv, lngv1, lnConst_pdft, lnbeta_cdft;

   lngv = LnGamma(df / 2);
   lngv1 = LnGamma((df + 1) / 2);
   lnConst_pdft = lngv1 - lngv - 0.5*log(Pi*df);
   lnbeta_cdft = lngv1 + lnghalf - lngv - log(df / 2);  /* log{ B((df+1)/2, 1/2) }  */

   pdf = 2 / scale * PDFt(z, 0, 1, df, lnConst_pdft)
      * CDFt(shape*z*sqrt((df + 1) / (df + z*z)), 0, 1, df + 1, lnbeta_cdft);

   return pdf;
}

double PDFSkewN(double x, double loc, double scale, double shape)
{
   double z = (x - loc) / scale, pdf = 2 / scale;

   pdf *= PDFNormal(z, 0, 1) * CDFNormal(shape*z);
   return pdf;
}

double Binomial(double n, int k, double *scale)
{
   /* calculates (n choose k), where n is any real number, and k is integer.
      If(*scale!=0) the result should be c+exp(*scale).
   */
   double c = 1, i, large = 1e200;

   *scale = 0;
   if ((int)k != k)
      error2("k is not a whole number in Binomial.");
   if (k == 0) return(1);
   if (n > 0 && (k<0 || k>n)) return (0);

   if (n > 0 && (int)n == n) k = min2(k, (int)n - k);
   for (i = 1; i <= k; i++) {
      c *= (n - k + i) / i;
      if (c > large) {
         *scale += log(c); c = 1;
      }
   }
   return(c);
}

double reflect(double x, double a, double b)
{
/* This returns a variable in the range (a, b) by reflecting x back into the range
*/
   int side = 0;  /* n is number of jumps over interval.  side=0 (left) or 1 (right). */
   double n, e = 0, smallv = 1e-200;    /* e is excess */

   if (b - a < smallv) {
      printf("\nimproper range x0 = %.9g (%.9g, %.9g)\n", x, a, b);
      exit(-1);
   }
   if (x < a) { e = a - x;  side = 0; }
   else if (x > b) { e = x - b;  side = 1; }
   if (e) {
      n = floor(e / (b - a));
      if (fmod(n, 2.0) > 0.1)   /* fmod should be 0 if n is even and 1 if n is odd. */
         side = 1 - side;       /* change side if n is odd */
      e -= n*(b - a);
      x = (side ? b - e : a + e);
   }

   /* If x lands on boundary after reflection, sample a point at random in the interval. */
   smallv = (b - a)*1e-9;
   while (x - a < smallv || b - x < smallv)  
      x = a + (b - a)*rndu();

   return(x);
}

double logPriorRatioGamma(double xnew, double x, double a, double b)
{
   /* This calculates the log of prior ratio when x has a gamma prior G(x; a, b) with mean a/b
      and x is updated from xold to xnew.
   */
   return (a - 1)*log(xnew / x) - b*(xnew - x);
}

int scanfile(FILE*fin, int *nrecords, int *nx, int *HasHeader, char line[], int ifields[])
{
   /* If the first line has letters, it is considered to be the header line, and HasHeader=0 is set.
   */
   int  i, lline = 1000000, nxline = 0, eof = 0, hastext;

   *nx = 0;  *HasHeader = 0;
   for (*nrecords = 0; ; ) {
      if (!fgets(line, lline, fin)) break;
      eof = feof(fin);
      if (*nrecords == 0 && strchr(line, '\n') == NULL)
         puts(" line too short or too long?");
      for (i = 0, hastext = 0; i < lline && line[i]; i++)
         if (line[i] != 'e' && line[i] != 'E' && isalpha(line[i])) { hastext = 1; break; }
      if (hastext) {
         if (*nrecords == 0) {
            *HasHeader = 1;
            printf("\nData file has a header line.\n");
         }
         else {
            printf("text found on line %d.", *nrecords + 1);
            error2("file format");
         }
      }
      nxline = splitline(line, MAXNFIELDS, ifields);

      if (nxline == 0)
         continue;
      if (*nrecords == 0)
         *nx = nxline;
      else if (*nx != nxline) {
         if (eof)
            break;
         else {
            printf("file format error: %d fields in line %d while %d fields in first line.",
               nxline, *nrecords + 1, *nx);
            error2("error in scanfile()");
         }
      }
      if (*nx > MAXNFIELDS) error2("raise MAXNFIELDS?");

      (*nrecords)++;
      /* printf("line # %3d:  %3d variables\n", *nrecords+1, nxline); */
   }
   rewind(fin);

   if (*HasHeader) {
      fgets(line, lline, fin);
      splitline(line, MAXNFIELDS, ifields);
   }
   if (*HasHeader)
      (*nrecords)--;

   return(0);
}

double Eff_IntegratedCorrelationTime(double x[], int n, double *mx, double *vx, double *rho1)
{
   /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
      sequence method.
      Note that this destroys x[].
   */
   double Tint = 1, rho0 = 0, rho, m = 0, s = 0;
   int  i, ir, minNr = 10, maxNr = 2000;

   /* if(n<1000) puts("chain too short for calculating Eff? "); */
   for (i = 0; i < n; i++) m += x[i];
   m /= n;
   for (i = 0; i < n; i++) x[i] -= m;
   for (i = 0; i < n; i++) s += x[i] * x[i];
   s = sqrt(s / n);
   for (i = 0; i < n; i++) x[i] /= s;

   if (mx) { *mx = m; *vx = s*s; }
   if (s / (fabs(m) + 1) < 1E-9)
      Tint = n;
   else {
      for (ir = 1; ir < min2(maxNr, n - minNr); ir++) {
         for (i = 0, rho = 0; i < n - ir; i++)
            rho += x[i] * x[i + ir];
         rho /= (n - ir);
         if (ir == 1) *rho1 = rho;
         if (ir > minNr && rho + rho0 < 0) break;
         Tint += rho * 2;
         rho0 = rho;
      }
   }
   return (1 / Tint);
}

int HPDinterval(double x[], int n, double HPD[2], double alpha)
{
   /* This calculates the HPD interval at the alpha level.
   */
   int jL0 = (int)(n*alpha / 2), jU0 = (int)(n*(1 - alpha / 2)), jL, jU, jLb = jL0;
   double w0 = x[jU0] - x[jL0], w = w0;
   int debug = 0;

   HPD[0] = x[jL0];
   HPD[1] = x[jU0];
   if (n < 3) return(-1);
   for (jL = 0, jU = jL + (jU0 - jL0); jU < n; jL++, jU++) {
      if (x[jU] - x[jL] < w) {
         jLb = jL;
         w = x[jU] - x[jL];
      }
   }
   HPD[0] = x[jLb];
   HPD[1] = x[jLb + jU0 - jL0];
   return(0);
}

char *printtime(char timestr[])
{
   /* print time elapsed since last call to starttimer()
   */
   time_t t;
   int h, m, s;

   t = time(NULL) - time_start;
   h = (int)t / 3600;
   m = (int)(t % 3600) / 60;
   s = (int)(t - (t / 60) * 60);
   if (h) sprintf(timestr, "%d:%02d:%02d", h, m, s);
   else   sprintf(timestr, "%2d:%02d", m, s);
   return(timestr);
}

int zero(double x[], int n)
{
   int i; for (i = 0; i < n; i++) x[i] = 0; return (0);
}

int ResetStepLengths(FILE *fout, double Pjump[], double finetune[], int nsteps)
{
   /* this abjusts the MCMC proposal step lengths, using equation 9 in
      Yang, Z. & Rodríguez, C. E. 2013 Searching for efficient Markov chain Monte Carlo proposal kernels. Proc. Natl .Acad. Sci. U.S.A. 110, 1930719312.
      PjumpOptimum = 0.3 is also from that paper.
   */
   int j, verybadstep = 0;
   double maxstep = 99;  /* max step length */

   if (noisy >= 3) {
      printf("\n(nsteps = %d)\nCurrent Pjump:    ", nsteps);
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", Pjump[j]);
      printf("\nCurrent finetune: ");
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", finetune[j]);
   }
   if (fout) {
      fprintf(fout, "\nCurrent Pjump:    ");
      for (j = 0; j < nsteps; j++)
         fprintf(fout, " %8.5f", Pjump[j]);
      fprintf(fout, "\nCurrent finetune: ");
      for (j = 0; j < nsteps; j++)
         fprintf(fout, " %8.5f", finetune[j]);
   }

   for (j = 0; j < nsteps; j++) {
      if (Pjump[j] < 0.001) {
         finetune[j] /= 100;
         verybadstep = 1;
      }
      else if (Pjump[j] > 0.999) {
         finetune[j] = min2(maxstep, finetune[j] * 100);
         verybadstep = 1;
      }
      else {
         finetune[j] *= tan(Pi / 2 * Pjump[j]) / tan(Pi / 2 * PjumpOptimum);
         finetune[j] = min2(maxstep, finetune[j]);
      }
   }

   if (noisy >= 3) {
      printf("\nNew     finetune: ");
      for (j = 0; j < nsteps; j++)
         printf(" %8.5f", finetune[j]);
      printf("\n\n");
   }
   if (fout) {
      fprintf(fout, "\nNew     finetune: ");
      for (j = 0; j < nsteps; j++)
         fprintf(fout, " %8.5f", finetune[j]);
      fprintf(fout, "\n");
   }

   return(verybadstep);
}

void starttimer(void)
{
   time_start = time(NULL);
}

int xtoy(double x[], double y[], int n)
{
   int i; for (i = 0; i < n; y[i] = x[i], i++) {}  return(0);
}

int DiscreteGamma(double freqK[], double rK[], double alpha, double beta, int K, int UseMedian)
{
   /* discretization of G(alpha, beta) with equal proportions in each category.
   */
   int i;
   double t, mean = alpha / beta, lnga1;

   if (UseMedian) {   /* median */
      for (i = 0; i < K; i++) rK[i] = QuantileGamma((i*2. + 1) / (2.*K), alpha, beta);
      for (i = 0, t = 0; i < K; i++) t += rK[i];
      for (i = 0; i < K; i++) rK[i] *= mean*K / t;   /* rescale so that the mean is alpha/beta. */
   }
   else {            /* mean */
      lnga1 = LnGamma(alpha + 1);
      for (i = 0; i < K - 1; i++) /* cutting points, Eq. 9 */
         freqK[i] = QuantileGamma((i + 1.0) / K, alpha, beta);
      for (i = 0; i < K - 1; i++) /* Eq. 10 */
         freqK[i] = IncompleteGamma(freqK[i] * beta, alpha + 1, lnga1);
      rK[0] = freqK[0] * mean*K;
      for (i = 1; i < K - 1; i++)  rK[i] = (freqK[i] - freqK[i - 1])*mean*K;
      rK[K - 1] = (1 - freqK[K - 2])*mean*K;
   }

   for (i = 0; i < K; i++) freqK[i] = 1.0 / K;

   return (0);
}

int PMatK80(double P[], double t, double kappa)
{
   /* PMat for JC69 and K80
   */
   int i, j;
   double e1, e2;

   if (t < -1e-6)
      printf("\nt = %.5f in PMatK80", t);

   e1 = expm1(-4 * t / (kappa + 2));
   if (fabs(kappa - 1) < 1e-20) {
      for (i = 0; i < 4; i++)
         for (j = 0; j < 4; j++)
            if (i == j) P[i * 4 + j] = 1. + 3 / 4.0 * e1;
            else        P[i * 4 + j] = -e1 / 4;
   }
   else {
      e2 = expm1(-2 * t*(kappa + 1) / (kappa + 2));
      for (i = 0; i < 4; i++)
         P[i * 4 + i] = 1 + (e1 + 2 * e2) / 4;
      P[0 * 4 + 1] = P[1 * 4 + 0] = P[2 * 4 + 3] = P[3 * 4 + 2] = (e1 - 2 * e2) / 4;
      P[0 * 4 + 2] = P[0 * 4 + 3] = P[2 * 4 + 0] = P[3 * 4 + 0] =
         P[1 * 4 + 2] = P[1 * 4 + 3] = P[2 * 4 + 1] = P[3 * 4 + 1] = -e1 / 4;
   }
   return (0);
}

int PMatTN93(double P[], double a1t, double a2t, double bt, double pi[])
{
   double T = pi[0], C = pi[1], A = pi[2], G = pi[3], Y = T + C, R = A + G;
   double e1, e2, e3, smallv = -1e-6;

   if (noisy && (a1t < smallv || a2t < smallv || bt < smallv))
      printf("\nat=%12.6f %12.6f  bt=%12.6f", a1t, a2t, bt);

   e1 = expm1(-bt);
   e2 = expm1(-(R*a2t + Y*bt));
   e3 = expm1(-(Y*a1t + R*bt));

   P[0 * 4 + 0] = 1 + (R*T*e1 + C*e3) / Y;
   P[0 * 4 + 1] = (R*e1 - e3)*C / Y;
   P[0 * 4 + 2] = -A*e1;
   P[0 * 4 + 3] = -G*e1;

   P[1 * 4 + 0] = (R*e1 - e3)*T / Y;
   P[1 * 4 + 1] = 1 + (R*C*e1 + T*e3) / Y;
   P[1 * 4 + 2] = -A*e1;
   P[1 * 4 + 3] = -G*e1;

   P[2 * 4 + 0] = -T*e1;
   P[2 * 4 + 1] = -C*e1;
   P[2 * 4 + 2] = 1 + Y*A / R*e1 + G / R*e2;
   P[2 * 4 + 3] = Y*G / R*e1 - G / R*e2;

   P[3 * 4 + 0] = -T*e1;
   P[3 * 4 + 1] = -C*e1;
   P[3 * 4 + 2] = Y*A / R*e1 - A / R*e2;
   P[3 * 4 + 3] = 1 + Y*G / R*e1 + A / R*e2;

   return(0);
}

int matout(FILE *fout, double x[], int n, int m)
{
   int i, j;
   fprintf(fout, "\n");
   for (i = 0; i < n; i++) {
      for (j = 0; j < m; j++)
         fprintf(fout, " %11.6f", x[i*m + j]);
      fprintf(fout, "\n");
   }
   return (0);
}

FILE *gfopen(char *filename, char *mode)
{
   FILE *fp;

   if (filename == NULL || filename[0] == 0)
      error2("file name empty.");

   fp = (FILE*)fopen(filename, mode);
   if (fp == NULL) {
      printf("\nerror when opening file %s\n", filename);
      if (!strchr(mode, 'r')) exit(-1);
      printf("tell me the full path-name of the file? ");
      scanf("%s", filename);
      if ((fp = (FILE*)fopen(filename, mode)) != NULL)  return(fp);
      puts("Can't find the file.  I give up.");
      exit(-1);
   }
   return(fp);
}

void error2(char * message)
{
   fprintf(stderr, "\nError: %s.\n", message); 
   exit(-1);
}

void SetSeed(int seed, int PrintSeed)
{
   /* Note seed is of type int with -1 meaning "please find a seed".
     z_rndu and w_rndu are of type unsigned int.
   */
   if (sizeof(int) != 4)
      error2("oh-oh, we are in trouble.  int not 32-bit?  rndu() assumes 32-bit int.");

   if (seed <= 0) {
      FILE *frand = fopen("/dev/urandom", "r");
      if (frand) {
         if (fread(&seed, sizeof(int), 1, frand) != 1)
            error2("failure to read white noise...");
         fclose(frand);
         seed = abs(seed * 2 - 1);
      }
      else {
         seed = abs(1234 * (int)time(NULL) + 1);
      }

      if (PrintSeed) {
         FILE *fseed;
         fseed = fopen("SeedUsed", "w");
         if (fseed == NULL) error2("can't open file SeedUsed.");
         fprintf(fseed, "%d\n", seed);
         fclose(fseed);
      }
   }

   z_rndu = (unsigned int)seed;
   w_rndu = (unsigned int)seed;
}

double rnduM0V1(void)
{
   /* uniform with mean 0 and variance 1 */
   return  1.732050807568877*(-1 + rndu() * 2);
}

double rndNormal(void)
{
   /* Standard normal variate, using the Box-Muller method (1958), improved by
      Marsaglia and Bray (1964).  The method generates a pair of N(0,1) variates,
      but only one is used.
      Johnson et al. (1994), Continuous univariate distributions, vol 1. p.153.
   */
   double u, v, s;

   for (; ;) {
      u = 2 * rndu() - 1;
      v = 2 * rndu() - 1;
      s = u*u + v*v;
      if (s > 0 && s < 1) break;
   }
   s = sqrt(-2 * log(s) / s);
   return (u*s);  /* (v*s) is the other N(0,1) variate, wasted. */
}

double rndTriangle(void)
{
   double u, z;
   /* Standard Triangle variate, generated using inverse CDF  */
   u = rndu();
   if (u > 0.5)
      z = sqrt(6.0) - 2.0*sqrt(3.0*(1.0 - u));
   else
      z = -sqrt(6.0) + 2.0*sqrt(3.0*u);
   return z;
}

double rndLaplace(void)
{
   /* Standard Laplace variate, generated using inverse CDF  */
   double u, r;
   u = rndu() - 0.5;
   r = log(1 - 2 * fabs(u)) * 0.70710678118654752440;
   return (u >= 0 ? -r : r);
}

double rndBactrian(void)
{
   /* This returns a variate from the 1:1 mixture of two normals N(-m, 1-m^2) and N(m, 1-m^2),
      which has mean 0 and variance 1.

      The value m = 0.95 is useful for generating MCMC proposals
   */
   double z = mBactrian + rndNormal()*sBactrian;
   if (rndu() < 0.5) z = -z;
   return (z);
}


double rndBactrianTriangle(void)
{
   /* This returns a variate from the 1:1 mixture of two Triangle Tri(-m, 1-m^2) and Tri(m, 1-m^2),
      which has mean 0 and variance 1.
   */
   double z = mBactrian + rndTriangle()*sBactrian;
   if (rndu() < 0.5) z = -z;
   return (z);
}

double rndBactrianLaplace(void)
{
   /* This returns a variate from the 1:1 mixture of two Laplace Lap(-m, 1-m^2) and Lap(m, 1-m^2),
      which has mean 0 and variance 1.
   */
   double z = mBactrian + rndLaplace()*sBactrian;
   if (rndu() < 0.5) z = -z;
   return (z);
}

int appendfile(FILE*fout, char*filename)
{
   FILE *fin = fopen(filename, "r");
   int ch, status = 0;

   if (fin == NULL) {
      printf("file %s not found!", filename);
      status = -1;
   }
   else {
      while ((ch = fgetc(fin)) != EOF)
         fputc(ch, fout);
      fclose(fin);
      fflush(fout);
   }
   return(status);
}

double rndu(void)
{
   /* 32-bit integer assumed.  From Ripley (1987) p. 46 or table 2.4 line 2. */
#if 0
   z_rndu = z_rndu * 69069 + 1;
   if (z_rndu == 0 || z_rndu == 4294967295)  z_rndu = 13;
   return z_rndu / 4294967295.0;
#else
   z_rndu = z_rndu * 69069 + 1;
   if (z_rndu == 0)  z_rndu = 12345671;
   return ldexp((double)z_rndu, -32);
#endif
}

int comparedouble(const void *a, const void *b)
{
   double aa = *(double*)a, bb = *(double*)b;
   return (aa > bb ? 1 : (aa < bb ? -1 : 0));
}

