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

/* treesub.c
 subroutines that operates on trees, inserted into other programs
 such as baseml, basemlg, codeml, and pamp.
 */



#include "mcmctree.h"
#include "tools.h"
#include "treesub.h"

#define MCMCTREE 1

extern char BASEs[], *EquateBASE[], AAs[], BINs[], CODONs[][4], nChara[], CharaMap[][64];
extern int NFunCall;
extern char BASEs[];
extern double PjumpOptimum;

extern int noisy;
extern int LASTROUND;
extern double _rateSite;
extern char *fossils[];
extern int npfossils[];

#ifdef  MCMCTREE
#define REALSEQUENCE
#define NODESTRUCTURE
#define LFUNCTIONS
#endif





#ifdef REALSEQUENCE

int hasbase(char *str)
{
   char *p = str, *eqdel = ".-?";
   while (*p)
      if (*p == eqdel[0] || *p == eqdel[1] || *p == eqdel[2] || isalpha(*p++))
         return(1);
   return(0);
}


int GetSeqFileType(FILE *fseq, int *NEXUSseq);
int IdenticalSeqs(void);
void RemoveEmptySequences(void);

int GetSeqFileType(FILE *fseq, int *format)
{
   /* NEXUSstart="begin data" and NEXUSdata="matrix" identify nexus file format.
    Modify if necessary.
    format: 0: alignment; 1: fasta; 2: nexus.
    */
   int  lline = 1000, ch, aligned;
   char fastastarter = '>';
   char line[1000], *NEXUSstart = "begin data", *NEXUSdata = "matrix", *p;
   char *ntax = "ntax", *nchar = "nchar";

   while (isspace(ch = fgetc(fseq)))
      ;
   ungetc(ch, fseq);
   if (ch == fastastarter) {
      *format = 1;
      ScanFastaFile(fseq, &com.ns, &com.ls, &aligned);
      if (aligned)
         return(0);
      else
         error2("The seq file appears to be in fasta format, but not aligned?");
   }
   if (fscanf(fseq, "%d%d", &com.ns, &com.ls) == 2) {
      *format = 0; return(0);
   }
   *format = 2;
   printf("\nseq file is not paml/phylip format.  Trying nexus format.");

   for (; ; ) {
      if (fgets(line, lline, fseq) == NULL) error2("seq err1: EOF");
      strcase(line, 0);
      if (strstr(line, NEXUSstart)) break;
   }
   for (; ; ) {
      if (fgets(line, lline, fseq) == NULL) error2("seq err2: EOF");
      strcase(line, 0);
      if ((p = strstr(line, ntax)) != NULL) {
         while (*p != '=') { if (*p == 0) error2("seq err"); p++; }
         sscanf(p + 1, "%d", &com.ns);
         if ((p = strstr(line, nchar)) == NULL) error2("expect nchar");
         while (*p != '=') { if (*p == 0) error2("expect ="); p++; }
         sscanf(p + 1, "%d", &com.ls);
         break;
      }
   }
   /* printf("\nns: %d\tls: %d\n", com.ns, com.ls);  */
   for (; ; ) {
      if (fgets(line, lline, fseq) == NULL) error2("seq err1: EOF");
      strcase(line, 0);
      if (strstr(line, NEXUSdata)) break;
   }
   return(0);
}

int PopupNEXUSComment(FILE *fseq)
{
   int ch, comment1 = ']';
   for (; ; ) {
      ch = fgetc(fseq);
      if (ch == EOF) error2("expecting ]");
      if (ch == comment1) break;
      if (noisy) putchar(ch);
   }
   return(0);
}


#if(MCMCTREE)

int ReadMorphology(FILE *fout, FILE *fin, int locus)
{
   int i, j;
   char line[1024], str[64];

   if ((data.zmorph[locus][0] = (double*)malloc((com.ns * 2 - 1)*com.ls * sizeof(double))) == NULL)
      error2("oom zmorph");
   if ((data.Rmorph[locus] = (double*)malloc(com.ls*com.ls * sizeof(double))) == NULL)
      error2("oom Rmorph");

   printf("Locus %d has morphological alignment\n", locus + 1);
   for (i = 1; i < com.ns * 2 - 1; i++) {
      data.zmorph[locus][i] = data.zmorph[locus][0] + i*com.ls;
   }
   for (i = 0; i < com.ns; i++) {
      fscanf(fin, "%s", com.spname[i]);
      printf("Reading data for species #%2d: %s     \r", i + 1, com.spname[i]);
      for (j = 0; j < com.ls; j++)
         fscanf(fin, "%lf", &data.zmorph[locus][i][j]);
   }

   for (i = 0; i < com.ns; i++) {
      fprintf(fout, "%-10s ", com.spname[i]);
      for (j = 0; j < com.ls; j++)
         fprintf(fout, " %8.5f", data.zmorph[locus][i][j]);
      FPN(fout);
   }

#if(0)
   fscanf(fin, "%s", str);
   fgets(line, 1024, fin);
   i = j = -1;
   if (strstr("Correlation", str)) {
      for (i = 0; i < com.ls; i++) {
         for (j = 0; j < com.ls; j++)
            if (fscanf(fin, "%lf", &data.Rmorph[locus][i*com.ls + j]) != 1) break;
         if (j < com.ls) break;
      }
   }
   if (i != com.ls || j != com.ls) {
      printf("\ndid not find a good R matrix.  Setting it to identity matrix I.\n");
      for (i = 0; i < com.ls; i++)
         for (j = 0; j < com.ls; j++)
            data.Rmorph[locus][i*com.ls + j] = (i == j);
   }
#endif
   return(0);
}

#endif

int ReadSeq(FILE *fout, FILE *fseq, int cleandata, int locus)
{
   /* read in sequence, translate into protein (CODON2AAseq), and
    This counts ngene but does not initialize lgene[].
    It also codes (transforms) the sequences.
    com.seqtype: 0=nucleotides; 1=codons; 2:AAs; 3:CODON2AAs; 4:BINs
    com.pose[] is used to store gene or site-partition labels.
    ls/3 gene marks for codon sequences.
    char opt_c[]="GIPM";
    G:many genes;  I:interlaved format;  P:patterns;  M:morphological characters

    Use cleandata=1 to clean up ambiguities.  In return, com.cleandata=1 if the
    data are clean or are cleaned, and com.cleandata=0 is the data are unclean.
    */
   char *p, *p1, eq = '.', comment0 = '[', *line;
   int format = 0;  /* 0: paml/phylip, 1: fasta; 2: NEXUS/nexus */
   int i, j, k, ch, noptline = 0, lspname = LSPNAME, miss = 0, nb;
   int lline = 10000, lt[NS], igroup, Sequential = 1, basecoding = 0;
   int n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);
   int gap = (n31 == 3 ? 3 : 10), nchar = (com.seqtype == AAseq ? 20 : 4);
   int h, b[3] = { 0 };
   char *pch = ((com.seqtype <= 1 || com.seqtype == CODON2AAseq) ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   char str[4] = "   ", tmp[32];
   char *NEXUSend = "end;";
   double lst;

#if(MCMCTREE)
   data.datatype[locus] = com.seqtype;
#endif
   str[0] = 0; h = -1; b[0] = -1; /* avoid warning */
   com.readpattern = 0;
   if (com.seqtype == 4) error2("seqtype==BINs, check with author");
   if (noisy >= 9 && (com.seqtype <= CODONseq || com.seqtype == CODON2AAseq)) {
      puts("\n\nAmbiguity character definition table:\n");
      for (i = 0; i < (int)strlen(BASEs); i++) {
         nb = (int)strlen(EquateBASE[i]);
         printf("%c (%d): ", BASEs[i], nb);
         for (j = 0; j < nb; j++)  printf("%c ", EquateBASE[i][j]);
         FPN(F0);
      }
   }
   GetSeqFileType(fseq, &format);

   if (com.ns > NS) error2("too many sequences.. raise NS?");
   if (com.ls%n31 != 0) {
      printf("\n%d nucleotides, not a multiple of 3!", com.ls); exit(-1);
   }
   if (noisy) printf("ns = %d  \tls = %d\n", com.ns, com.ls);

   for (j = 0; j < com.ns; j++) {
      if (com.spname[j]) free(com.spname[j]);
      com.spname[j] = (char*)malloc((lspname + 1) * sizeof(char));
      for (i = 0; i < lspname + 1; i++) com.spname[j][i] = 0;
      if ((com.z[j] = (unsigned char*)realloc(com.z[j], com.ls * sizeof(unsigned char))) == NULL)
         error2("oom z");
   }
   com.rgene[0] = 1;   com.ngene = 1;
   lline = max2(lline, com.ls / n31*(n31 + 1) + lspname + 50);
   if ((line = (char*)malloc(lline * sizeof(char))) == NULL) error2("oom line");

   /* first line */
   if (format == 0) {
      if (!fgets(line, lline, fseq)) error2("ReadSeq: first line");
      com.readpattern = (strchr(line, 'P') || strchr(line, 'p'));
#if(MCMCTREE)
      if (strchr(line, 'M') || strchr(line, 'm')) {
         data.datatype[locus] = MORPHC;
         data.zpopvar[locus] = 0; data.ldetRm[locus] = 0;                /* MdR */
         sscanf(line, "%s %lf %lf", tmp, &data.zpopvar[locus], &data.ldetRm[locus]);  /* MdR */
         printf("MdR: got %f %f, locus: %d\n", data.zpopvar[locus], data.ldetRm[locus], locus);
         /* TODO: Do some error checking. MdR. */
      }
#endif
   }
#if(MCMCTREE)
   if (data.datatype[locus] == MORPHC) { /* morhpological data */
      ReadMorphology(fout, fseq, locus);
      return(0);
   }
   else
#endif
      if (!com.readpattern) {
         if (com.pose) free(com.pose);
         if ((com.pose = (int*)malloc(com.ls / n31 * sizeof(int))) == NULL)
            error2("oom pose");
         for (j = 0; j < com.ls / n31; j++) com.pose[j] = 0;      /* gene #1, default */
      }
      else {
         if (com.pose) free(com.pose);
         com.pose = NULL;
      }
      if (format) goto readseq;

      for (j = 0; j < lline && line[j] && line[j] != '\n'; j++) {
         if (!isalnum(line[j])) continue;
         line[j] = (char)toupper(line[j]);
         switch (line[j]) {
         case 'G': noptline++;   break;
         case 'C': basecoding = 1; break;
         case 'S': Sequential = 1; break;
         case 'I': Sequential = 0; break;
         case 'P':               break;  /* already dealt with. */
         default:
            printf("\nBad option '%c' in first line of seqfile\n", line[j]);
            exit(-1);
         }
      }
      if (strchr(line, 'C')) {   /* protein-coding DNA sequences */
         if (com.seqtype == 2) error2("option C?");
         if (com.seqtype == 0) {
            if (com.ls % 3 != 0 || noptline < 1)  error2("option C?");
            com.ngene = 3;
            for (i = 0; i < 3; i++) com.lgene[i] = com.ls / 3;
         }
         noptline--;
      }

      /* option lines */
      for (j = 0; j < noptline; j++) {
         for (ch = 0; ; ) {
            ch = (char)fgetc(fseq);
            if (ch == comment0)
               PopupNEXUSComment(fseq);
            if (isalnum(ch)) break;
         }

         ch = (char)toupper(ch);
         switch (ch) {
         case ('G'):
            if (basecoding) error2("sequence data file: incorrect option format, use GC?\n");
#if(defined MCMCTREE || defined BPP)
            error2("sequence data file: option G should not be used.");
#endif
            if (fscanf(fseq, "%d", &com.ngene) != 1) error2("expecting #gene here..");
            if (com.ngene > NGENE) error2("raise NGENE?");

            fgets(line, lline, fseq);
            if (!blankline(line)) {    /* #sites in genes on the 2nd line */
               for (i = 0, p = line; i < com.ngene; i++) {
                  while (*p && !isalnum(*p)) p++;
                  if (sscanf(p, "%d", &com.lgene[i]) != 1) break;
                  while (*p && isalnum(*p)) p++;
               }
               /* if ngene is large and some lgene is on the next line */
               for (; i < com.ngene; i++)
                  if (fscanf(fseq, "%d", &com.lgene[i]) != 1) error2("EOF at lgene");

               for (i = 0, k = 0; i < com.ngene; i++)
                  k += com.lgene[i];
               if (k != com.ls / n31) {
                  matIout(F0, com.lgene, 1, com.ngene);
                  printf("\n%6d != %d", com.ls / n31, k);
                  puts("\nOption G: total length over genes is not correct");
                  if (com.seqtype == 1) {
                     puts("Note: gene length is in number of codons.");
                  }
                  puts("Sequence length in number of nucleotides.");
                  exit(-1);
               }
               if (!com.readpattern)
                  for (i = 0, k = 0; i < com.ngene; k += com.lgene[i], i++)
                     for (j = 0; j < com.lgene[i]; j++)
                        com.pose[k + j] = i;

            }
            else {                   /* site marks on later line(s)  */
               if (com.readpattern)
                  error2("option PG: use number of patterns in each gene and not site marks");
               for (k = 0; k < com.ls / n31; ) {
                  if (com.ngene > 9)  fscanf(fseq, "%d", &ch);
                  else {
                     do ch = fgetc(fseq); while (!isdigit(ch));
                     ch = ch - (int)'1' + 1;  /* assumes 1,2,...,9 are consecutive */
                  }
                  if (ch<1 || ch>com.ngene)
                  {
                     printf("\ngene mark %d at %d?\n", ch, k + 1);  exit(-1);
                  }
                  com.pose[k++] = ch - 1;
               }
               if (!fgets(line, lline, fseq)) error2("sequence file, gene marks");
            }
            break;
         default:
            printf("Bad option '%c' in option lines in seqfile\n", line[0]);
            exit(-1);
         }
      }

   readseq:
      /* read sequence */
      if (Sequential) {    /* sequential */
         if (noisy) printf("Reading sequences, sequential format..\n");
         for (j = 0; j < com.ns; j++) {
            lspname = LSPNAME;
            for (i = 0; i < 2 * lspname; i++) line[i] = '\0';
            if (!fgets(line, lline, fseq)) error2("EOF?");
            if (blankline(line)) {
               if (PopEmptyLines(fseq, lline, line))
               {
                  printf("error in sequence data file: empty line (seq %d)\n", j + 1); exit(-1);
               }
            }
            p = line + (line[0] == '=' || line[0] == '>');
            while (isspace(*p)) p++;
            if ((ch = (int)(strstr(p, "  ") - p)) < lspname && ch > 0) lspname = ch;
            strncpy(com.spname[j], p, lspname);
            k = (int)strlen(com.spname[j]);
            p += (k < lspname ? k : lspname);

            for (; k > 0; k--) /* trim spaces */
               if (!isgraph(com.spname[j][k]))   com.spname[j][k] = 0;
               else    break;

               if (noisy >= 2) printf("Reading seq #%2d: %s       %s", j + 1, com.spname[j], (noisy > 3 ? "\n" : "\n"));
               for (k = 0; k < com.ls; p++) {
                  while (*p == '\n' || *p == '\0') {
                     p = fgets(line, lline, fseq);
                     if (p == NULL)
                     {
                        printf("\nEOF at site %d, seq %d\n", k + 1, j + 1); exit(-1);
                     }
                  }
                  *p = (char)toupper(*p);
                  if ((com.seqtype == BASEseq || com.seqtype == CODONseq) && *p == 'U')
                     *p = 'T';
                  p1 = strchr(pch, *p);
                  if (p1 && p1 - pch >= nchar)
                     miss = 1;
                  if (*p == eq) {
                     if (j == 0) error2("Error in sequence data file: . in 1st seq.?");
                     com.z[j][k] = com.z[0][k];  k++;
                  }
                  else if (p1)
                     com.z[j][k++] = *p;
                  else if (isalpha(*p)) {
                     printf("\nError in sequence data file: %c at %d seq %d.\n", *p, k + 1, j + 1);
                     puts("Make sure to separate the sequence from its name by 2 or more spaces.");
                     exit(0);
                  }
                  else if (*p == (char)EOF) error2("EOF?");
               }           /* for(k) */
               if (strchr(p, '\n') == NULL) /* pop up line return */
                  while ((ch = fgetc(fseq)) != '\n' && ch != EOF);
         }   /* for (j,com.ns) */
      }
      else { /* interlaved */
         if (noisy) printf("Reading sequences, interlaved format..\n");
         FOR(j, com.ns) lt[j] = 0;  /* temporary seq length */
         for (igroup = 0; ; igroup++) {
            /*
             printf ("\nreading block %d ", igroup+1);  matIout(F0,lt,1,com.ns);*/

            FOR(j, com.ns) if (lt[j] < com.ls) break;
            if (j == com.ns) break;
            FOR(j, com.ns) {
               if (!fgets(line, lline, fseq)) {
                  printf("\nerr reading site %d, seq %d group %d\nsites read in each seq:",
                     lt[j] + 1, j + 1, igroup + 1);
                  error2("EOF?");
               }
               if (!hasbase(line)) {
                  if (j) {
                     printf("\n%d, seq %d group %d", lt[j] + 1, j + 1, igroup + 1);
                     error2("empty line.");
                  }
                  else
                     if (PopEmptyLines(fseq, lline, line) == -1) {
                        printf("\n%d, seq %d group %d", lt[j] + 1, j + 1, igroup + 1);
                        error2("EOF?");
                     }
               }
               p = line;
               if (igroup == 0) {
                  lspname = LSPNAME;
                  while (isspace(*p)) p++;
                  if ((ch = (int)(strstr(p, "  ") - p)) < lspname && ch > 0)
                     lspname = ch;
                  strncpy(com.spname[j], p, lspname);
                  k = (int)strlen(com.spname[j]);
                  p += (k < lspname ? k : lspname);

                  for (; k > 0; k--)   /* trim spaces */
                     if (!isgraph(com.spname[j][k]))
                        com.spname[j][k] = 0;
                     else
                        break;
                  if (noisy >= 2) printf("Reading seq #%2d: %s     \n", j + 1, com.spname[j]);
               }
               for (; *p && *p != '\n'; p++) {
                  if (lt[j] == com.ls) break;
                  *p = (char)toupper(*p);
                  if ((com.seqtype == BASEseq || com.seqtype == CODONseq) && *p == 'U')
                     *p = 'T';
                  p1 = strchr(pch, *p);
                  if (p1 && p1 - pch >= nchar)
                     miss = 1;
                  if (*p == eq) {
                     if (j == 0) {
                        printf("err: . in 1st seq, group %d.\n", igroup);
                        exit(-1);
                     }
                     com.z[j][lt[j]] = com.z[0][lt[j]];
                     lt[j]++;
                  }
                  else if (p1)
                     com.z[j][lt[j]++] = *p;
                  else if (isalpha(*p)) {
                     printf("\nerr: unrecognised character %c at %d seq %d block %d.",
                        *p, lt[j] + 1, j + 1, igroup + 1);
                     exit(-1);
                  }
                  else if (*p == (char)EOF) error2("EOF");
               }         /* for (*p) */
            }            /* for (j,com.ns) */

            if (noisy > 2) {
               printf("\nblock %3d:", igroup + 1);
               for (j = 0; j < com.ns; j++) printf(" %6d", lt[j]);
            }

         }               /* for (igroup) */
      }
      if (format == 2) {  /* NEXUS format: pop up extra lines until "end;"  */
         for (; ; ) {
            if (fgets(line, lline, fseq) == NULL) error2("seq err1: EOF");
            strcase(line, 0);
            if (strstr(line, NEXUSend)) break;
         }
      }
      free(line);

      /*** delete empty sequences ******************/
#if(0)
      { int ns1 = com.ns, del[100] = { 0 };
      for (i = 0; i < com.ns; i++) {
         for (h = 0; h < com.ls; h++)   if (com.z[i][h] != '?') break;
         if (h == com.ls) {
            del[i] = 1;
            ns1--;
         }
      }
      fprintf(fnew, "\n\n%4d %6d\n", ns1, com.ls);
      for (i = 0; i < com.ns; i++) {
         if (!del[i]) {
            fprintf(fnew, "\n%-40s  ", com.spname[i]);
            for (h = 0; h < com.ls; h++) fprintf(fnew, "%c", com.z[i][h]);
         }
      }
      fflush(fnew);
      }
#endif


      if (!miss)
         com.cleandata = 1;
      else if (cleandata) {  /* forced removal of ambiguity characters */
         if (noisy > 2)  puts("\nSites with gaps or missing data are removed.");
         if (fout) {
            fprintf(fout, "\nBefore deleting alignment gaps\n");
            fprintf(fout, " %6d %6d\n", com.ns, com.ls);
            printsma(fout, com.spname, com.z, com.ns, com.ls, com.ls, gap, com.seqtype, 0, 0, NULL);
         }
         RemoveIndel();
         if (fout) fprintf(fout, "\nAfter deleting gaps. %d sites\n", com.ls);
      }

      if (fout && !com.readpattern) {/* verbose=1, listing sequences again */
         fprintf(fout, " %6d %6d\n", com.ns, com.ls);
         printsma(fout, com.spname, com.z, com.ns, com.ls, com.ls, gap, com.seqtype, 0, 0, NULL);
      }

      if (n31 == 3) com.ls /= n31;

      /* IdenticalSeqs(); */


      if (noisy >= 2) printf("\nSequences read..\n");
      if (com.ls == 0) {
         puts("no sites. Got nothing to do");
         return(1);
      }

#if (defined MCMCTREE)
      /* Check and remove empty sequences.  */
      if (com.cleandata == 1)
         RemoveEmptySequences();

#endif


      /***** GIBBON, 2017.10.8 **/
      /* return(0); */

      if (!com.readpattern) {
         PatternWeight();
      }
      else {  /*  read pattern counts */
         com.npatt = com.ls;
         if ((com.fpatt = (double*)realloc(com.fpatt, com.npatt * sizeof(double))) == NULL)
            error2("oom fpatt");
         for (h = 0, lst = 0; h < com.npatt; h++) {
            fscanf(fseq, "%lf", &com.fpatt[h]);
            lst += com.fpatt[h];
            if (com.fpatt[h]<0 || com.fpatt[h]>1e66)
               printf("fpatth[%d] = %.6g\n", h + 1, com.fpatt[h]);
         }
         if (lst > 1.00001) {
            com.ls = (int)lst;
            if (noisy) printf("\n%d site patterns read, %d sites\n", com.npatt, com.ls);
         }
         if (com.ngene == 1) {
            com.lgene[0] = com.ls;
            com.posG[0] = 0;
            com.posG[1] = com.npatt;
         }
         else {
            for (j = 0, com.posG[0] = 0; j < com.ngene; j++)
               com.posG[j + 1] = com.posG[j] + com.lgene[j];

            for (j = 0; j < com.ngene; j++) {
               com.lgene[j] = (j == 0 ? 0 : com.lgene[j - 1]);
               for (h = com.posG[j]; h < com.posG[j + 1]; h++)
                  com.lgene[j] += (int)com.fpatt[h];
            }
         }
      }

      EncodeSeqs();

      if (fout) {
         fprintf(fout, "\nPrinting out site pattern counts\n\n");
         printPatterns(fout);

         if (com.verbose >= 2) {
            fprintf(fout, "\nSite-to-pattern map: ");
            for (i = 0; i < com.ls; i++)
               fprintf(fout, " %2d", com.pose[i] + 1);
            fprintf(fout, "\n");
         }
      }

      return (0);
}




void RemoveEmptySequences(void)
{
   /* this removes empty sequences (? or - only) and adjust com.ns
    */
   int j, h, nsnew;
   char emptyseq[NS];

   for (j = 0; j < com.ns; j++) {
      emptyseq[j] = 1;
      for (h = 0; h < com.ls*(com.seqtype == 1 ? 3 : 1); h++)
         if (com.z[j][h] != '?' && com.z[j][h] != '-') {
            emptyseq[j] = 0;
            break;
         }
   }
   for (j = 0, nsnew = 0; j < com.ns; j++) {
      if (emptyseq[j]) {
         printf("seq #%3d: %-30s is removed\n", j + 1, com.spname[j]);
         free(com.z[j]);
         free(com.spname[j]);
         continue;
      }
      com.z[nsnew] = com.z[j];
      com.spname[nsnew] = com.spname[j];
      nsnew++;
   }
   for (j = nsnew; j < com.ns; j++) {
      com.z[j] = NULL;
      com.spname[j] = NULL;
   }
   com.ns = nsnew;
   if (com.ns <= 0)
      error2("\nThis locus has nothing left after empty sequences are removed.  Use cleandata = 0.\n");
}


int printPatterns(FILE *fout)
{
   int j, h, n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);
   int gap = (n31 == 3 ? 3 : 10); //, n=(com.seqtype==AAseq?20:4);

   fprintf(fout, "\n%10d %10d  P", com.ns, com.npatt*n31);
   if (com.ngene > 1) {
      fprintf(fout, " G\nG %d  ", com.ngene);
      for (j = 0; j < com.ngene; j++)
         fprintf(fout, "%7d", com.posG[j + 1] - com.posG[j]);
   }
   FPN(fout);

   if (com.seqtype == 1 && com.cleandata) {
      ; /* nothing is printed out for yn00, as the coding is different. */
   }
   else
      printsma(fout, com.spname, com.z, com.ns, com.npatt, com.npatt, gap, com.seqtype, 1, 0, NULL);
   if (com.ls > 1.0001) {
      fprintf(fout, "\n");
      for (h = 0; h < com.npatt; h++) {
         fprintf(fout, " %4.0f", com.fpatt[h]);
         if ((h + 1) % 15 == 0) FPN(fout);
      }
      fprintf(fout, "\n\n");
   }
   return(0);
}



void EncodeSeqs(void)
{
   /* This encodes sequences and set up com.TipMap[][], called after sites are collapsed
      into patterns.
   */
   int n = com.ncode, nA, is, h, i, j, k, ic, indel = 0, ch, b[3];
   char *pch = ((com.seqtype == 0 || com.seqtype == 1) ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   unsigned char c[4] = "", str[4] = "   ";

   if (com.seqtype != 1) {
      for (is = 0; is < com.ns; is++) {
         for (h = 0; h < com.npatt; h++) {
            ch = com.z[is][h];
            com.z[is][h] = (char)(k = (int)(strchr(pch, ch) - pch));
            if (k < 0) {
               printf("strange character %c in seq %d site %d\n", ch, is + 1, h + 1);
               exit(-1);
            }
         }
      }
   }
}


int print1site(FILE*fout, int h)
{
   /* This print out one site in the sequence data, com.z[].  It may be the h-th
   site in the original data file or the h-th pattern.  The data are coded.
   naa > 1 if the codon codes for more than one amino acid.
   */
   char *pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   char compatibleAAs[20] = "";
   int n = com.ncode, i, b, aa = 0;

   for (i = 0; i < com.ns; i++) {
      b = com.z[i][h];
      if (com.seqtype == 0 || com.seqtype == 2)
         fprintf(fout, "%c", pch[b]);
   }
   return(0);
}


void SetMapAmbiguity(int seqtype, int ModelAA2Codon)
{
   /* This sets up CharaMap, the map from the ambiguity characters to resolved characters.
   */
   int n = com.ncode, nc, i, j, i0, i1, i2, nb[3], ib[3][4], iaa, ic;
   char *pch = (seqtype == 0 ? BASEs : (seqtype == 2 ? AAs : BINs));
   char *pbases = (seqtype == 0 ? BASEs : NULL);
   char **pEquateBASE = (seqtype == 0 ? EquateBASE : NULL);
   char debug = 0;

   for (j = 0; j < n; j++) {  /* basic characters, coded according to the definition in pch. */
      nChara[j] = (char)1;
      CharaMap[j][0] = (char)j;
   }
   if (seqtype != 1 && ModelAA2Codon == 0) {  /* nucleotides or amino acids */
      for (j = n, pch += n; *pch; j++, pch++) {
         if (com.seqtype == 0 || com.seqtype == 5) {  /* ambiguities are allowed for those 2 types */
            nChara[j] = (char)strlen(pEquateBASE[j]);
            for (i = 0; i < nChara[j]; i++)
               CharaMap[j][i] = (char)(strchr(pbases, pEquateBASE[j][i]) - pbases);
         }
         else {  /* for non-nucleotide characters, ambiguity characters must be ? or -. */
            nChara[j] = (char)n;
            for (i = 0; i < n; i++)
               CharaMap[j][i] = (char)i;
         }
         if (debug) {
            printf("character %c (%d): ", pbases[j], nChara[j]);
            for (i = 0; i < nChara[j]; i++)
               printf("%c", pbases[(int)CharaMap[j][i]]);
            printf("\n");
         }
      }
   }
}

int PatternWeight(void)
{
   /* This collaps sites into patterns, for nucleotide, amino acid, or codon sequences.
      Site patterns are represented by 0-ended character strings.  The code adds 1 to the character
      during transposing so that it works for both unoded and coded data.
      com.pose[i] has labels for genes as input and maps sites to patterns in return.
      com.fpatt has site-pattern counts.
      This deals with multiple genes/partitions, and uses com.ngene, com.lgene[], com.posG[] etc.
      Sequences z[ns][ls] are copied into patterns zt[ls*lpatt], and bsearch is used
      twice to avoid excessive copying of com.fpatt[] and com.pose[], the first round to count
      npatt and generate the site patterns and the second round to generate fpatt[] & com.pose[].
   */
   int maxnpatt = com.ls, h, ip, l, u, j, k, same, ig, *poset;
   // int gap = (com.seqtype==CODONseq ? 3 : 10);
   int n31 = (com.seqtype == CODONseq ? 3 : 1);
   int lpatt = com.ns*n31 + 1;   /* extra 0 used for easy debugging, can be voided */
   int *p2s;  /* point patterns to sites in zt */
   char timestr[36];
   unsigned char *p, *zt;
   double nc = (com.seqtype == 1 ? 64 : com.ncode) + !com.cleandata + 1;
   int debug = 0;

   /* (A) Collect and sort patterns.  Get com.npatt, com.lgene, com.posG.
   Move sequences com.z[ns][ls] into sites zt[ls*lpatt].
   Use p2s to map patterns to sites in zt to avoid copying.
   */
   if (noisy) printf("Counting site patterns.. %s\n", printtime(timestr));

   if ((com.seqtype == 1 && com.ns < 5) || (com.seqtype != 1 && com.ns < 7))
      maxnpatt = (int)(pow(nc, (double)com.ns) + 0.5) * com.ngene;
   if (maxnpatt > com.ls) maxnpatt = com.ls;
   p2s = (int*)malloc(maxnpatt * sizeof(int));
   zt = (char*)malloc(com.ls*lpatt * sizeof(char));
   if (p2s == NULL || zt == NULL)  error2("oom p2s or zt");
   memset(zt, 0, com.ls*lpatt * sizeof(char));
   for (j = 0; j < com.ns; j++)
      for (h = 0; h < com.ls; h++)
         for (k = 0; k < n31; k++)
            zt[h*lpatt + j*n31 + k] = (unsigned char)(com.z[j][h*n31 + k] + 1);

   for (ig = 0; ig < com.ngene; ig++) com.lgene[ig] = 0;
   for (ig = 0, com.npatt = 0; ig < com.ngene; ig++) {
      com.posG[ig] = l = u = ip = com.npatt;
      for (h = 0; h < com.ls; h++) {
         if (com.pose[h] != ig) continue;
         if (debug) printf("\nh %3d %s", h, zt + h*lpatt);

         /* bsearch in existing patterns.  Knuth 1998 Vol3 Ed2 p.410
         ip is the loc for match or insertion.  [l,u] is the search interval.
         */
         same = 0;
         if (com.lgene[ig]++ != 0) {  /* not 1st pattern? */
            for (l = com.posG[ig], u = com.npatt - 1; ; ) {
               if (u < l) break;
               ip = (l + u) / 2;
               k = strcmp(zt + h*lpatt, zt + p2s[ip] * lpatt);
               if (k < 0)        u = ip - 1;
               else if (k > 0)   l = ip + 1;
               else { same = 1;  break; }
            }
         }
         if (!same) {
            if (com.npatt > maxnpatt)
               error2("npatt > maxnpatt");
            if (l > ip) ip++;        /* last comparison in bsearch had k > 0. */
                                     /* Insert new pattern at ip.  This is the expensive step. */
            if (ip < com.npatt)
               memmove(p2s + ip + 1, p2s + ip, (com.npatt - ip) * sizeof(int));
            p2s[ip] = h;
            com.npatt++;
         }

         if (debug) {
            printf(": %3d (%c ilu %3d%3d%3d) ", com.npatt, (same ? '.' : '+'), ip, l, u);
            for (j = 0; j < com.npatt; j++)
               printf(" %s", zt + p2s[j] * lpatt);
         }
         if (noisy && ((h + 1) % 10000 == 0 || h + 1 == com.ls))
            printf("\rCompressing, %6d patterns at %6d / %6d sites (%.1f%%), %s",
               com.npatt, h + 1, com.ls, (h + 1.) * 100 / com.ls, printtime(timestr));

      }     /* for (h)  */
      if (noisy) FPN(F0);
   }        /* for (ig) */

            /* (B) count pattern frequencies and collect pose[] */
   com.posG[com.ngene] = com.npatt;
   for (j = 0; j < com.ngene; j++)
      if (com.lgene[j] == 0)
         error2("some genes do not have any sites?");
   for (j = 1; j < com.ngene; j++)
      com.lgene[j] += com.lgene[j - 1];

   com.fpatt = (double*)realloc(com.fpatt, com.npatt * sizeof(double));
   poset = (int*)malloc(com.ls * sizeof(int));
   if (com.fpatt == NULL || poset == NULL) error2("oom poset");
   memset(com.fpatt, 0, com.npatt * sizeof(double));

   for (ig = 0; ig < com.ngene; ig++) {
      for (h = 0; h < com.ls; h++) {
         if (com.pose[h] != ig) continue;
         for (same = 0, l = com.posG[ig], u = com.posG[ig + 1] - 1; ; ) {
            if (u < l) break;
            ip = (l + u) / 2;
            k = strcmp(zt + h*lpatt, zt + p2s[ip] * lpatt);
            if (k < 0)        u = ip - 1;
            else if (k > 0)   l = ip + 1;
            else { same = 1;  break; }
         }
         com.fpatt[ip]++;
         poset[h] = ip;
         if (noisy && ((h + 1) % 10000 == 0 || h + 1 == com.ls))
            printf("\rCollecting fpatt[] & pose[], %6d patterns at %6d / %6d sites (%.1f%%), %s",
               com.npatt, h + 1, com.ls, (h + 1.) * 100 / com.ls, printtime(timestr));
      }     /* for (h)  */
      if (noisy) FPN(F0);
   }        /* for (ig) */

   if (com.seqtype == CODONseq && com.ngene == 3 && com.lgene[0] == com.ls / 3)
      puts("\nCheck option G in data file?\n");

   for (j = 0; j < com.ns; j++) {
      com.z[j] = (unsigned char*)realloc(com.z[j], com.npatt*n31 * sizeof(unsigned char));
      for (ip = 0, p = com.z[j]; ip < com.npatt; ip++)
         for (k = 0; k < n31; k++)
            *p++ = (unsigned char)(zt[p2s[ip] * lpatt + j*n31 + k] - 1);
   }
   memcpy(com.pose, poset, com.ls * sizeof(int));
   free(poset);  free(p2s);  free(zt);
   return (0);
}



void AddFreqSeqGene(int js, int ig, double pi0[], double pi[]);


void Chi2FreqHomo(double f[], int ns, int nc, double X2G[2])
{
   /* This calculates a chi-square like statistic for testing that the base
    or amino acid frequencies are identical among sequences.
    f[ns*nc] where ns is #sequences (rows) and nc is #states (columns).
    */
   int i, j;
   double mf[64] = { 0 }, smallv = 1e-50;

   X2G[0] = X2G[1] = 0;
   for (i = 0; i < ns; i++)
      for (j = 0; j < nc; j++)
         mf[j] += f[i*nc + j] / ns;

   for (i = 0; i < ns; i++) {
      for (j = 0; j < nc; j++) {
         if (mf[j] > smallv) {
            X2G[0] += square(f[i*nc + j] - mf[j]) / mf[j];
            if (f[i*nc + j])
               X2G[1] += 2 * f[i*nc + j] * log(f[i*nc + j] / mf[j]);
         }
      }
   }
}

int InitializeBaseAA(FILE *fout)
{
   /* Count site patterns (com.fpatt) and calculate base or amino acid frequencies
    in genes and species.  This works on raw (uncoded) data.
    Ambiguity characters in sequences are resolved by iteration.
    For frequencies in each species, they are resolved within that sequence.
    For average base frequencies among species, they are resolved over all
    species.

    This routine is called by baseml and aaml.  codonml uses another
    routine InitializeCodon()
    */
   char *pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   char indel[] = "-?";
   int wname = 30, h, js, k, ig, nconstp, n = com.ncode;
   int irf, nrf = 20;
   double pi0[20], t, lmax = 0, X2G[2], *pisg;  /* freq for species & gene, for X2 & G */

   if (noisy) printf("Counting frequencies..");
   if (fout)  fprintf(fout, "\nFrequencies..");
   if ((pisg = (double*)malloc(com.ns*n * sizeof(double))) == NULL)
      error2("oom pisg");
   for (h = 0, nconstp = 0; h < com.npatt; h++) {
      for (js = 1; js < com.ns; js++)
         if (com.z[js][h] != com.z[0][h])  break;
      if (js == com.ns && com.z[0][h] != indel[0] && com.z[0][h] != indel[1])
         nconstp += (int)com.fpatt[h];
   }
   for (ig = 0, zero(com.pi, n); ig < com.ngene; ig++) {
      if (com.ngene > 1)
         fprintf(fout, "\n\nGene %2d (len %4d)", ig + 1, com.lgene[ig] - (ig == 0 ? 0 : com.lgene[ig - 1]));
      fprintf(fout, "\n%*s", wname, "");
      for (k = 0; k < n; k++) fprintf(fout, "%7c", pch[k]);

      /* The following block calculates freqs in each species for each gene.
       Ambiguities are resolved in each species.  com.pi and com.piG are
       used for output only, and are not be used later with missing data.
       */
      zero(com.piG[ig], n);
      zero(pisg, com.ns*n);
      for (js = 0; js < com.ns; js++) {
         fillxc(pi0, 1.0 / n, n);
         for (irf = 0; irf < nrf; irf++) {
            zero(com.pi, n);
            AddFreqSeqGene(js, ig, pi0, com.pi);
            t = sum(com.pi, n);
            if (t < 1e-10) {
               printf("Some sequences are empty.\n");
               fillxc(com.pi, 1.0 / n, n);
            }
            else
               abyx(1 / t, com.pi, n);
            if (com.cleandata || com.cleandata || (t = distance(com.pi, pi0, n)) < 1e-8)
               break;
            xtoy(com.pi, pi0, n);
         }   /* for(irf) */
         fprintf(fout, "\n%-*s", wname, com.spname[js]);
         for (k = 0; k < n; k++) fprintf(fout, "%8.5f", com.pi[k]);
         if (com.ncode == 4 && com.ngene == 1) fprintf(fout, " GC = %5.3f", com.pi[1] + com.pi[3]);
         for (k = 0; k < n; k++) com.piG[ig][k] += com.pi[k] / com.ns;
         xtoy(com.pi, pisg + js*n, n);
      }    /* for(js,ns) */
      if (com.ngene > 1) {
         fprintf(fout, "\n\n%-*s", wname, "Mean");
         for (k = 0; k < n; k++) fprintf(fout, "%7.4f", com.piG[ig][k]);
      }

      Chi2FreqHomo(pisg, com.ns, n, X2G);

      fprintf(fout, "\n\nHomogeneity statistic: X2 = %.5f G = %.5f ", X2G[0], X2G[1]);

      /* fprintf(frst1,"\t%.5f", X2G[1]); */

   }  /* for(ig) */
   if (noisy) printf("\n");

   /* If there are missing data, the following block calculates freqs
    in each gene (com.piG[]), as well as com.pi[] for the entire sequence.
    Ambiguities are resolved over entire data sets across species (within
    each gene for com.piG[]).  These are used in ML calculation later.
    */
   if (com.cleandata) {
      for (ig = 0, zero(com.pi, n); ig < com.ngene; ig++) {
         t = (ig == 0 ? com.lgene[0] : com.lgene[ig] - com.lgene[ig - 1]) / (double)com.ls;
         for (k = 0; k < n; k++)  com.pi[k] += com.piG[ig][k] * t;
      }
   }
   else {
      for (ig = 0; ig < com.ngene; ig++) {
         xtoy(com.piG[ig], pi0, n);
         for (irf = 0; irf < nrf; irf++) {  /* com.piG[] */
            zero(com.piG[ig], n);
            for (js = 0; js < com.ns; js++)
               AddFreqSeqGene(js, ig, pi0, com.piG[ig]);
            t = sum(com.piG[ig], n);
            if (t < 1e-10) {
               puts("empty sequences?");
               for (k = 0; k < n; k++)  com.piG[ig][k] = 1.0 / n;
               t = 1;
            }
            abyx(1 / t, com.piG[ig], n);
            if (distance(com.piG[ig], pi0, n) < 1e-8) break;
            xtoy(com.piG[ig], pi0, n);
         }         /* for(irf) */
      }            /* for(ig) */
      zero(pi0, n);
      for (k = 0; k < n; k++) for (ig = 0; ig < com.ngene; ig++)
         pi0[k] += com.piG[ig][k] / com.ngene;
      for (irf = 0; irf < nrf; irf++) {  /* com.pi[] */
         zero(com.pi, n);
         for (ig = 0; ig < com.ngene; ig++)  for (js = 0; js < com.ns; js++)
            AddFreqSeqGene(js, ig, pi0, com.pi);

         t = sum(com.pi, n);
         if (t < 1e-10) {
            for (k = 0; k < n; k++)  com.pi[k] = 1.0 / n;
            t = 1;
         }
         abyx(1 / sum(com.pi, n), com.pi, n);
         if (distance(com.pi, pi0, n) < 1e-8) break;
         xtoy(com.pi, pi0, n);
      }            /* for(ig) */
   }
   fprintf(fout, "\n\n%-*s", wname, "Average");
   for (k = 0; k < n; k++) fprintf(fout, "%8.5f", com.pi[k]);
   if (!com.cleandata) fputs("\n(Ambiguity characters are used to calculate freqs.)\n", fout);

   fprintf(fout, "\n\n# constant sites: %6d (%.2f%%)",
      nconstp, (double)nconstp*100. / com.ls);

   if (com.model == 0 || (com.seqtype == BASEseq && com.model == 1)) {
      fillxc(com.pi, 1. / n, n);
      for (ig = 0; ig < com.ngene; ig++)
         xtoy(com.pi, com.piG[ig], n);
   }
   if (com.seqtype == BASEseq && com.model == 5) { /* T92 model */
      com.pi[0] = com.pi[2] = (com.pi[0] + com.pi[2]) / 2;
      com.pi[1] = com.pi[3] = (com.pi[1] + com.pi[3]) / 2;
      for (ig = 0; ig < com.ngene; ig++) {
         com.piG[ig][0] = com.piG[ig][2] = (com.piG[ig][0] + com.piG[ig][2]) / 2;
         com.piG[ig][1] = com.piG[ig][3] = (com.piG[ig][1] + com.piG[ig][3]) / 2;
      }
   }

   /* this is used only for REV & REVu in baseml and model==3 in aaml */
   if (com.seqtype == AAseq) {
      for (k = 0, t = 0; k < n; k++)  t += (com.pi[k] > 0);
      if (t <= 4)
         puts("\n\a\t\tAre these a.a. sequences?");
   }
   if (com.cleandata && com.ngene == 1) {
      for (h = 0, lmax = -(double)com.ls*log((double)com.ls); h < com.npatt; h++)
         if (com.fpatt[h] > 1) lmax += com.fpatt[h] * log((double)com.fpatt[h]);
   }
   if (fout) {
      if (lmax) fprintf(fout, "\nln Lmax (unconstrained) = %.6f\n", lmax);
      fflush(fout);
   }

   free(pisg);
   return(0);
}


void AddFreqSeqGene(int js, int ig, double pi0[], double pi[])
{
   /* This adds the character counts in sequence js in gene ig to pi,
    using pi0, by resolving ambiguities.  The data are coded.  com.cleandata==1 or 0.
    This is for nucleotide and amino acid sequences only.
    */
    //char *pch = (com.seqtype==0 ? BASEs : (com.seqtype==2 ? AAs : BINs));
   int k, h, b, n = com.ncode;
   double t;

   if (com.cleandata) {
      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++)
         pi[com.z[js][h]] += com.fpatt[h];
   }
   else {
      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
         b = com.z[js][h];
         if (b < n)
            pi[b] += com.fpatt[h];
         else {
            /*
             if(com.seqtype==BASEseq) {
             NucListall(BASEs[b], &nb, ib);
             for(k=0,t=0; k<nb; k++) t += pi0[ib[k]];
             for(k=0; k<nb; k++)
             pi[ib[k]] += pi0[ib[k]]/t * com.fpatt[h];
             }
             */
            if (com.seqtype == BASEseq) {
               for (k = 0, t = 0; k < nChara[b]; k++)
                  t += pi0[(int)CharaMap[b][k]];
               for (k = 0; k < nChara[b]; k++)
                  pi[(int)CharaMap[b][k]] += pi0[(int)CharaMap[b][k]] / t * com.fpatt[h];
            }
            else if (com.seqtype == AAseq)  /* unrecognized AAs are treated as "?". */
               for (k = 0; k < n; k++) pi[k] += pi0[k] * com.fpatt[h];
         }
      }
   }
}


int RemoveIndel(void)
{
   /* Remove ambiguity characters and indels in the untranformed sequences,
    Changing com.ls and com.pose[] (site marks for multiple genes).
    For codonml, com.ls is still 3*#codons
    Called at the end of ReadSeq, when com.pose[] are still site marks.
    All characters in com.z[][] not found in the character string pch are
    considered ambiguity characters and are removed.
    */
   int  n = com.ncode, h, k, j, js, lnew, nindel, n31 = 1;
   char b, *miss;  /* miss[h]=1 if site (codon) h is missing, 0 otherwise */
   char *pch = ((com.seqtype <= 1 || com.seqtype == CODON2AAseq) ? BASEs : (com.seqtype == 2 ? AAs : BINs));

   if (com.seqtype == CODONseq || com.seqtype == CODON2AAseq) {
      n31 = 3; n = 4;
   }

   if (com.ls%n31) error2("ls in RemoveIndel.");
   if ((miss = (char*)malloc(com.ls / n31 * sizeof(char))) == NULL)
      error2("oom miss");
   for (h = 0; h < com.ls / n31; h++)
      miss[h] = 0;
   for (js = 0; js < com.ns; js++) {
      for (h = 0, nindel = 0; h < com.ls / n31; h++) {
         for (k = 0; k < n31; k++) {
            b = (char)toupper(com.z[js][h*n31 + k]);
            for (j = 0; j < n; j++)
               if (b == pch[j]) break;
            if (j == n) {
               miss[h] = 1; nindel++;
            }
         }
      }
      if (noisy > 2 && nindel)
         printf("\n%6d ambiguity characters in seq. %d", nindel, js + 1);
   }
   if (noisy > 2) {
      for (h = 0, k = 0; h < com.ls / n31; h++)  if (miss[h]) k++;
      printf("\n%d sites are removed. ", k);
      if (k < 1000)
         for (h = 0; h < com.ls / n31; h++)  if (miss[h]) printf(" %2d", h + 1);
   }

   for (h = 0, lnew = 0; h < com.ls / n31; h++) {
      if (miss[h]) continue;
      for (js = 0; js < com.ns; js++) {
         for (k = 0; k < n31; k++)
            com.z[js][lnew*n31 + k] = com.z[js][h*n31 + k];
      }
      com.pose[lnew] = com.pose[h];
      lnew++;
   }
   com.ls = lnew*n31;
   free(miss);
   return (0);
}


int Site2Pattern(FILE *fout)
{
   int h;
   fprintf(fout, "\n\nMapping site to pattern (i.e. site %d has pattern %d):\n",
      com.ls - 1, com.pose[com.ls - 2] + 1);
   for (h = 0; h < com.ls; h++) {
      fprintf(fout, "%6d", com.pose[h] + 1);
      if ((h + 1) % 10 == 0) FPN(fout);
   }
   FPN(fout);
   return (0);
}


#endif



int print1seq(FILE*fout, unsigned char *z, int ls, int pose[])
{
   /* This prints out one sequence, and the sequences are encoded.
    z[] contains patterns if (pose!=NULL)
    This uses com.seqtype.
    */
   int h, hp, gap = 10;
   char *pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   // char str[4]="";
   // int nb = (com.seqtype==CODONseq?3:1);

   for (h = 0; h < ls; h++) {
      hp = (pose ? pose[h] : h);
      if (com.seqtype != CODONseq) {
         fprintf(fout, "%c", pch[z[hp]]);
         if ((h + 1) % gap == 0) fputc(' ', fout);
      }
      else
         fprintf(fout, "%s ", CODONs[z[hp]]);
   }
   return(0);
}

void printSeqs(FILE *fout, unsigned char *z[], unsigned char *spnames[], int ns, int ls, int npatt, double fpatt[], int *pose, char keep[], int format)
{
   /* Print sequences into fout, using paml (format=0 or 1) or NEXUS (format=2)
      formats.
      Use pose=NULL if called before site patterns are collapsed.
      keep[] marks the sequences to be printed.  Use NULL for keep if all sequences
      are to be printed.
      Sequences may (com.cleandata==1) and may not (com.cleandata=0) be coded.
      z[] has site patterns if pose!=NULL.
      This uses com.seqtype, and ls is the number of codons for codon seqs.
      See notes in print1seq()

      format = 0,1: PAML sites or patterns
      2: NEXUS Nexus format.
      3: NEXUS Nexus JC69 format

      This is used by evolver.  Check and merge with printsma().
   */
   int h, i, j, ls1, n31 = (com.seqtype == 1 ? 3 : 1), nskept = ns, wname = 10;
   char *dt = (com.seqtype == AAseq ? "protein" : "dna");
   char *pch = (com.seqtype == 0 ? BASEs : AAs);

   ls1 = (format == 1 ? npatt : ls);
   if (keep)
      for (j = 0; j < ns; j++) nskept -= !keep[j];
   if (format == 0 || format == 1)
      fprintf(fout, "\n\n%d %d %s\n\n", nskept, ls1*n31, (format == 1 ? " P" : ""));
   else if (format == 2 || format == 3) {  /* NEXUS format */
      fprintf(fout, "\nbegin data;\n");
      fprintf(fout, "   dimensions ntax=%d nchar=%d;\n", nskept, ls1*n31);
      fprintf(fout, "   format datatype=%s missing=? gap=-;\n   matrix\n", dt);
   }

   for (j = 0; j < ns; j++, FPN(fout)) {
      if (keep && !keep[j]) continue;
      fprintf(fout, "%s%-*s  ", (format == 2 || format == 3 ? "      " : ""), wname, spnames[j]);
      if (format == 3) {     /* NEXUS Nexus JC69 format */
         for (h = 0, ls1 = 0; h < com.npatt; h++)
            for (i = 0; i < (int)com.fpatt[h]; i++) {
               fprintf(fout, "%c", pch[z[j][h]]);
               if (++ls1 % 10 == 0) fprintf(fout, " ");
            }
      }
      else
         print1seq(fout, z[j], (format == 1 ? npatt : ls), pose);
   }
   if (format == 2 || format == 3) fprintf(fout, "   ;\nend;");
   else if (format == 1) {
      for (h = 0, FPN(fout); h < com.npatt; h++) {
         /* fprintf(fout," %12.8f", com.fpatt[h]/(double)ls); */
         fprintf(fout, " %4.0f", com.fpatt[h]);
         if ((h + 1) % 15 == 0) FPN(fout);
      }
   }
   fprintf(fout, "\n\n");
   fflush(fout);
}


#define gammap(x,alpha) (alpha*(1-pow(x,-1.0/alpha)))
/* DistanceREV () used to be here, moved to pamp.
 */

#if (defined BASEML || defined BASEMLG || defined MCMCTREE || defined PROBTREE || defined YULETREE)

double SeqDivergence(double x[], int model, double alpha, double *kappa)
{
   /* alpha=0 if no gamma
    return -1 if in error.
    Check DistanceF84() if variances are wanted.
    */
   int i, j;
   double p[4], Y, R, a1, a2, b, P1, P2, Q, fd, tc, ag, GC;
   double d, smallv = 1e-10 / com.ls, largek = 999, larged = 9;

   if (testXMat(x)) {
      matout(F0, x, 4, 4);
      printf("\nfrequency matrix error, setting distance to large d");
      return(larged);
   }
   for (i = 0, fd = 1, zero(p, 4); i < 4; i++) {
      fd -= x[i * 4 + i];
      FOR(j, 4) { p[i] += x[i * 4 + j] / 2;  p[j] += x[i * 4 + j] / 2; }
   }
   P1 = x[0 * 4 + 1] + x[1 * 4 + 0];
   P2 = x[2 * 4 + 3] + x[3 * 4 + 2];
   Q = x[0 * 4 + 2] + x[0 * 4 + 3] + x[1 * 4 + 2] + x[1 * 4 + 3] + x[2 * 4 + 0] + x[2 * 4 + 1] + x[3 * 4 + 0] + x[3 * 4 + 1];
   if (fd < smallv)
      return(0);
   if (P1 < smallv) P1 = 0;
   if (P2 < smallv) P2 = 0;
   if (Q < smallv) Q = 0;
   Y = p[0] + p[1];    R = p[2] + p[3];  tc = p[0] * p[1]; ag = p[2] * p[3];

   switch (model) {
   case (JC69):
      FOR(i, 4) p[i] = .25;
   case (F81):
      for (i = 0, b = 0; i < 4; i++)  b += p[i] * (1 - p[i]);
      if (1 - fd / b <= 0) return (larged);

      if (alpha <= 0) return (-b*log(1 - fd / b));
      else return  (-b*gammap(1 - fd / b, alpha));
   case (K80):
      a1 = 1 - 2 * (P1 + P2) - Q;   b = 1 - 2 * Q;
      /*      if (a1<=0 || b<=0) return (-1); */
      if (a1 <= 0 || b <= 0) return (larged);
      if (alpha <= 0) { a1 = -log(a1);  b = -log(b); }
      else { a1 = -gammap(a1, alpha); b = -gammap(b, alpha); }
      a1 = .5*a1 - .25*b;  b = .25*b;
      if (b > smallv) *kappa = a1 / b; else *kappa = largek;
      return (a1 + 2 * b);
   case (F84):
      if (Y < smallv || R < smallv) {
         *kappa = -1;  d = larged;
      }
      else {
         a1 = (2 * (tc + ag) + 2 * (tc*R / Y + ag*Y / R)*(1 - Q / (2 * Y*R)) - P1 - P2) / (2 * tc / Y + 2 * ag / R);
         b = 1 - Q / (2 * Y*R);
         /*      if (a1<=0 || b<=0) return (-1); */
         if (a1 <= 0 || b <= 0) return (larged);
         if (alpha <= 0) { a1 = -log(a1); b = -log(b); }
         else { a1 = -gammap(a1, alpha); b = -gammap(b, alpha); }
         a1 = .5*a1;  b = .5*b;
         *kappa = a1 / b - 1;
         *kappa = max2(*kappa, -.5);
         d = 4 * b*(tc*(1 + *kappa / Y) + ag*(1 + *kappa / R) + Y*R);
      }
      return d;
   case (HKY85):         /* HKY85, from Rzhetsky & Nei (1995 MBE 12, 131-51) */
      if (0 && Y < smallv || R < smallv) {
         *kappa = -1;  d = larged;
      }
      else {
         *kappa = largek;
         a1 = 1 - Y*P1 / (2 * tc) - Q / (2 * Y);
         a2 = 1 - R*P2 / (2 * ag) - Q / (2 * R);
         b = 1 - Q / (2 * Y*R);
         if (a1 <= 0 || a2 <= 0 || b <= 0) return (larged);
         if (alpha <= 0) { a1 = -log(a1); a2 = -log(a2); b = -log(b); }
         else { a1 = -gammap(a1, alpha); a2 = -gammap(a2, alpha); b = -gammap(b, alpha); }
         a1 = -R / Y*b + a1 / Y;
         a2 = -Y / R*b + a2 / R;
         if (b > 0) *kappa = min2((a1 + a2) / (2 * b), largek);
         d = 2 * (p[0] * p[1] + p[2] * p[3])*(a1 + a2) / 2 + 2 * Y*R*b;
      }
      return d;
   case (T92):
      *kappa = largek;
      GC = p[1] + p[3];
      a1 = 1 - Q - (P1 + P2) / (2 * GC*(1 - GC));   b = 1 - 2 * Q;
      if (a1 <= 0 || b <= 0) return (larged);
      if (alpha <= 0) { a1 = -log(a1); b = -log(b); }
      else { a1 = -gammap(a1, alpha); b = -gammap(b, alpha); }
      if (Q > 0) *kappa = 2 * a1 / b - 1;
      return 2 * GC*(1 - GC)*a1 + (1 - 2 * GC*(1 - GC)) / 2 * b;
   case (TN93):         /* TN93  */
      if (0 && Y < smallv || R < smallv || tc < smallv || ag < smallv) {
         *kappa = -1;  d = larged;
      }
      else {
         a1 = 1 - Y*P1 / (2 * tc) - Q / (2 * Y);
         a2 = 1 - R*P2 / (2 * ag) - Q / (2 * R);
         b = 1 - Q / (2 * Y*R);
         /*      if (a1<=0 || a2<=0 || b<=0) return (-1); */
         if (a1 <= 0 || a2 <= 0 || b <= 0) return (larged);
         if (alpha <= 0) { a1 = -log(a1); a2 = -log(a2); b = -log(b); }
         else { a1 = -gammap(a1, alpha); a2 = -gammap(a2, alpha); b = -gammap(b, alpha); }
         a1 = .5 / Y*(a1 - R*b);  a2 = .5 / R*(a2 - Y*b);  b = .5*b;
         *kappa = largek;
         /*
         printf("\nk1&k2 = %.6f %.6f\n", a1/b,a2/b);
         */
         if (b > 0) *kappa = min2((a1 + a2) / (2 * b), largek);
         d = 4 * p[0] * p[1] * a1 + 4 * p[2] * p[3] * a2 + 4 * Y*R*b;
      }
      return d;
   }
   return (-1);
}


double DistanceIJ(int is, int js, int model, double alpha, double *kappa)
{
   /* Distance between sequences is and js.
    See DistanceMatNuc() for more details.
    */
   char b0, b1;
   int h, n = 4, missing = 0;
   double x[16], sumx, larged = 9;

   zero(x, 16);
   if (com.cleandata && com.seqtype == 0) {
      for (h = 0; h < com.npatt; h++)
         x[com.z[is][h] * n + com.z[js][h]] += com.fpatt[h];
   }
   else {
      for (h = 0; h < com.npatt; h++) {
         b0 = com.z[is][h];
         b1 = com.z[js][h];
         if (b0 < n && b1 < n)
            x[b0*n + b1] += com.fpatt[h];
         else
            missing = 1;
      }
   }
   sumx = sum(x, 16);

   if (sumx <= 0) return(larged);    /* questionable??? */
   abyx(1. / sum(x, 16), x, 16);
   return SeqDivergence(x, model, alpha, kappa);
}




extern int nR;
extern double Cijk[], Root[];

int QTN93(int model, double Q[], double kappa1, double kappa2, double pi[])
{
   int i, j;
   double T = pi[0], C = pi[1], A = pi[2], G = pi[3], Y = T + C, R = A + G, scalefactor;

   if (model == JC69 || model == F81) kappa1 = kappa2 = com.kappa = 1;
   else if (com.model < TN93)       kappa2 = kappa1;
   if (model == F84) { kappa2 = 1 + kappa1 / R; kappa1 = 1 + kappa1 / Y; }
   scalefactor = 1 / (2 * T*C*kappa1 + 2 * A*G*kappa2 + 2 * Y*R);

   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) Q[i * 4 + j] = (i == j ? 0 : 1);
   Q[0 * 4 + 1] = Q[1 * 4 + 0] = kappa1;
   Q[2 * 4 + 3] = Q[3 * 4 + 2] = kappa2;
   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) Q[i * 4 + j] *= pi[j] * scalefactor;
   for (i = 0; i < 4; i++) { Q[i * 4 + i] = 0;  Q[i * 4 + i] = -sum(Q + i * 4, 4); }

   return (0);
}

int RootTN93(int model, double kappa1, double kappa2, double pi[],
   double *scalefactor, double Root[])
{
   double T = pi[0], C = pi[1], A = pi[2], G = pi[3], Y = T + C, R = A + G;

   if (model == JC69 || model == F81) kappa1 = kappa2 = com.kappa = 1;
   else if (com.model < TN93)       kappa2 = kappa1;
   if (model == F84) { kappa2 = 1 + kappa1 / R; kappa1 = 1 + kappa1 / Y; }

   *scalefactor = 1 / (2 * T*C*kappa1 + 2 * A*G*kappa2 + 2 * Y*R);

   Root[0] = 0;
   Root[1] = -(*scalefactor);
   Root[2] = -(Y + R*kappa2) * (*scalefactor);
   Root[3] = -(Y*kappa1 + R) * (*scalefactor);
   return (0);
}


int eigenTN93(int model, double kappa1, double kappa2, double pi[],
   int *nR, double Root[], double Cijk[])
{
   /* initialize Cijk[] & Root[], which are the only part to be changed
    for a new substitution model
    for JC69, K80, F81, F84, HKY85, TN93
    Root: real Root divided by v, the number of nucleotide substitutions.
    */
   int i, j, k, nr;
   double scalefactor, U[16], V[16], t;
   double T = pi[0], C = pi[1], A = pi[2], G = pi[3], Y = T + C, R = A + G;

   if (model == JC69 || model == F81) kappa1 = kappa2 = com.kappa = 1;
   else if (com.model < TN93)       kappa2 = kappa1;
   RootTN93(model, kappa1, kappa2, pi, &scalefactor, Root);

   *nR = nr = 2 + (model == K80 || model >= F84) + (model >= HKY85);
   U[0 * 4 + 0] = U[1 * 4 + 0] = U[2 * 4 + 0] = U[3 * 4 + 0] = 1;
   U[0 * 4 + 1] = U[1 * 4 + 1] = 1 / Y;   U[2 * 4 + 1] = U[3 * 4 + 1] = -1 / R;
   U[0 * 4 + 2] = U[1 * 4 + 2] = 0;  U[2 * 4 + 2] = G / R;  U[3 * 4 + 2] = -A / R;
   U[2 * 4 + 3] = U[3 * 4 + 3] = 0;  U[0 * 4 + 3] = C / Y;  U[1 * 4 + 3] = -T / Y;

   xtoy(pi, V, 4);
   V[1 * 4 + 0] = R*T;   V[1 * 4 + 1] = R*C;
   V[1 * 4 + 2] = -Y*A;  V[1 * 4 + 3] = -Y*G;
   V[2 * 4 + 0] = V[2 * 4 + 1] = 0;  V[2 * 4 + 2] = 1;   V[2 * 4 + 3] = -1;
   V[3 * 4 + 0] = 1;  V[3 * 4 + 1] = -1;   V[3 * 4 + 2] = V[3 * 4 + 3] = 0;

   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
      Cijk[i * 4 * nr + j*nr + 0] = U[i * 4 + 0] * V[0 * 4 + j];
      switch (model) {
      case JC69:
      case F81:
         for (k = 1, t = 0; k < 4; k++) t += U[i * 4 + k] * V[k * 4 + j];
         Cijk[i * 4 * nr + j*nr + 1] = t;
         break;
      case K80:
      case F84:
         Cijk[i * 4 * nr + j*nr + 1] = U[i * 4 + 1] * V[1 * 4 + j];
         for (k = 2, t = 0; k < 4; k++) t += U[i * 4 + k] * V[k * 4 + j];
         Cijk[i * 4 * nr + j*nr + 2] = t;
         break;
      case HKY85:   case T92:   case TN93:
         for (k = 1; k < 4; k++)  Cijk[i * 4 * nr + j*nr + k] = U[i * 4 + k] * V[k * 4 + j];
         break;
      default:
         error2("model in eigenTN93");
      }
   }
   return(0);
}


#endif


#ifdef NODESTRUCTURE


double CountTrees(int ns, int rooted)
{
   double i, ntree = 1;

   if (ns > 70) error2("ns too large in NumberofTrees().");
   for (i = 4; i <= ns + rooted; i++)
      ntree *= 2 * i - 5;
   return (ntree);
}

double CountLHs(int ns)
{
   double i, nLH = 1;

   if (ns > 70) error2("ns too large in NumberofLHs().");
   for (i = 3; i <= ns; i++)
      nLH *= i*(i - 1) / 2;
   return (nLH);
}



void BranchToNode(void)
{
   /* tree.root need to be specified before calling this
    */
   int i, from, to;

   tree.nnode = tree.nbranch + 1;
   for (i = 0; i < tree.nnode; i++)
   {
      nodes[i].father = nodes[i].ibranch = -1;  nodes[i].nson = 0;
   }
   for (i = 0; i < tree.nbranch; i++) {
      from = tree.branches[i][0];
      to = tree.branches[i][1];
      nodes[from].sons[nodes[from].nson++] = to;
      nodes[to].father = from;
      nodes[to].ibranch = i;
   }
   /*  nodes[tree.root].branch=0;  this breaks method=1 */
}

void NodeToBranchSub(int inode);

void NodeToBranchSub(int inode)
{
   int i, ison;

   for (i = 0; i < nodes[inode].nson; i++) {
      tree.branches[tree.nbranch][0] = inode;
      tree.branches[tree.nbranch][1] = ison = nodes[inode].sons[i];
      nodes[ison].ibranch = tree.nbranch++;
      if (nodes[ison].nson > 0)  NodeToBranchSub(ison);
   }
}

void NodeToBranch(void)
{
   tree.nbranch = 0;
   NodeToBranchSub(tree.root);
   if (tree.nnode != tree.nbranch + 1)
      error2("nnode != nbranch + 1?");
}


void ClearNode(int inode)
{
   /* a source of confusion. Try not to use this routine.
    */
   nodes[inode].father = nodes[inode].ibranch = -1;
   nodes[inode].nson = 0;
   nodes[inode].branch = nodes[inode].age = 0;
   /* nodes[inode].label = -1; */
   /* nodes[inode].branch = 0; clear node structure only, not branch lengths */
   /* for(i=0; i<com.ns; i++) nodes[inode].sons[i]=-1; */
}

int ReadTreeB(FILE *ftree, int popline)
{
   char line[255];
   int nodemark[NS * 2 - 1] = { 0 }; /* 0: absent; 1: father only (root); 2: son */
   int i, j, state = 0, YoungAncestor = 0;

   if (com.clock) {
      puts("\nbranch representation of tree might not work with clock model.");
      getchar();
   }

   fscanf(ftree, "%d", &tree.nbranch);
   for (j = 0; j < tree.nbranch; j++) {
      for (i = 0; i < 2; i++) {
         if (fscanf(ftree, "%d", &tree.branches[j][i]) != 1) state = -1;
         tree.branches[j][i]--;
         if (tree.branches[j][i]<0 || tree.branches[j][i]>com.ns * 2 - 1)
            error2("ReadTreeB: node numbers out of range");
      }
      nodemark[tree.branches[j][1]] = 2;
      if (nodemark[tree.branches[j][0]] != 2) nodemark[tree.branches[j][0]] = 1;
      if (tree.branches[j][0] < com.ns)  YoungAncestor = 1;

      printf("\nBranch #%3d: %3d -> %3d", j + 1, tree.branches[j][0] + 1, tree.branches[j][1] + 1);

   }
   if (popline) fgets(line, 254, ftree);
   for (i = 0, tree.root = -1; i < tree.nbranch; i++)
      if (nodemark[tree.branches[i][0]] != 2) tree.root = tree.branches[i][0];
   if (tree.root == -1) error2("root err");
   for (i = 0; i < com.ns; i++)
      if (nodemark[i] == 0) {
         matIout(F0, nodemark, 1, com.ns);
         error2("branch specification of tree");
      }

   if (YoungAncestor) {
      puts("\nAncestors in the data?  Take care.");
      if (!com.cleandata) {
         puts("This kind of tree does not work with unclean data.");
         getchar();
      }
   }

   BranchToNode();
   return (state);
}


int OutTreeB(FILE *fout)
{
   int j;
   char *fmt[] = { " %3d..%-3d", " %2d..%-2d" };
   FOR(j, tree.nbranch)
      fprintf(fout, fmt[0], tree.branches[j][0] + 1, tree.branches[j][1] + 1);
   return (0);
}

int GetTreeFileType(FILE *ftree, int *ntree, int *pauptree, int shortform);

int GetTreeFileType(FILE *ftree, int *ntree, int *pauptree, int shortform)
{
   /* paupstart="begin trees" and paupend="translate" identify paup tree files.
    paupch=";" will be the last character before the list of trees.
    Modify if necessary.
    */
   int i, k, lline = 32000, ch = 0, paupch = ';';
   char line[32000];
   char *paupstart = "begin tree", *paupend = "translate";

   *pauptree = 0;
   k = fscanf(ftree, "%d%d", &i, ntree);
   if (k == 2) {
      if (i == com.ns)  return(0);                 /* old paml style */
      else              error2("Number of sequences different in tree and seq files.");
   }
   else if (k == 1) { *ntree = i; return(0); }           /* phylip & molphy style */
   while (ch != '(' && !isalnum(ch) && ch != EOF)  ch = fgetc(ftree);  /* treeview style */
   if (ch == '(') { *ntree = -1; ungetc(ch, ftree); return(0); }

   puts("\n# seqs in tree file does not match.  Read as the nexus format.");
   for (; ; ) {
      if (fgets(line, lline, ftree) == NULL) error2("tree err1: EOF");
      strcase(line, 0);
      if (strstr(line, paupstart)) { *pauptree = 1; *ntree = -1; break; }
   }
   if (shortform) return(0);
   for (; ; ) {
      if (fgets(line, lline, ftree) == NULL) error2("tree err2: EOF");
      strcase(line, 0);
      if (strstr(line, paupend)) break;
   }
   for (; ; ) {
      if ((ch = fgetc(ftree)) == EOF) error2("tree err3: EOF");
      if (ch == paupch) break;
   }
   if (fgets(line, lline, ftree) == NULL) error2("tree err4: EOF");

   return(0);
}

int PopPaupTree(FILE *ftree);
int PopPaupTree(FILE *ftree)
{
/* This reads out the string in front of the tree in the nexus format,
   typically "tree PAUP_1 = [&U]" with "[&U]" optional
*/
   int ch;

   for (; ;) {
      ch = fgetc(ftree);
      if (ch == '(')
      {
         ungetc(ch, ftree); return(0);
      }
      else if (ch == EOF || ch == '/')
         return(-1);
   }
   return(0);
}


void DownTreeCladeLabel(int inode, int label)
{
/* This goes down the tree to change $ labels (stored in nodes[].label2) into
   # labels (stored in nodes[].label).  To deal with nested clade labels,
   branches within a clade are labeled by negative numbers initially, and
   converted to positive labels at the end of the algorithm.
   nodes[].label and nodes[].label2 are initialized to -1 before this routine is called.
*/
   int i;

   if (nodes[inode].label2 != -1)
      label = (int)nodes[inode].label2;
   if (inode != tree.root && nodes[inode].label == -1)
      nodes[inode].label = label;
   for (i = 0; i < nodes[inode].nson; i++)
      DownTreeCladeLabel(nodes[inode].sons[i], label);
}

int IsNameNumber(char line[])
{
   /* returns 0 if line has species number; 1 if it has name.
    It uses com.ns.
    */
   int isname = 1, alldigits = 1, n;
   char *p = line;

   while (*p)
      if (!isdigit(*p++)) { alldigits = 0; break; }
   if (alldigits) {
      n = atoi(line);
      if (n >= 1 && n <= com.ns) isname = 0;
   }
   return(isname);
}


int ReadUntil(FILE *fin, char line[], char delimiters[], int maxline)
{
   /* read into line[] from fin until any character in delimiters[] or maxline,
      & trim spaces at the end.
   */
   int i;

   for (i = 0; i < maxline; i++) {  /* read notes into line[] */
      line[i] = (char)fgetc(fin);
      if (line[i] == EOF)
         error2("EOF when reading node label");
      if (strchr(delimiters, (int)line[i])) {
         ungetc(line[i], fin);
         line[i] = '\0';
         break;
      }
   }
   for (; i > 0; i--) { /* trim spaces*/
      if (isgraph(line[i]))
         break;
      else
         line[i] = 0;
   }
   return(i + 1);
}



/* Read a newick tree from ftree.
Branch lengths are read in nodes[].branch, and branch (node) labels
(integers) are preceeded by # and read in nodes[].label.  If the clade label
$ is used, the label is read nodes[].label2 first and then moved into
nodes[].label in the routine DownTreeCladeLabel().

Calibration information for mcmctree may be read into nodes[].branch and nodes[].label,
as well as nodes[].NodeStr, and is processed inside mcmctree.
*haslength is set to 1 (branch lengths), 2 (calibration info) or 3 (both).
However, the bit for calibrations is set only if the symbols > < exist and not for
calibrations specified using L, U, G, etc, which will be stored in nodes[].annotation
and processed using ProcessNodeAnnotation() in mcmctree.

This assumes that com.ns is known.  Names are considered case-sensitive, with trailing spaces ignored.

copyname = 0: species numbers and names are both accepted, but names have to match the names
in com.spname[], which are from the sequence data file.
Used by baseml and codeml, for example.
1: species names are copied into com.spname[], but species numbers are accepted.
Used by evolver for simulation, in which case no species names were read before.
2: the tree must have species names, which are copied into com.spname[].

isname = 0:   species number; 1: species name;
*/
int ReadTreeN(FILE *ftree, int *haslength, int copyname, int popline)
{
   int cnode, cfather = -1;  /* current node and father */
   int inode = 0;  /* node number that will have the next branch length */
   int cladeLabels = 0, i, k, level = 0, isname, ch = ' ', tip = 0;
   /* delimiters for names & labels (labels if not in quotes) */
   char check[NS], delimitersN[] = "(),:;#[", delimitersL[] = "(),:;", quote[2][4] = { "\"\'[", "\"\']" };
   int lline = 32000, namelabel;  /* namelabel = 0 for name 1 for label */
   char line[32000], *pch;
   char debug = 0;

   /*
   (((A1, B1) ' #1 B(0.1, 0.3)', (A2, B2) #1), C) 'B(0.9, 1.1)';
   ((A)H#H[&gamma = 0:3];H#H)R;
   ((A : 0:02; (B : 0:01)H#H1[&gamma=0:3] : 0:01)S : 0:03; (H#H1 : 0:02;C : 0:03)T : 0:02)R :0:03;
   */
   if (com.ns <= 0)  error2("you should specify # seqs in the tree file.");

   for (i = 0; i < com.ns; i++) check[i] = 0;
   *haslength = 0;
   tree.nnode = com.ns;  tree.nbranch = 0;
   for (i = 0; i < 2 * com.ns - 1; i++) {
      nodes[i].father = nodes[i].ibranch = -1;
      nodes[i].nson = 0;  nodes[i].label = nodes[i].label2 = -1;  nodes[i].branch = 0;
      nodes[i].age = 0;  /* For TipDate models this is set later for each tree. */
      nodes[i].latent = 0;
      nodes[i].annotation = NULL;
   }
   while (isspace(ch))
      ch = fgetc(ftree);  /* skip spaces */
   ungetc(ch, ftree);
   if (isdigit(ch)) {
      ReadTreeB(ftree, popline); return(0);
   }
   if (PopPaupTree(ftree) == -1) return(-1);

   for (; ; ) {
      ch = fgetc(ftree);
      if (ch == EOF) return(-1);
      else if (ch == ';') {
         if (level != 0) error2("; in treefile");
         else            break;
      }
      else if (ch == ',') namelabel = 0;
      else if (!isgraph(ch))
         continue;
      else if (ch == '(') {       /* left (  */
         level++;
         namelabel = 0;
         cnode = tree.nnode++;
         if (tree.nnode > 2 * com.ns - 1)
            error2("check #seqs and tree: perhaps too many '('?");
         if (cfather >= 0) {
            if (nodes[cfather].nson >= MAXNSONS) {
               printf("there are at least %d daughter nodes, raise MAXNSONS?", nodes[cfather].nson);
               exit(-1);
            }
            nodes[cfather].sons[nodes[cfather].nson++] = cnode;
            nodes[cnode].father = cfather;
            tree.branches[tree.nbranch][0] = cfather;
            tree.branches[tree.nbranch][1] = cnode;
            nodes[cnode].ibranch = tree.nbranch++;
         }
         else
            tree.root = cnode;
         cfather = cnode;
      }
      /* treating : and > in the same way is risky. */
      else if (ch == ')') {
         level--;  namelabel = 1;  inode = cfather;  cfather = nodes[cfather].father;
      }
      else if (ch == ':') {
         *haslength = 1;
         fscanf(ftree, "%lf", &nodes[inode].branch);
      }
      else if (namelabel == 0) { /* read species name or number for tip */
         if (level <= 0) error2("expecting ; in the tree file");
         ungetc(ch, ftree);
         ReadUntil(ftree, line, delimitersN, lline);
         isname = IsNameNumber(line);
         if (isname == 0) {     /* number */
            if (copyname == 2) error2("Use names in tree.");
            sscanf(line, "%d", &cnode);
            cnode--;
         }
         else {                 /* name */
            if (copyname == 0) {
               for (i = 0; i < com.ns; i++)
                  if (!strcmp(line, com.spname[i])) break;
               if ((cnode = i) == com.ns) {
                  printf("\nSpecies %s?\n", line); exit(-1);
               }
            }
            else {
               if (tip > com.ns - 1) {
                  error2("error in tree: too many species in tree");
               }
               strcpy(com.spname[cnode = tip++], line);
            }
         }
         nodes[cnode].father = cfather;
         if (nodes[cfather].nson >= MAXNSONS)
            error2("too many daughter nodes, raise MAXNSONS");

         nodes[cfather].sons[nodes[cfather].nson++] = cnode;
         tree.branches[tree.nbranch][0] = cfather;
         tree.branches[tree.nbranch][1] = cnode;
         nodes[cnode].ibranch = tree.nbranch++;
         inode = cnode;
         check[cnode]++;
         namelabel = 1;  /* it is possile that a name is followed by a label */
      }
      else {     /* label or annotation (node can be tip or internal node) */
         if (strchr(quote[0], ch)) {
            k = ReadUntil(ftree, line, quote[1], lline);
            ch = fgetc(ftree);   /* pop off ' "  */
         }
         else {
            ungetc(ch, ftree);
            k = ReadUntil(ftree, line, delimitersL, lline);
         }
         nodes[inode].annotation = (char*)malloc((k + 1) * sizeof(char));
         if (nodes[inode].annotation == NULL) error2("oom annotation");
         strcpy(nodes[inode].annotation, line);
         /* This line is used in MCcoal.  ProcessNodeAnnotation is used to process node annotation elsewhere. */
         if (line[0] == '#')        /* branch labels used in baseml and codeml */
            sscanf(line + 1, "%lf", &nodes[inode].label);
         else if (line[0] == '$')   /* clade labels used to label branches in baseml and codeml */
            sscanf(line + 1, "%lf", &nodes[inode].label2);

         else if (line[0] == '@' || line[0] == '=') {
            sscanf(line + 1, "%lf", &nodes[inode].age);
         }
      }
   }   /* for( ; ; ) */

   if (popline)
      fgets(line, lline, ftree);
   for (i = 0; i < com.ns; i++) {
      if (check[i] > 1) {
         printf("\nSeq/species #%d (%s) occurs more than once in the tree\n", i + 1, com.spname[i]); exit(-1);
      }
      else if (check[i] < 1) {
         printf("\nSeq #%d (%s) is missing in the tree\n", i + 1, com.spname[i]);
         exit(-1);
      }
   }
   if (tree.nbranch > 2 * com.ns - 2) {
      printf("nbranch %d", tree.nbranch); puts("too many branches in tree?");
   }
   if (tree.nnode != tree.nbranch + 1) {
      printf("\nnnode%6d != nbranch%6d + 1\n", tree.nnode, tree.nbranch);
      exit(-1);
   }

   if (debug) PrintTree(1);
   return (0);
}

int OutSubTreeN(FILE *fout, int inode, int spnames, int printopt)
{
   int i, dad, nsib;

   if (inode > com.ns * 2 - 1)
      error2("inode large?");
   dad = nodes[inode].father;
   nsib = (inode == tree.root ? 0 : nodes[dad].nson);

   if (inode != tree.root && inode == nodes[dad].sons[0])
      fputc('(', fout);

   for (i = 0; i < nodes[inode].nson; i++)
      OutSubTreeN(fout, nodes[inode].sons[i], spnames, printopt);

   if (nodes[inode].nson == 0) { /* inode is tip */
      if (spnames) {
         if (printopt&PrNodeNum) fprintf(fout, "%d_", inode + 1);
         fprintf(fout, "%s", com.spname[inode]);
      }
      else
         fprintf(fout, "%d", inode + 1);
   }
   if ((printopt & PrNodeNum) && nodes[inode].nson)
      fprintf(fout, " %d ", inode + 1);

#if (defined BPP)
   if ((printopt & PrNodeStr) && inode < com.ns && nodes[inode].label > 0)
      fprintf(fout, " [&theta=%.6f]", nodes[inode].label);
#else
   if ((printopt & PrLabel) && nodes[inode].label > 0)
      fprintf(fout, " #%.6f", nodes[inode].label);
   /* fprintf(fout, " [&label=%.6f]", nodes[inode].label); */
#endif

   if ((printopt & PrAge) && nodes[inode].age)
      fprintf(fout, " @%.6f", nodes[inode].age);

   /*  Add branch labels to be read by Rod Page's TreeView.  */
#if (defined CODEML)
   if ((printopt & PrOmega) && inode != tree.root)
      fprintf(fout, " #%.6g ", nodes[inode].omega);
#elif (defined (EVOLVER) || defined (MCMCTREE) || defined (BPP))
   if ((printopt & PrNodeStr) && inode >= com.ns && nodes[inode].annotation)
      fprintf(fout, " %s", nodes[inode].annotation);
    else if ((printopt & PrNodeStr) && nodes[inode].latent && nodes[inode].annotation)
         fprintf(fout, " %s", nodes[inode].annotation);
#endif

   if ((printopt & PrBranch) && (inode != tree.root || nodes[inode].branch > 0))
      fprintf(fout, ": %.6f", nodes[inode].branch);
   if (nsib == 0)            /* root */
      fputc(';', fout);
   else if (inode == nodes[dad].sons[nsib - 1])  /* last sib */
      fputc(')', fout);
   else                     /* not last sib */
      fprintf(fout, ", ");

   return (0);
}


int OutTreeN(FILE *fout, int spnames, int printopt)
{
   OutSubTreeN(fout, tree.root, spnames, printopt);
   return(0);
}


int DeRoot(void)
{
   /* This cnages the bifurcation at the root into a trifurcation, but setting one of
    the sons to be the new root.  The new root is the first son that is not a tip.
    tree.nnode is updated, but the routine does not re-number the nodes, so the new
    node labels do not go from ns, ns + 1, ..., as they normally should.
    */
   int i, ison, sib, root = tree.root;

   if (nodes[root].nson != 2) error2("in DeRoot?");

   ison = nodes[root].sons[i = 0];
   if (nodes[ison].nson == 0)
      ison = nodes[root].sons[i = 1];
   sib = nodes[root].sons[1 - i];
   nodes[sib].branch += nodes[ison].branch;
   nodes[sib].father = tree.root = ison;
   nodes[tree.root].father = -1;
   nodes[tree.root].sons[nodes[tree.root].nson++] = sib;  /* sib added as the last child of the new root */
   nodes[tree.root].branch = 0;
   tree.nnode--;  /* added 2007/4/9 */
   return(0);
}

int PruneSubTreeN(int inode, int keep[])
{
   /* This prunes tips from the tree, using keep[com.ns].  Removed nodes in the
    big tree has nodes[].father=-1 and nodes[].nson=0.
    Do not change nodes[inode].nson and nodes[inode].sons[] until after the
    node's descendent nodes are all processed.  So when a son is deleted,
    only the father node's nson is changed, but not
    */
   int i, j, ison, father = nodes[inode].father, nson0 = nodes[inode].nson;

   nodes[inode].label = 0;
   for (i = 0; i < nson0; i++)
      PruneSubTreeN(nodes[inode].sons[i], keep);

   /* remove inode because of no descendents.
    Note this does not touch the father node */
   if (inode < com.ns && keep[inode] == 0)
      nodes[inode].father = -1;
   else if (inode >= com.ns) {
      for (i = 0, nodes[inode].nson = 0; i < nson0; i++) {
         ison = nodes[inode].sons[i];
         if (nodes[ison].father != -1)
            nodes[inode].sons[nodes[inode].nson++] = nodes[inode].sons[i];
      }
      if (nodes[inode].nson == 0)
         nodes[inode].father = -1;
   }

   /* remove inode if it has a single descendent ison */
   if (inode >= com.ns && nodes[inode].nson == 1 && inode != tree.root) {
      ison = nodes[inode].sons[0];
      nodes[ison].father = father;
      nodes[ison].branch += nodes[inode].branch;
      nodes[ison].label++;  /* records # deleted nodes for branch ison */
      for (j = 0; j < nodes[father].nson; j++) {
         if (nodes[father].sons[j] == inode)
         {
            nodes[father].sons[j] = ison;  break;
         }
      }
      nodes[inode].nson = 0;
      nodes[inode].father = -1;
   }
   else if (nodes[inode].nson == 1 && inode == tree.root) { /* move down root if root has 1 descendent */
      nodes[inode].father = -1;
      nodes[inode].nson = 0;
      ison = nodes[tree.root].sons[0];
      tree.root = ison;
      nodes[tree.root].father = -1;
      nodes[tree.root].branch = 0;
   }

   /*
    printf("\nVisiting inode %d\n", inode);
    FOR(i, tree.nnode) printf(" %2d", i);  FPN(F0);
    FOR(i, tree.nnode) printf(" %2.0f", nodes[i].label); FPN(F0);
    */
   return(0);
}


int GetSubTreeN(int keep[], int space[])
{
   /* This removes some tips to generate the subtree.  Branch lengths are
    preserved by summing them up when some nodes are removed.
    The algorithm use post-order tree traversal to remove tips and nodes.  It
    then switches to the branch representation to renumber nodes.
    space[] can be NULL.  If not, it returns newnodeNO[], which holds the
    new node numbers; for exmaple, newnodeNO[12]=5 means that old node 12 now
    becomes node 5.

    The routine does not change com.ns or com.spname[], which have to be updated
    outside.

    CHANGE OF ROOT happens if the root in the old tree had >=3 sons, but has 2
    sons in the new tree and if (!com.clock).  In that case, the tree is derooted.

    This routine does not work if a current seq is ancestral to some others
    and if that sequence is removed. (***check this comment ***)

    Different formats for keep[] are used.  Suppose the current tree is for
    nine species: a b c d e f g h i.

    (A) keep[]={1,0,1,1,1,0,0,1,0} means that a c d e h are kept in the tree.
    The old tip numbers are not changed, so that OutTreeN(?,1,?) gives the
    correct species names or OutTreeN(?,0,?) gives the old species numbers.

    (B) keep[]={1,0,2,3,4,0,0,5,0} means that a c d e h are kept in the tree, and
    they are renumbered 0 1 2 3 4 and all the internal nodes are renumbered
    as well to be consecutive.  Note that the positive numbers have to be
    consecutive natural numbers.

    keep[]={5,0,2,1,4,0,0,3,0} means that a c d e h are kept in the tree.
    However, the order of the sequences are changed to d c h e a, so that the
    numbers are now 0 1 2 3 4 for d c h e a.  This is useful when the subtree
    is extracted from a big tree for a subset of the sequence data, while the
    species are odered d c h e a in the sequence data file.
    This option can be used to renumber the tips in the complete tree.
    */
   int nsnew, i, j, k, nnode0 = tree.nnode, sumnumber = 0, newnodeNO[2 * NS - 1], ison, sib;
   int unrooted = (nodes[tree.root].nson >= 3);  /* com.clock is not checked here! */
   double *branch0;
   int debug = 0;

   if (debug) { FOR(i, com.ns) printf("%-30s %2d\n", com.spname[i], keep[i]); }
   for (i = 0, nsnew = 0; i < com.ns; i++)
      if (keep[i]) { nsnew++; sumnumber += keep[i]; }
   if (nsnew < 2)  return(-1);

   /* mark removed nodes in the big tree by father=-1 && nson=0.
    nodes[].label records the number of nodes collapsed.
    */
   PruneSubTreeN(tree.root, keep);
   /* If unrooted tree has a bifurcation at the new root, collapse root.  */
   if (unrooted && nodes[tree.root].nson == 2) {
      ison = nodes[tree.root].sons[i = 0];
      if (nodes[ison].nson == 0)
         ison = nodes[tree.root].sons[i = 1];
      sib = nodes[tree.root].sons[1 - i];

      nodes[sib].branch += nodes[ison].branch;
      nodes[sib].label += nodes[ison].label + 2;
      nodes[sib].father = tree.root = ison;
      nodes[tree.root].father = -1;
      nodes[tree.root].sons[nodes[tree.root].nson++] = sib;  /* sib added as the last child of the new root */
      nodes[tree.root].branch = 0;
   }
   if (debug) PrintTree(1);

   for (i = 0, k = 1; i < tree.nnode; i++) if (nodes[i].father != -1) k++;
   tree.nnode = k;
   NodeToBranch();

   /* to renumber the nodes */
   if (sumnumber > nsnew) {
      if (sumnumber != nsnew*(nsnew + 1) / 2)
         error2("keep[] not right in GetSubTreeN");

      if ((branch0 = (double*)malloc(nnode0 * sizeof(double))) == NULL) error2("oom#");
      FOR(i, nnode0) branch0[i] = nodes[i].branch;
      FOR(i, nnode0) newnodeNO[i] = -1;
      FOR(i, com.ns) if (keep[i]) newnodeNO[i] = keep[i] - 1;

      newnodeNO[tree.root] = k = nsnew;
      tree.root = k++;
      for (; i < nnode0; i++) {
         if (nodes[i].father == -1) continue;
         for (j = 0; j < tree.nbranch; j++)
            if (i == tree.branches[j][1]) break;
         if (j == tree.nbranch)
            error2("strange here");
         newnodeNO[i] = k++;
      }
      for (j = 0; j < tree.nbranch; j++) FOR(i, 2)
         tree.branches[j][i] = newnodeNO[tree.branches[j][i]];
      BranchToNode();
      for (i = 0; i < nnode0; i++) {
         if (newnodeNO[i] > -1)
            nodes[newnodeNO[i]].branch = branch0[i];
      }
      free(branch0);
   }

   if (space) memmove(space, newnodeNO, (com.ns * 2 - 1) * sizeof(int));
   return (0);
}


void PrintTree(int timebranches)
{
   int i, j;

   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%6s%6s", "dad", "node");
   if (timebranches)  printf("%10s%10s%10s", "age", "branch", "label");
   printf(" %7s%7s", "nson:", "sons");
   for (i = 0; i < tree.nnode; i++) {
      printf("\n%6d%6d", nodes[i].father, i);
      if (timebranches)
         printf(" %9.6f %9.6f %9.0f", nodes[i].age, nodes[i].branch, nodes[i].label);

      printf("%7d: ", nodes[i].nson);
      for (j = 0; j < nodes[i].nson; j++)  printf(" %2d", nodes[i].sons[j]);
      if (nodes[i].annotation) printf(" annotation: |%s| ", nodes[i].annotation);
   }
   FPN(F0);
   OutTreeN(F0, 0, 0); FPN(F0);
   OutTreeN(F0, 1, 0); FPN(F0);
   OutTreeN(F0, 1, 1); FPN(F0);
}


void PointconPnodes(void)
{
   /* This points the nodes[com.ns+inode].conP to the right space in com.conP.
      The space is different depending on com.cleandata (0 or 1)
      This routine updates internal nodes com.conP only.
      End nodes (com.conP0) are updated in InitConditionalPNode().
   */
   int nintern = 0, i;

   for (i = 0; i < tree.nbranch + 1; i++)
      if (nodes[i].nson > 0)  /* more thinking */
         nodes[i].conP = com.conP + (size_t)com.ncode*com.npatt*nintern++;
}


int SetxInitials(int np, double x[], double xb[][2])
{
   /* This forces initial values into the boundary of the space
    */
   int i;

   for (i = com.ntime; i < np; i++) {
      if (x[i] < xb[i][0] * 1.005) x[i] = xb[i][0] * 1.05;
      if (x[i] > xb[i][1] / 1.005) x[i] = xb[i][1] / 1.05;
   }
   for (i = 0; i < com.np; i++) {
      if (x[i] < xb[i][0]) x[i] = xb[i][0] * 1.2;
      if (x[i] > xb[i][1]) x[i] = xb[i][1] * .8;
   }
   return(0);
}


#if(defined(BASEML) || defined(CODEML) || defined(MCMCTREE))

int GetTipDate(double *TipDate, double *TipDate_TimeUnit)
{
   /* This scans species/sequence names to collect the sampling dates.  The last field of
      the name contains the date.
      Divergence times are rescaled by using TipDate_TimeUnit.
   */
   int i, j, indate, ndates = 0;
   double young = -1, old = -1;
   char *p;

   *TipDate = 0;
   for (i = 0, ndates = 0; i < stree.nspecies; i++) {
      stree.nodes[i].age = 0;
      j = (int)strlen(stree.nodes[i].name);
      for (indate = 0, p = stree.nodes[i].name + j - 1; j >= 0; j--, p--) {
         if (isdigit(*p) || *p == '.') indate = 1;
         else if (indate)
            break;
      }
      sscanf(p + 1, "%lf", &stree.nodes[i].age);
      if (stree.nodes[i].age <= 0)
         error2("Tip date <= 0: somehow i am using positive numbers only for tip date");
      else
         ndates++;

      if (i == 0)
         young = old = stree.nodes[i].age;
      else {
         old = min2(old, stree.nodes[i].age);
         young = max2(young, stree.nodes[i].age);
      }
   }
   if (ndates == 0) {
      if (*TipDate_TimeUnit == -1) *TipDate_TimeUnit = 1;
      return(0);
   }
   if (ndates != stree.nspecies)
      printf("TipDate model requires date for each sequence.");

   *TipDate = young;
   if (*TipDate_TimeUnit <= 0)
      *TipDate_TimeUnit = (young - old)*2.5;
   if (young - old < 1e-100)
      error2("TipDate: all sequences are of the same age?");
   for (i = 0; i < stree.nspecies * 2 - 1; i++) {
      if (i < stree.nspecies || stree.nodes[i].fossil) {
         stree.nodes[i].age = (young - stree.nodes[i].age) / *TipDate_TimeUnit;
         if (stree.nodes[i].age < 1e-100) stree.nodes[i].age = 0;
      }
   }
   com.tipDateOld = 0;
   for (i = 0; i < stree.nspecies * 2 - 1; i++) {
      if (i < stree.nspecies || stree.nodes[i].fossil) {
        com.tipDateOld = max2(com.tipDateOld, stree.nodes[i].age);
      }
   }


   if (noisy) printf("\nTipDate model\nDate range: (%.2f, %.2f) => (0, %.2f). TimeUnit = %.2f.\n",
      young, old, (young - old) / *TipDate_TimeUnit, *TipDate_TimeUnit);



   return(0);
}

#endif


int NSameBranch(char partition1[], char partition2[], int nib1, int nib2, int IBsame[])
{
   /* counts the number of correct (identical) bipartitions.
    nib1 and nib2 are the numbers of interior branches in the two trees
    correctIB[0,...,(correctbranch-1)] lists the correct interior branches,
    that is, interior branches in tree 1 that is also in tree 2.
    IBsame[i]=1 if interior branch i is correct.
    */
   int i, j, k = 0, nsamebranch;

#if(1)
   for (i = 0, nsamebranch = 0; i < nib1; i++)
      for (j = 0, IBsame[i] = 0; j < nib2; j++) {
         if (strcmp(partition1 + i*(com.ns + 1), partition2 + j*(com.ns + 1)) == 0) {
            nsamebranch++;  IBsame[i] = 1;  break;
         }
      }
#else
   for (i = 0, nsamebranch = 0; i < nib1; i++)
      for (j = 0, IBsame[i] = 0; j < nib2; j++) {
         for (k = 0; k < com.ns; k++)
            if (partition1[i*(com.ns + 1) + k] != partition2[j*(com.ns + 1) + k]) break;
         if (k == com.ns) {
            nsamebranch++;  IBsame[i] = 1;  break;
         }
      }
#endif
   return (nsamebranch);
}


int AddSpecies(int is, int ib)
{
   /* Add species (is) to tree at branch ib.  The tree currently has
    is+1-1 species.  Interior node numbers are increased by 2 to make
    room for the new nodes.
    if(com.clock && ib==tree.nbranch), the new species is added as an
    outgroup to the rooted tree.
    */
   int i, j, it;

   if (ib > tree.nbranch + 1 || (ib == tree.nbranch && !com.clock)) return(-1);

   if (ib == tree.nbranch && com.clock) {
      FOR(i, tree.nbranch) FOR(j, 2)
         if (tree.branches[i][j] >= is) tree.branches[i][j] += 2;
      it = tree.root;  if (tree.root >= is) it += 2;
      FOR(i, 2) tree.branches[tree.nbranch + i][0] = tree.root = is + 1;
      tree.branches[tree.nbranch++][1] = it;
      tree.branches[tree.nbranch++][1] = is;
   }
   else {
      FOR(i, tree.nbranch) FOR(j, 2)
         if (tree.branches[i][j] >= is) tree.branches[i][j] += 2;
      it = tree.branches[ib][1];
      tree.branches[ib][1] = is + 1;
      tree.branches[tree.nbranch][0] = is + 1;
      tree.branches[tree.nbranch++][1] = it;
      tree.branches[tree.nbranch][0] = is + 1;
      tree.branches[tree.nbranch++][1] = is;
      if (tree.root >= is) tree.root += 2;
   }
   BranchToNode();
   return (0);
}


#endif   /* ifdef NODESTRUCTURE */


#ifdef LFUNCTIONS

void InitializeNodeScale(void)
{
   /* This allocates memory to hold scale factors for nodes and also decide on the
    nodes for scaling by calling SetNodeScale().
    The scaling node is chosen before the iteration by counting the number of
    nodes visited in the post-order tree travesal algorithm (see the routine
    SetNodeScale).
    See Yang (2000 JME 51:423-432) for details.
    The memory required is  com.NnodeScale*com.npatt*sizeof(double).
    */
   int i, nS;

   if (com.clock >= 5) return;

   com.NnodeScale = 0;
   com.nodeScale = (char*)realloc(com.nodeScale, tree.nnode * sizeof(char));
   if (com.nodeScale == NULL) error2("oom");
   for (i = 0; i < tree.nnode; i++) com.nodeScale[i] = (char)0;
   SetNodeScale(tree.root);
   nS = com.NnodeScale*com.npatt;
   if (com.conPSiteClass) nS *= com.ncatG;
   if (com.NnodeScale) {
      if ((com.nodeScaleF = (double*)realloc(com.nodeScaleF, nS * sizeof(double))) == NULL)
         error2("oom nscale");
      memset(com.nodeScaleF, 0, nS * sizeof(double));

      if (noisy) {
         printf("\n%d node(s) used for scaling (Yang 2000 J Mol Evol 51:423-432):\n", com.NnodeScale);
         for (i = 0; i < tree.nnode; i++)
            if (com.nodeScale[i]) printf(" %2d", i + 1);
         FPN(F0);
      }
   }
}


int SetNodeScale(int inode)
{
   /* This marks nodes for applying scaling factors when calculating f[h].
    */
   int i, ison, d = 0, every;

   if (com.seqtype == 0)       every = 100;   /* baseml */
   else if (com.seqtype == 1)  every = 15;    /* codonml */
   else                     every = 50;    /* aaml */

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      d += (nodes[ison].nson ? SetNodeScale(ison) : 1);
   }
   if (inode != tree.root && d > every) {
      com.nodeScale[inode] = 1;
      d = 1;
      com.NnodeScale++;
   }
   return(d);
}

int NodeScale(int inode, int pos0, int pos1)
{
   /* scale to avoid underflow
    */
   int h, k, j, n = com.ncode;
   double t, smallw = 1e-12;

   for (j = 0, k = 0; j < tree.nnode; j++)   /* k-th node for scaling */
      if (j == inode) break;
      else if (com.nodeScale[j]) k++;

      for (h = pos0; h < pos1; h++) {
         for (j = 0, t = 0; j < n; j++)
            if (nodes[inode].conP[h*n + j] > t)
               t = nodes[inode].conP[h*n + j];

         if (t < 1e-300) {
            for (j = 0; j < n; j++)
               nodes[inode].conP[h*n + j] = 1;  /* both 0 and 1 fine */
            com.nodeScaleF[k*com.npatt + h] = -800;  /* this is problematic? */
         }
         else {
            for (j = 0; j < n; j++)
               nodes[inode].conP[h*n + j] /= t;
            com.nodeScaleF[k*com.npatt + h] = log(t);
         }
      }
      return(0);
}

static double *dfsites;



void print_lnf_site(int h, double logfh)
{
}

double lfundG(double x[], int np)
{
   /* likelihood function for site-class models.
      This deals with scaling for nodes to avoid underflow if(com.NnodeScale).
      The routine calls fx_r() to calculate com.fhK[], which holds log{f(x|r)}
      when scaling or f(x|r) when not.  Scaling factors are set and used for each
      site class (ir) to calculate log(f(x|r).  When scaling is used, the routine
      converts com.fhK[] into f(x|r), by collecting scaling factors into lnL.
      The rest of the calculation then becomes the same and relies on f(x|r).
      Check notes in fx_r.
      This is also used for NSsites models in codonml.
      Note that scaling is used between fx_r() and ConditionalPNode()
      When this routine is used under the multiple-gene site-class model, note
      that right now it assumes one set of com.freqK[] for the different genes,
      which may be an issue.
   */
   int h, ir, it, FPE = 0;
   double lnL = 0, fh = 0, t;

   NFunCall++;
   fx_r(x, np);

   for (h = 0; h < com.npatt; h++) {
      if (com.fpatt[h] <= 0 && com.print >= 0) continue;
      if (com.NnodeScale) { /* com.fhK[] has log{f(x|r}.  Note the scaling for nodes */
         for (ir = 1, it = 0; ir < com.ncatG; ir++) /* select term for scaling */
            if (com.fhK[ir*com.npatt + h] > com.fhK[it*com.npatt + h]) it = ir;
         t = com.fhK[it*com.npatt + h];
         for (ir = 0, fh = 0; ir < com.ncatG; ir++)
            fh += com.freqK[ir] * exp(com.fhK[ir*com.npatt + h] - t);
         fh = t + log(fh);
      }
      else {
         for (ir = 0, fh = 0; ir < com.ncatG; ir++)
            fh += com.freqK[ir] * com.fhK[ir*com.npatt + h];
         if (fh <= 0) {
            if (!FPE) {
               FPE = 1;  matout(F0, x, 1, np);
               printf("\nlfundG: h=%4d  fhK=%9.6e\ndata: ", h + 1, fh);
               print1site(F0, h);
               FPN(F0);
            }
            fh = 1e-300;
         }
         fh = log(fh);
      }
      lnL -= fh*com.fpatt[h];
      if (LASTROUND == 2) dfsites[h] = fh;
      if (com.print < 0) print_lnf_site(h, fh);
   }

   return(lnL);
}


int SetPSiteClass(int iclass, double x[])
{
   /* This sets parameters for the iclass-th site class
    This is used by ConditionalPNode() and also updateconP in both algorithms
    For method=0 and 1.
    */
   int k = com.nrgene + !com.fix_kappa;
   double *pkappa = NULL, *xcom = x + com.ntime, mr;

   _rateSite = com.rK[iclass];
   return (0);
}

extern int prt, Locus, Ir;


int fx_r(double x[], int np)
{
   /* This calculates f(x|r) if(com.NnodeScale==0) or log{f(x|r)}
    if(com.NnodeScale>0), that is, the (log) probability of observing data x
    at a site, given the rate r or dN/dS ratio for the site.  This is used by
    the discrete-gamma models in baseml and codeml as well as the NSsites models
    in codeml.
    The results are stored in com.fhK[com.ncatG*com.npatt].
    This deals with underflows with large trees using global variables
    com.nodeScale and com.nodeScaleF[com.NnodeScale*com.npatt].
    */
   int  h, ir, i, k, ig, FPE = 0;
   double fh, smallw = 1e-12; /* for testing site class with w=0 */


   for (ig = 0; ig < com.ngene; ig++) { /* alpha may differ over ig */
      if (com.Mgene > 1 || com.nalpha > 1)
         SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, x);
      for (ir = 0; ir < com.ncatG; ir++) {
         if (ir && com.conPSiteClass) {  /* shift com.nodeScaleF & conP */
            if (com.NnodeScale)
               com.nodeScaleF += (size_t)com.npatt*com.NnodeScale;
            for (i = com.ns; i < tree.nnode; i++)
               nodes[i].conP += (size_t)(tree.nnode - com.ns)*com.ncode*com.npatt;
         }
         SetPSiteClass(ir, x);
         ConditionalPNode(tree.root, ig, x);

         for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
            if (com.fpatt[h] <= 0 && com.print >= 0) continue;
            for (i = 0, fh = 0; i < com.ncode; i++)
               fh += com.pi[i] * nodes[tree.root].conP[h*com.ncode + i];
            if (fh <= 0) {
               if (fh < -1e-10 /* && !FPE */) { /* note that 0 may be o.k. here */
                  FPE = 1; matout(F0, x, 1, np);
                  printf("\nfx_r: h = %d  r = %d fhK = %.5e ", h + 1, ir + 1, fh);
                  if (com.seqtype == 0 || com.seqtype == 2) {
                     printf("Data: ");
                     print1site(F0, h);
                     FPN(F0);
                  }
               }
               fh = 1e-300;
            }
            if (!com.NnodeScale)
               com.fhK[ir*com.npatt + h] = fh;
            else {
               com.fhK[ir*com.npatt + h] = log(fh);
               for (k = 0; k < com.NnodeScale; k++)
                  com.fhK[ir*com.npatt + h] += com.nodeScaleF[k*com.npatt + h];
            }
         }  /* for (h) */
      }     /* for (ir) */

      if (com.conPSiteClass) {  /* shift pointers conP back */
         if (com.NnodeScale)
            com.nodeScaleF -= (com.ncatG - 1)*com.NnodeScale*(size_t)com.npatt;
         for (i = com.ns; i < tree.nnode; i++)
            nodes[i].conP -= (size_t)(com.ncatG - 1)*(tree.nnode - com.ns)*com.ncode*com.npatt;
      }
   }  /* for(ig) */
   return(0);
}


double lfun(double x[], int np)
{
   /* likelihood function for models of one rate for all sites including
    Mgene models.
    */
   int  h, i, k, ig, FPE = 0;
   double lnL = 0, fh;

   NFunCall++;

   for (ig = 0; ig < com.ngene; ig++) {
      if (com.Mgene > 1)
         SetPGene(ig, 1, 1, 0, x);
      ConditionalPNode(tree.root, ig, x);

      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
         if (com.fpatt[h] <= 0 && com.print >= 0) continue;
         for (i = 0, fh = 0; i < com.ncode; i++){
            fh += com.pi[i] * nodes[tree.root].conP[h*com.ncode + i];}
         if (fh <= 0) {
            if (fh < -1e-5 && noisy) {
               printf("\nfh = %.6f negative\n", fh);
               exit(-1);
            }
            if (!FPE) {
               FPE = 1;  matout(F0, x, 1, np);
               printf("lfun: h=%4d  fh=%9.6e\nData: ", h + 1, fh);
               print1site(F0, h);
               FPN(F0);
            }
            fh = 1e-80;
         }
         fh = log(fh);
         for (k = 0; k < com.NnodeScale; k++)
            fh += com.nodeScaleF[k*com.npatt + h];

         lnL -= fh*com.fpatt[h];
         if (LASTROUND == 2) dfsites[h] = fh;
         if (com.print < 0)
            print_lnf_site(h, fh);
      }
   }
   return (lnL);
}

#endif         /* #ifdef LFUNCTIONS */



#if (defined MCMCTREE)

int ProcessNodeAnnotation(int *haslabel)
{
   /* This processes fossil calibration information that has been read into
      nodes[].annotation, and copy the info into stree.nodes[].
      It uses both stree and nodes[], before nodes[] is destroyed.
      This is called before sequence alignments at loci are read.

      Further nodes:
      Other distributions such as G, SN, ST must be specified using the format 'G{alpha, beta}',
      say, and are processed here.  Simple bounds can also be specified using the format
      'L{0.5}', 'U{1.0}', or 'B{0.5, 1.0}', in which case they are processed here.
      I kept this complexity, as the notation of < and > is easy to understand, and because
      ReadTreeN can read node labels such as #, $ anyway.

      stree.duplication:
      stree.nodes[i].label = 0:  usual node, which may and may not have calibrations.
                            -1:  node i is the driving node, whose age is shared by other nodes.
      stree.nodes[j].label = i:  mirrored node j shares the age of node i, and can't have calibration.
   */
   int s = stree.nspecies, i, j, k, nfossiltype = 7;
   int *nodelabel, nmirror, mismatch;
   char *pch, *pch2, *braces="({";
   double tailL = 0.025, tailR = 0.025, p_LOWERBOUND = 0.1, c_LOWERBOUND = 1.0;

   nodelabel = (int*)malloc(s * sizeof(int));
   if (nodelabel == NULL) error2("oom");
   for (i = 0; i < s - 1; i++)  nodelabel[i] = stree.nodes[s + i].label = 0;

   for (i = s; i < s * 2 - 1; i++) {
      if (nodes[i].annotation == NULL)
         continue;
      if (pch = strchr(nodes[i].annotation, '#')) {    /* for models for duplicated genes */
         *haslabel = 1;
         sscanf(pch + 1, "%d", &stree.nodes[i].label);
         nodelabel[i - s] = stree.nodes[i].label;
      }

      if (strchr(nodes[i].annotation, '>') && strchr(nodes[i].annotation, '<')) {
         stree.nfossil++;
         stree.nodes[i].fossil = BOUND_F;
         stree.nodes[i].pfossil[2] = tailL;
         stree.nodes[i].pfossil[3] = tailR;
         if (pch = strchr(nodes[i].annotation, '>'))
            sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[0]);
         if (pch = strchr(nodes[i].annotation, '<'))
            sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[1]);
         if (stree.nodes[i].pfossil[0] > stree.nodes[i].pfossil[1]) {
            printf("fossil bounds (%.4f, %.4f)", stree.nodes[i].pfossil[0], stree.nodes[i].pfossil[1]);
            error2("fossil bounds in tree incorrect");
         }
      }
      else if (pch = strchr(nodes[i].annotation, '>')) {
         stree.nfossil++;
         stree.nodes[i].fossil = LOWER_F;
         sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[0]);
         stree.nodes[i].pfossil[1] = p_LOWERBOUND;
         stree.nodes[i].pfossil[2] = c_LOWERBOUND;
         stree.nodes[i].pfossil[3] = tailL;
      }
      else if (pch = strchr(nodes[i].annotation, '<')) {
         stree.nfossil++;
         stree.nodes[i].fossil = UPPER_F;
         sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[1]);
         stree.nodes[i].pfossil[2] = tailR;
      }
      else {
         for (j = 1; j < nfossiltype + 1; j++)
            if ((pch = strstr(nodes[i].annotation, fossils[j]))) break;
         if (j == nfossiltype + 1) continue;

         stree.nfossil++;
         stree.nodes[i].fossil = j;
         pch2 = strchr(pch, braces[0]);
         if (pch2 == NULL) pch2 = strchr(pch, braces[1]);
         if (pch2 == NULL)
           error2("did not find left brace");
         pch = pch2 + 1;

         switch (j) {
         case (LOWER_F):
            /* truncated Cauchy default prior L(tL, p, c) */
            stree.nodes[i].pfossil[1] = p_LOWERBOUND;
            stree.nodes[i].pfossil[2] = c_LOWERBOUND;
            stree.nodes[i].pfossil[3] = tailL;
            sscanf(pch, "%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1],
               &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3]);
            break;
         case (UPPER_F):
            stree.nodes[i].pfossil[2] = tailR;
            sscanf(pch, "%lf,%lf", &stree.nodes[i].pfossil[1], &stree.nodes[i].pfossil[2]);
            break;
         case (BOUND_F):
            stree.nodes[i].pfossil[2] = tailL;
            stree.nodes[i].pfossil[3] = tailR;
            sscanf(pch, "%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1],
               &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3]);
            if (stree.nodes[i].pfossil[0] > stree.nodes[i].pfossil[1]) {
               printf("fossil bounds (%.4f, %.4f)", stree.nodes[i].pfossil[0], stree.nodes[i].pfossil[1]);
               error2("fossil bounds in tree incorrect");
            }
            break;
         case (GAMMA_F):
            sscanf(pch, "%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1]);
            break;
         case (SKEWN_F):
            sscanf(pch, "%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1], &stree.nodes[i].pfossil[2]);
            break;
         case (SKEWT_F):
            sscanf(pch, "%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1], &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3]);
            break;
         case (S2N_F):
            sscanf(pch, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1],
               &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3], &stree.nodes[i].pfossil[4],
               &stree.nodes[i].pfossil[5], &stree.nodes[i].pfossil[6]);
            break;
         }
      }
      stree.nodes[i].usefossil = 1;
      nodes[i].branch = 0;
      free(nodes[i].annotation);
      nodes[i].annotation = NULL;
   }

   /* With duplication dating, fossil calibrations are copied to the main node, while mirrored nodes
      have no calibrations.  If i is the main node, mirrored node j has stree.nodes[j].label = i.
   */
   if (stree.duplication) {
      if (*haslabel == 0) error2("need node labels for duplication dating...");
      for (i = s; i < s * 2 - 1; i++) printf("%4d", nodelabel[i - s]);  FPN(F0);

      for (i = s; i < s * 2 - 1; i++) {
         if (nodelabel[i - s] == 0)  continue;
         nmirror = 0;
         for (j = i + 1; j < s * 2 - 1; j++) {
            if (nodelabel[i - s] != nodelabel[j - s]) continue;
            nmirror++;
            stree.duplication++;
            /* node j mirros i. */
            if (stree.nodes[i].fossil && stree.nodes[j].fossil) {
               /* if both & j have calibrations, they must be the same */
               mismatch = 0;
               if (stree.nodes[i].fossil != stree.nodes[j].fossil) mismatch = 1;
               for (k = 0; k < npfossils[stree.nodes[j].fossil]; k++)
                  if (fabs(stree.nodes[i].pfossil[k] - stree.nodes[j].pfossil[k]) > 1e-6) mismatch = 1;
               if (mismatch) {
                  printf("\nnodes %2d %2d share age but have different calibrations? ", i + 1, j + 1);
                  error2("err");
               }
               stree.nfossil--;
            }
            if (stree.nodes[i].fossil == 0 && stree.nodes[j].fossil) { /* copy fossil from j to i. */
               stree.nodes[i].usefossil = 1;
               stree.nodes[i].fossil = stree.nodes[j].fossil;
               for (k = 0; k < npfossils[stree.nodes[j].fossil]; k++) {
                  stree.nodes[i].pfossil[k] = stree.nodes[j].pfossil[k];
                  stree.nodes[j].pfossil[k] = -1;
               }
            }
            nodelabel[j - s] = 0;
            stree.nodes[j].usefossil = 0;
            stree.nodes[j].fossil = 0;
            stree.nodes[j].label = i;
            stree.nodes[i].label = -1;

            if (noisy) printf("node %2d is mirrored by node %d.\n", i, j);
         }
         if (nmirror == 0) {
            printf("\nnode %2d has label %d, but is not mirrored by any other node? ", i, nodelabel[i]);
            error2("error");
         }
      }

      for (i = s; i < s * 2 - 1; i++) printf("%4d", nodelabel[i - s]);  FPN(F0);
      for (i = s; i < s * 2 - 1; i++) printf("%4d", stree.nodes[i].label);  FPN(F0);

      /* stree.duplication has the number of mirrored nodes, so that s - 1 - stree.duplication is the
      number of distinct ages on the tree */
      stree.duplication--;
      if (noisy) printf("\nduplication dating: %2d node ages are mirrored\n", stree.duplication);
   }
   free(nodelabel);
   return(0);
}

#endif


/* routines for dating analysis of heterogeneous data */
#if (defined BASEML || defined CODEML || defined MCMCTREE)

int GenerateGtree(int locus);

int ReadTreeSeqs(FILE *fout)
{
   /* This reads the combined species tree, the fossil calibration information, and sequence
      data at each locus.  stree.nodes[].pfossil[] has tL, tU for bounds or alpha and beta
      for the gamma prior.

      This also constructs the gene tree at each locus, by pruning the master species tree.
   */
   FILE *fseq, *ftree;
   int haslength, i, j, locus, clean0 = com.cleandata;
   double tailL = 0.025, tailR = 0.025, p_LOWERBOUND = 0.1, c_LOWERBOUND = 1.0;

   ftree = gfopen(com.treef, "r");

   /* read master species tree and process fossil calibration info */
   fscanf(ftree, "%d%d", &stree.nspecies, &i);
   com.ns = stree.nspecies;
   if (com.ns > NS) error2("raise NS?");
   /* to read master species names into stree.nodes[].name */
   if (noisy) puts("Reading master tree.");
   for (j = 0; j < stree.nspecies; j++)
      com.spname[j] = stree.nodes[j].name;
   nodes = nodes_t;

   ReadTreeN(ftree, &haslength, 1, 1);
   OutTreeN(F0, 1, 0); FPN(F0);
   /* OutTreeN(F0,1,1); FPN(F0); */

   if (com.clock == 5 || com.clock == 6)
      for (i = 0; i < tree.nnode; i++)  nodes[i].branch = nodes[i].label = 0;
   for (i = 0; i < tree.nnode; i++)
      if (nodes[i].label < 0) nodes[i].label = 0;  /* change -1 into 0 */

   /* copy master tree into stree */
   if (tree.nnode != 2 * com.ns - 1)
      error2("check and think about multifurcating trees.");
   stree.nnode = tree.nnode;  stree.nbranch = tree.nbranch;
   stree.root = tree.root;    stree.nfossil = 0;
   for (i = 0; i < stree.nnode; i++) {
     stree.nodes[i].latent = 0;
      stree.nodes[i].father = nodes[i].father;
      stree.nodes[i].nson = nodes[i].nson;
      stree.nodes[i].label = -1;
      if (nodes[i].nson != 0 && nodes[i].nson != 2)
         error2("master tree has to be binary.");
      for (j = 0; j < stree.nodes[i].nson; j++)
         stree.nodes[i].sons[j] = nodes[i].sons[j];
   }

   /***** (((A1 , B1) '#1 B(0.1,0.3)', (A2 , B2) #1), C ) 'B(0.9, 1.1)';  *****/
#if (defined MCMCTREE)
   if (!com.TipDate)
      ProcessNodeAnnotation(&i);
#else
   ProcessNodeAnnotation(&i);
#endif

   /* read sequences at each locus, construct gene tree by pruning stree */
   data.ngene = com.ndata;
   com.ndata = 1;
   fseq = gfopen(com.seqf, "r");
   if ((gnodes = (struct TREEN**)malloc(sizeof(struct TREEN*)*data.ngene)) == NULL)
      error2("oom");

   printf("\nReading sequence data..  %d loci\n", data.ngene);
   for (locus = 0; locus < data.ngene; locus++) {
      fprintf(fout, "\n\n*** Locus %d ***\n", locus + 1);
      printf("\n\n*** Locus %d ***\n", locus + 1);

      com.cleandata = (char)clean0;
      for (j = 0; j < stree.nspecies; j++)
         com.spname[j] = NULL; /* points to nowhere */
      ReadSeq(fout, fseq, clean0, locus);               /* allocates com.spname[] */

      data.ns[locus] = com.ns;
      data.ls[locus] = com.ls;
#if(MCMCTREE)
      if (data.datatype[locus] == MORPHC)
         ;
      else
#endif
      {
         if (com.seqtype == 0 || com.seqtype == 2)
            InitializeBaseAA(fout);
         fflush(fout);

#if(!defined MCMCTREE)
         if ((com.seqtype == 0 || com.seqtype == 2) && com.model == 0) {
            PatternWeightJC69like();
            if (fout) {
               fprintf(fout, "\n\nPrinting out site pattern counts\n");
               printPatterns(fout);
            }
         }
#endif

         xtoy(com.pi, data.pi[locus], com.ncode);
         data.cleandata[locus] = (char)com.cleandata;
         data.npatt[locus] = com.npatt;
         data.fpatt[locus] = com.fpatt; com.fpatt = NULL;
         for (i = 0; i < com.ns; i++) {
            data.z[locus][i] = com.z[i];
            com.z[i] = NULL;
         }
         printf("%3d patterns, %s\n", com.npatt, (com.cleandata ? "clean" : "messy"));
      }

      GenerateGtree(locus);      /* free com.spname[] */
   }  /* for(locus) */
   if (noisy) printf("\n");
   for (i = 0, com.cleandata = 1; i < data.ngene; i++)
      if (data.cleandata[i] == 0)
         com.cleandata = 0;

   fclose(ftree); fclose(fseq);
   SetMapAmbiguity(com.seqtype, 0);

#if(defined (MCMCTREE))
   if (com.TipDate)
      GetTipDate(&com.TipDate, &com.TipDate_TimeUnit);

    if (com.latentData)
      ReadLatentSeqs();
#endif

   return(0);
}


int GenerateGtree(int locus)
{
   /* construct the gene tree at locus by pruning tips in the master species
    tree.  com.spname[] have names of species at the current locus (probably read
    from the sequence alignment at the locus).  They are used by the routine to compare
    with stree.nodes[].name to decide which species to keep for the locus.
    See GetSubTreeN() for more details.
    */
   int ns = data.ns[locus], i, j, ipop[NS], keep[NS], newnodeNO[2 * NS - 1];

   for (j = 0; j < stree.nspecies; j++) keep[j] = 0;
   for (i = 0; i < ns; i++) {
      for (j = 0; j < stree.nspecies; j++)
         if (!strcmp(com.spname[i], stree.nodes[j].name)) break;
      if (j == stree.nspecies) {
         printf("species %s not found in master tree\n", com.spname[i]);
         exit(-1);
      }
      if (keep[j]) {
         printf("\nspecies %s occurs twice in locus %d", com.spname[i], locus + 1);
         error2("\ngiving up...");
      }
      keep[j] = i + 1;  ipop[i] = j;  /* seq j in alignment is species i in master tree. */
      free(com.spname[i]);
   }

   /* copy master species tree and then prune it. */
   copySptree();
   GetSubTreeN(keep, newnodeNO);
   com.ns = ns;

   for (i = 0; i < stree.nnode; i++)
      if (newnodeNO[i] != -1) nodes[newnodeNO[i]].ipop = i;
   /* printGtree(0);  */

   gnodes[locus] = (struct TREEN*)malloc((ns * 2 - 1) * sizeof(struct TREEN));
   if (gnodes[locus] == NULL) error2("oom gtree");
   memcpy(gnodes[locus], nodes, (ns * 2 - 1) * sizeof(struct TREEN));
   data.root[locus] = tree.root;

   return(0);
}


int printGtree(int printBlength)
{
   int i, j;

   for (i = 0; i < com.ns; i++)
      com.spname[i] = stree.nodes[nodes[i].ipop].name;
   for (i = 0; i < tree.nnode; i++)
      if (i != tree.root)
         nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%7s%7s %8s %7s%7s", "father", "node", "(ipop)", "nson:", "sons");
   for (i = 0; i < tree.nnode; i++) {
      printf("\n%7d%7d   (%2d) %7d  ",
         nodes[i].father + 1, i + 1, nodes[i].ipop + 1, nodes[i].nson);
      for (j = 0; j < nodes[i].nson; j++) printf(" %2d", nodes[i].sons[j] + 1);
   }
   FPN(F0); OutTreeN(F0, 0, 0); FPN(F0); OutTreeN(F0, 1, 0); FPN(F0);
   if (printBlength) { OutTreeN(F0, 1, 1); FPN(F0); }
   return(0);
}


void copySptree(void)
{
   /* This copies stree into nodes = nodes_t, for printing or editing
    */
   int i, j;

   nodes = nodes_t;
   com.ns = stree.nspecies;   tree.root = stree.root;
   tree.nnode = stree.nnode;  tree.nbranch = stree.nbranch;
   for (i = 0; i < stree.nnode; i++) {
      /* this is used by mcmctree */
      if (i < com.ns) com.spname[i] = stree.nodes[i].name;

      /* The following may be needed by bpp.  Check carefully. */
      /*
       if(i<com.ns) strcpy(com.spname[i], stree.nodes[i].name);
       */
      nodes[i].father = stree.nodes[i].father;
      nodes[i].nson = stree.nodes[i].nson;
      for (j = 0; j < nodes[i].nson; j++)
         nodes[i].sons[j] = stree.nodes[i].sons[j];
      nodes[i].fossil = stree.nodes[i].fossil;
      nodes[i].age = stree.nodes[i].age;
      nodes[i].latent = stree.nodes[i].latent;
      nodes[i].latentAge = stree.nodes[i].latentAge;
      if (i != tree.root)
         nodes[i].branch = stree.nodes[nodes[i].father].age - stree.nodes[i].age;
   }
}

void printStree(void)
{
   int i, j, k, maxlname=10;

   for (i = 0; i < stree.nspecies; i++)
      if ((k = (int)strlen(stree.nodes[i].name)) > maxlname)
         maxlname = k;

   printf("\n************\nSpecies tree\nns = %d  nnode = %d", stree.nspecies, stree.nnode);
   printf("\n%7s%7s  %-*s %9s %8s%16s\n", "father", "node", maxlname+2, "name", "time", "sons", "fossil");
   for (i = 0; i < stree.nnode; i++) {
      printf("%7d%7d  %-*s %8.5f",
         stree.nodes[i].father + 1, i + 1, maxlname+6, stree.nodes[i].name, stree.nodes[i].age);
      if (stree.nodes[i].nson)
         printf("  (%2d %2d)", stree.nodes[i].sons[0] + 1, stree.nodes[i].sons[1] + 1);

#ifdef MCMCTREE
      if ((k = stree.nodes[i].fossil)) {
         printf("  %s ( ", fossils[k]);
         for (j = 0; j < npfossils[k]; j++) {
            printf("%6.4f", stree.nodes[i].pfossil[j + (k == UPPER_F)]);
            printf("%s", (j == npfossils[k] - 1 ? " ) " : ", "));
         }
      }
#endif

      printf("\n");
   }
   copySptree();
   FPN(F0); OutTreeN(F0, 0, 0); FPN(F0); OutTreeN(F0, 1, 0);  FPN(F0);
   OutTreeN(F0, 1, 1); FPN(F0);
}


#endif

void ReadLatentSeqs(void) {
  FILE *flatent;
  int ch = ' ';
  int lline = 32000;
  char line[32000];
  int numLat = 0;
  //memset(com.latent, 0, NS);
  flatent = gfopen(com.latentf, "r");

  regex_t regString, regWhiteSp;
  int regCompiled;

  /* Regex for word with possible whitespace on either side */
  regCompiled = regcomp(&regString, "^([[:space:]]*)([^[:space:]]+)([[:space:]]*)$" , REG_EXTENDED);

  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }

  size_t ngroups = regString.re_nsub + 1;
  regmatch_t *groups = malloc(ngroups * sizeof(regmatch_t));

  /* Regex for any amount of whitespace */
  regCompiled = regcomp(&regWhiteSp, "^([[:space:]]*)$" , REG_EXTENDED);

  if (regCompiled == 1) {
    fprintf(stderr, "Regular expression did not compile.\n");
    exit(1);
  }

  /* Loops through each line in the latent file */
  while(fgets(line, lline, flatent)) {

    /* If there is a single word on the line*/
    if (regexec(&regString, line, ngroups, groups, 0) == 0 ) {
      int i = 0;

      /* Find the species the name on the line matches */
      while (i < com.ns) {

        int length = max2(groups[2].rm_eo - groups[2].rm_so, strlen(stree.nodes[i].name));
        int strCompRet = strncmp(line + groups[2].rm_so, com.spname[i], length);

        if (strCompRet == 0) {
          /* Record which sequences are latent*/
          stree.nodes[i].latent = 1;
          stree.nodes[i].latentAge = stree.nodes[i].age;
          //printf("%s\n", line + groups[2].rm_so );
          break;
        }
        i++;
      }

      if (i == com.ns) {
        fprintf(stderr, "\nError: %s from latent file not found in species tree.\n", line);
        exit(-1);
      }
      numLat++;

    } else if (regexec(&regWhiteSp, line, 0, NULL, 0) == 0 ) ;

    else {
      error2("Latent file in the incorrect format. \nThere should be a single sequence name on each line.\n");
    }
  }

  /* Last line is in file goes through the loop twice.*/
  if (regexec(&regString, line, 0, NULL, 0))
    numLat--;

  if (numLat ==0 || numLat == com.ns) {
    error2("Number of latent sequences should be positive and less than the number of species. ");
  }

  regfree(&regString);
  regfree(&regWhiteSp);
  fclose(flatent);
}
