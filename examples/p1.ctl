seed = 1

seqfile = p1_mcmc.fa 
treefile = p1_mcmc.txt
mcmcfile = mcmc.txt
outfile = out.txt

seqtype = 0 * 0 is nucleotide data
usedata = 1 * usedata = 1 produces the posterior distribution, usedate = 0 produces the prior distribution

ndata = 1
clock = 1 * Strict clock, required for HIVtree
TipDate = 1 1000

RootAge = G(8,60) * Root age prior
model = 4 * Substitution model, 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 1.5 * Gamma model of rate variation, turned off if 0, otherwise it will be on
ncatG = 5 * Number of rate categories in discrete gamma for the +Gamma model of rate variation

cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

alpha_gamma = 4 8 * Gamma prior on alpha parameter for +G model of rate variation
rgene_gamma = 2 200 * Gamma prior on the substitution rate, note this must be choosen in the time transformed time units
kappa_gamma = 8 1 * Gamma prior on kappa parameter in DNA substitution models with kappa

print = 1 * Determines amount of output of program
burnin = 5000 * Length of the burnin
sampfreq = 2 * Sample every other iteration
nsample = 70000 * Total number of samples after the burnin

latentFile = p1_seqsL * File containing the names of the latent sequences
latentBound = 10.911 * Oldest possible age of latent sequence in the transformed time units used in HIVtree
