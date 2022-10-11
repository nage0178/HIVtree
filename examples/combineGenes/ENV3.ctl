seed = 49

seqfile = ../257_ENV3.fa 
treefile = ../ENV3_mcmc.txt
mcmcfile = mcmc.txt
outfile = out.txt

seqtype = 0
usedata = 1

ndata = 1
clock = 1
TipDate = 1 1000

*RootAge = G(365,100)
RootAge = G(.25, 110)
model = 4
alpha = 1.5
ncatG = 5

cleandata = 0

alpha_gamma = 4 8 
rgene_gamma = 2 200
sigma2_gamma = 2 2 1

print = 1
burnin = 5000
sampfreq = 2
nsample = 50000
kappa_gamma = 8 1

latentFile = ../ENV3_seqsL
latentBound = 3.921
