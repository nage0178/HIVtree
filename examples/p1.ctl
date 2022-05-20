seed = 1

seqfile = p1_mcmc.fa 
treefile = p1_mcmc.txt
mcmcfile = mcmc.txt
outfile = out.txt

seqtype = 0
usedata = 1

ndata = 1
clock = 1
TipDate = 1 1000

RootAge = G(8,60)
model = 4
alpha = 1.5
ncatG = 5

cleandata = 0

alpha_gamma = 4 8 

rgene_gamma = 2 200

print = 1
burnin = 5000
sampfreq = 2
nsample = 70000
kappa_gamma = 8 1

latentFile = p1_seqsL
latentBound = 10.911
