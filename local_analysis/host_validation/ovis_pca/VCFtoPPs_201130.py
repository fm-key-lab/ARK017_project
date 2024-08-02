#! /usr/bin/env python2

import re, sys, os, gzip
from operator import itemgetter
import numpy as np
arg = sys.argv

def column(mat, i):
    return [row[i] for row in mat]

def PLextract(vcfline, a_ref, a_alt):
    avec = ["AA", "AC", "CC", "GA", "GC", "GG", "TA", "TC", "TG", "TT"]
    if vcfline[-2] == "./.":
        return [0, 0, 0]
    else:
        tnum1 = [i for i,aval in enumerate(avec) if aval == a_ref*2]
        tnum1.extend([i for i,aval in enumerate(avec) if aval == a_ref+a_alt or aval == a_alt+a_ref])
        tnum1.extend([i for i,aval in enumerate(avec) if aval == a_alt*2])
        PLvec = [int(val) for val in vcfline[-1].split(",")]
        tPLs = [PLvec[t] for t in tnum1]
        minPL = min(tPLs)
        return [val - minPL for val in tPLs]


## A function to calculate GPs from PLs and priors
def GPcalc(PLvec, priorvec):
    tgl1s = [10**(-1*float(val)/10) for val in PLvec]  ## Calculate genotype likelihoods
    tgl2s = [val / sum(tgl1s) for val in tgl1s]  ## Normalize genotype likelihoods
    tgp1s = [tgl2 * priorvec[i] for i,tgl2 in enumerate(tgl2s)]  ## Calculate genotype probability
    ## Calculate genotype probability
    return [val / sum(tgp1s) for val in tgp1s]


argvec = [arg[i].split("=") for i in range(1,len(arg))]

inFile = [argval[1] for argval in argvec if argval[0] == "--in"][0]
posFile = [argval[1] for argval in argvec if argval[0] == "--pos"][0]
outFile = [argval[1] for argval in argvec if argval[0] == "--out"][0]

r1 = os.getcwd() + "/"; os.chdir(r1)

## Check if prior is provided; otherwise take the GATK default
priors = [0.4995, 0.0010, 0.4995]
if "--prior" in [val for val in column(argvec, 0)]:
    priors = [[float(v) for v in argval[1].split(",")] for argval in argvec if argval[0] == "--prior"][0]


## First, retrieve SNP information together with alleles
with open(posFile, "r") as F1:
    pvs = [[val for val in line.strip().split()] for line in F1]

F1.close()

## Then, import per-chromosome information from VCF
## This is already filtered to "CHROM POS REF ALT GT APL"
with open(inFile, "r") as F1:
    vcfs = [[val for val in line.strip().split()] for line in F1]

F1.close()

## Get Phred-scaled genotype likelihood information for all sites in the posFile
pvec2 = column(vcfs, 1)
count = 0; nmax = len(vcfs); nullvec = [0, 0, 0]
PLs = []
for snp in pvs:
    if count >= nmax:
        PLs.append([val for val in nullvec]); continue
    while snp[1] != pvec2[count]:
        count += 1
        if count >= nmax:
            break
    if count >= nmax:
        PLs.append([val for val in nullvec]); continue
    elif snp[1] == pvec2[count]:
        PLs.append([val for val in PLextract(vcfs[count], snp[-2], snp[-1])])
        continue
    else:
        PLs.append([val for val in nullvec]); continue


## Calculate GPs from PLs and priors
GPs = [GPcalc(PL, priorvec = priors) for PL in PLs]

## Write out PLs
F1 = open(outFile, "w")
F1.writelines('\n'.join([','.join(["{0:.4f}".format(val) for val in GP]) for GP in GPs]) + "\n")
F1.close()


