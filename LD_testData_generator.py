
import os, sys, subprocess, argparse, re
import numpy as np

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-ni', type=int, metavar='num_ind', required=True, help='Number of individuals')
parser.add_argument('-ns',type=int,metavar='num_sites',required=True,help='number of sites')
parser.add_argument('-p',type=int,metavar='Ploidy',required=True,help='2 for diploid; 4 for tetraploid')
parser.add_argument('-f',type=str,metavar='freq_model',required=True,help='rare, intermediate, or random')
parser.add_argument('-o',type=str,metavar='out_dir',required=True,help='output directory')
args = parser.parse_args()


exec "newgeno = open('%sLDTest_%dind_%dKsites_ploidy%d_%sfreqs.genotypes.txt','w')" % (args.o,args.ni,args.ns/1000,args.p,args.f)
exec "newpos = open('%sLDTest_%dind_%dKsites_ploidy%d_%sfreqs.positions.txt','w')" % (args.o,args.ni,args.ns/1000,args.p,args.f)

pos=0
for j in range(0,args.ns):
	if args.f == 'random':
		p = np.random.random()
	elif args.f == 'rare':
		rx=np.random.random()
		if rx<0.5:
			p = np.random.uniform(0.0,0.2)
		else:
			p = np.random.uniform(0.8,1.0)
	elif args.f == 'intermediate':
		p = np.random.uniform(0.3,0.7)
	else:
		print('invalid allele frequency model specified') 
		break

	pos+=10
	gg=np.random.binomial(args.p,p,args.ni)

	for i in gg:
		newgeno.write("%s\t" % str(i))
	newgeno.write("\n")
	newpos.write("%s\n" % str(pos))
