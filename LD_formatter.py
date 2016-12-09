import os, sys, subprocess, argparse, re

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-T', type=str, metavar='table_file', required=True, help='path to directory containing table files (from -VariantsToTable) of scaffolds containing just scaff, pos, and GT fields')
parser.add_argument('-MF',type=float,metavar='min_fraction',required=True,default=0.75,help='minimum fraction of individuals with genotype calls in order for site to be included in output')
parser.add_argument('-P',type=float,metavar='Ploidy',required=True,help='2 for diploid; 4 for tetraploid')
args = parser.parse_args()


with open(args.T,"rU") as table:
	prefix=args.T[:-6]
	posFile=open(prefix+".positions.txt","w")
	GTfile=open(prefix+".genotypes.txt","w")
	numsites=0
	for i,line in enumerate(table):
		if i==0:
			line=line.strip("\n")
			line=line.split("\t")
			numind=len(line[2:])
		if i>0:
			ref='0'
			line=line.strip("\n")
			line=line.split("\t")
			pos=line[1]+"\n"
			GT=''
			alt=False
			numobs=0			
			for j,gt in enumerate(line[2:]):
				gt=gt.split("/")
				if gt[0]=='.':
					GT+='9\t'
				else:
					gtc=int(args.P)
					if ref=='0':
						ref=gt[0]
					for g in gt:
						if g == ref:
							gtc-=1
						else:
							alt=True
					GT+=str(gtc)+"\t"
					numobs+=1
			if i%100000==0:
				print i	
			GT.strip("\t")
			GT+="\n"
			if alt==True and float(numobs)/float(numind)>args.MF:
				numsites+=1
				GTfile.write(GT)
				posFile.write(pos)

	print "number of individuals = ", numind
	print "number of sites = ", numsites


