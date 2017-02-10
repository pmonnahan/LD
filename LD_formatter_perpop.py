import os, sys, subprocess, argparse, re

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-t', type=str, metavar='table_dir', required=True, help='path to directory containing table file(s) (from -VariantsToTable) containing just scaff, pos, and GT fields')
parser.add_argument('-mf',type=float,metavar='min_fraction',required=True,default=0.75,help='minimum fraction of individuals with genotype calls in order for site to be included in output')
parser.add_argument('-p',type=float,metavar='Ploidy',required=True,help='2 for diploid; 4 for tetraploid')
parser.add_argument('-o',type=str,metavar='output',required=True,help='absolute path to directory for all population files')
parser.add_argument('-mip',type=int,metavar='MinimumIndividualsInPopulation',required=True,help='what is the minimum number of individuals in a population necessary to be used in LD code')


args = parser.parse_args()

pops = []
oldpop = ""
indkey = []
popcount = []
indcount = 0
sitecount = []

if os.path.exists(args.o) is False:
	os.mkdir(args.o)
outdir = args.o

table_files = []
for fil in os.listdir(args.t):
	if fil.endswith(".table"):
		table_files.append(fil)


for file in table_files:

	if os.path.exists(args.o + '/LDinput_' + file + '/') is False:
		os.mkdir(args.o + '/LDinput_' + file + '/')
	outdir = args.o + '/LDinput_' + file + '/'
	infofile=open(outdir + file + 'INFO.txt', 'w')

	with open(args.t + file,"rU") as table:
		prefix=args.t[:-6]
		# posFile=open(prefix+".positions.txt","w")
		# GTfile=open(prefix+".genotypes.txt","w")
		numsites=0
		for i,line in enumerate(table):
			if i==0:

				line=line.strip("\n")
				line=line.split("\t")
				numind=len(line[2:])
				scaf = line[0]

				for k,j in enumerate(line[2:]):
					curpop=j.split("_")[0]
					if k==0:
						oldpop=curpop
						pops.append(curpop)
						indkey.append(curpop)
						indcount=1
						exec "f%dgeno = open('%s%s.genotypes.txt','w')" % (len(pops),outdir,curpop)
						exec "f%dpos = open('%s%s.positions.txt','w')" % (len(pops),outdir,curpop)
					elif curpop != oldpop:
						popcount.append(indcount)
						sitecount.append(0)
						pops.append(curpop)
						indkey.append(curpop)
						indcount=1
						oldpop=curpop
						exec "f%dgeno = open('%s%s.genotypes.txt','w')" % (len(pops),outdir,curpop)
						exec "f%dpos = open('%s%s.positions.txt','w')" % (len(pops),outdir,curpop)
					elif curpop == oldpop:
						indcount+=1
						indkey.append(curpop)

				popcount.append(indcount)
				sitecount.append(0)
				exec "f%dgeno = open('%s%s.genotypes.txt','w')" % (len(pops),outdir,curpop)
				exec "f%dpos = open('%s%s.positions.txt','w')" % (len(pops),outdir,curpop)



			if i>0:
				line=line.strip("\n")
				line=line.split("\t")
				pos=line[1]
				

				for num,popp in enumerate(pops):
					ref='0'
					GT=''
					alt=False
					numobs=0
					if popcount[num]>=int(args.mip):
						for ind in range(0,popcount[num]):
							gt = line[2+sum(popcount[:num])+ind].split("/")
							if gt[0]=='.':
								GT+='9\t'
							else:
								gtc=int(args.p)
								if ref=='0':
									ref=gt[0]
								for g in gt:
									if g == ref:
										gtc-=1
									else:
										alt=True
								GT+=str(gtc)+"\t"
								numobs+=1

						GT.strip("\t")

						if alt==True and float(numobs)/float(popcount[num])>args.mf:
							exec "f%dgeno.write('%s')" % (num+1,GT)
							exec 'f%dgeno.write("""\n""")' % (num+1)
							exec "f%dpos.write('%s')" % (num+1,str(pos))
							exec 'f%dpos.write("""\n""")' % (num+1)
							sitecount[num]+=1


				if i%100000==0:
					print(i)	


	for h,g in enumerate(pops):
		st=g+"\t"+str(popcount[h])+"\t"+str(sitecount[h])+"\n"
		infofile.write(st)



