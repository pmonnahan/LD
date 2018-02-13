#########################
# Author: Patrick Monnahan
# Purpose:  This is a wrapper that submits jobs to the NBI slurm cluster that run LD.cpp for each population
#
#!!Must change library path based on the gcc version used in the compilation
#command for compilation is g++ -std=c++11 -Wall -fno-use-linker-plugin -o ../x86_64/LD LD.cpp
#and gcc version is 5.3.1

import os, sys, argparse, subprocess

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Takes input folder from LD_formatter_perPop.py and runs the LD.cpp program for each population present in the input folder.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='directory containing .positions and .genos files for each population as well as an INFO file containing the number of sites and individuals for each population')
parser.add_argument('-p', type=str, metavar='ploidy', required=True, help='2 for diploids; 4 for tets')
parser.add_argument('-mf', type=str, metavar='minimum-fraction-of-genotyped-individuals', default='0.9', help='Minimum fraction of individuals in a population that must be genotyped')
parser.add_argument('-maf', type=str, metavar='MinimumAlleleFrequency', default='0.05', help='minimum minor allele frequency for a site to be considered in a population')
parser.add_argument('-P', type=str, metavar='print', required=True, help='if true, then print shell scripts and not execute them')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Directory for output files')
parser.add_argument('-mem', type=str, metavar='memory', default='10000', help='number of Mb requested to run each job')

args = parser.parse_args()

pops=[]
numind=[]
numsites=[]


for file in os.listdir(args.i):
	if file.endswith('.tableINFO.txt'):
		infofile=open(args.i+file,'r')
		for line in infofile:
			line=line.strip("\n")
			line=line.split("\t")
			pops.append(line[0])
			numind.append(line[1])
			numsites.append(line[2])

if args.o.endswith("/") is False:
	args.o += "/"
if os.path.exists(args.o) is False:
    os.mkdir(args.o)
    os.mkdir(args.o + 'OandE/')
outputdir = args.o

infofile = open(outputdir + "RunParameters.txt",'w')
infofile.write("input directory = " + args.i + "\n" +
				"ploidy = " + args.p + "\n" +
				"minimum fraction of genotyped individuals = " + args.mf + "\n" +
				"minimum allele frequency = " + args.maf + "\n" +
				"output directory = " + args.o + "\n")


for j,pop in enumerate(pops):

	shfile = open('LD.sh','w')

	shfile.write('#!/bin/bash\n'+
					'#SBATCH -J LD.'+pop+'.sh'+'\n'+
					'#SBATCH -e '+outputdir+'OandE/LD.'+pop+'.err'+'\n'+
					'#SBATCH -o '+outputdir+'OandE/LD.'+pop+'.out'+'\n'+
					'#SBATCH -p nbi-medium\n'+
					'#SBATCH -n 1\n'+
					'#SBATCH -t 1-00:00\n'+
					'#SBATCH --mem='+args.mem+'\n'+
					'LD_LIBRARY_PATH=/nbi/software/testing/bin/core/../..//gcc/5.1.0/x86_64/lib64:$LD_LIBRARY_PATH\n'+
					'export LD_LIBRARY_PATH\n'+
					'/nbi/software/testing/ldr2/0.1/x86_64/LD '+args.i+pop+'.genotypes.txt '+args.i+pop+'.positions.txt ' + numind[j]+' '+numsites[j]+' '+args.mf+' '+args.maf+' '+outputdir+pop+' '+args.p+'\n')
	shfile.close()


	if args.P == 'false':
		cmd = ('sbatch LD.sh')
		p = subprocess.Popen(cmd, shell=True)
		sts = os.waitpid(p.pid, 0)[1]
	elif args.P == 'true':
		file = open('LD.sh','r')
		data = file.read()
		print(data)
	os.remove('LD.sh')
