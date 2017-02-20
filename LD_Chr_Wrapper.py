import os
import argparse
import subprocess

# create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Takes input folder from LD_formatter_perPop.py and runs the LD.cpp program for each population present in the input folder.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='directory containing .positions and .genos files for each population as well as an INFO file containing the number of sites and individuals for each population')
parser.add_argument('-p', type=str, metavar='ploidy', required=True, help='2 for diploids; 4 for tets')
parser.add_argument('-mf', type=str, metavar='minimum-fraction-of-genotyped-individuals', default='0.9', help='Minimum fraction of individuals in a population that must be genotyped')
parser.add_argument('-maf', type=str, metavar='MinimumAlleleFrequency', default='0.05', help='minimum minor allele frequency for a site to be considered in a population')
parser.add_argument('-ws', type=str, metavar='WindowSize', default='1000000', help='window size within which to calculate r2 between snps of given distance apart (specified by -d argument)')
parser.add_argument('-d', type=int, metavar='SNPdistance', default=10000, help='will calculate r2 for SNPs between d and d - 1000 bp apart')
parser.add_argument('-P', type=str, metavar='print', required=True, help='if true, then print shell scripts and not execute them')
parser.add_argument('-mem', type=str, metavar='memory', default='10000', help='number of Mb requested to run each job')

args = parser.parse_args()

pops = []
numind = []
numsites = []

dub = args.d
dlb = args.d - 1000

if args.i.endswith("/") is False:
    args.i += "/"

# !!Must change library path based on the gcc version used in the compilation
# command for compilation is g++ -std=c++11 -Wall -fno-use-linker-plugin -o ../x86_64/LD LD.cpp
# and gcc version is 5.1.0
for direct in os.listdir(args.i):
    if direct.endswith('.table'):
        for file in os.listdir(args.i + direct + "/"):
            if file.endswith('.tableINFO.txt'):
                infofile = open(args.i + direct + "/" + file, 'r')
                for line in infofile:
                    line = line.strip("\n")
                    line = line.split("\t")
                    pops.append(line[0])
                    numind.append(line[1])
                    numsites.append(line[2])
        inputdir = args.i + direct + "/"
        outputdir = args.i + direct + "/" + 'output_ws' + str(int(args.ws) / 1000) + "kb_d" + str(args.d) + "_mf" + args.mf + "_maf" + args.maf + "/"

        if args.P == 'false':
            if os.path.exists(outputdir) is False:
                os.mkdir(outputdir)
            if os.path.exists(outputdir + 'OandE/') is False:
                os.mkdir(outputdir + 'OandE/')

            infofile = open(outputdir + "RunParameters.txt", 'w')
            infofile.write("input directory = " + args.i + "\n" +
                           "ploidy = " + args.p + "\n" +
                           "minimum fraction of genotyped individuals = " + args.mf + "\n" +
                           "minimum allele frequency = " + args.maf + "\n" +
                           "output directory = " + args.o + "\n" +
                           "window size = " + args.ws + "\n" +
                           "r2 distance = " + str(args.d) + "\n")


        for j, pop in enumerate(pops):

            shfile = open('LD.sh', 'w')

            shfile.write('#!/bin/bash\n' +
                         '#SBATCH -J LD.' + pop + '.sh' + '\n' +
                         '#SBATCH -e ' + outputdir + 'OandE/LD_Chr.' + pop + '.err' + '\n' +
                         '#SBATCH -o ' + outputdir + 'OandE/LD_Chr.' + pop + '.out' + '\n' +
                         '#SBATCH -p nbi-short\n' +
                         '#SBATCH -n 1\n' +
                         '#SBATCH -t 0-01:00\n' +
                         '#SBATCH --mem=' + args.mem + '\n' +
                         'LD_LIBRARY_PATH=/nbi/software/testing/bin/core/../..//gcc/5.1.0/x86_64/lib64:$LD_LIBRARY_PATH\n' +
                         'export LD_LIBRARY_PATH\n' +
                         '/nbi/software/testing/ldr2/0.1/x86_64/LD_Chr ' + inputdir + pop + '.genotypes.txt ' + inputdir + pop + '.positions.txt ' + numind[j] + ' ' + numsites[j] + ' ' + args.mf + ' ' + args.maf + ' ' + outputdir + pop + ' ' + args.p + ' ' + args.ws + ' ' + str(dub) + ' ' + str(dlb) + '\n')
            shfile.close()


            if args.P == 'false':
                cmd = ('sbatch LD.sh')
                p = subprocess.Popen(cmd, shell=True)
                sts = os.waitpid(p.pid, 0)[1]
            elif args.P == 'true':
                file = open('LD.sh', 'r')
                data = file.read()
                print(data)
            os.remove('LD.sh')
