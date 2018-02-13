########################################
# Author: Patrick Monnahan
# Purpose: This code is meant to parse a table file created from GATK's VariantsToTable function that contains multiple individuals from multiple populations.  It will parse individuals by population and produce two files per population.  One file contains the positions of each site and the other file contains the genotypes for each individual in the population recoded as allele counts.  These inputs can then be used by LD.cpp or LD_Chrom.cpp to calculate average correlations among genotype values across sites
########################################




import os
import argparse

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-t', type=str, metavar='table_dir', required=True, help='path to directory containing table file(s) (from -VariantsToTable) containing just scaff, pos, and GT fields')
parser.add_argument('-mf', type=float, metavar='min_fraction', required=True, default=0.75, help='minimum fraction of individuals with genotype calls in order for site to be included in output')
parser.add_argument('-p', type=float, metavar='Ploidy', required=True, help='2 for diploid; 4 for tetraploid')
parser.add_argument('-o', type=str, metavar='output', required=True, help='absolute path to directory for all population files')
parser.add_argument('-mip', type=int, metavar='MinimumIndividualsInPopulation', required=True, help='what is the minimum number of individuals in a population necessary to be used in LD code')


args = parser.parse_args()

if os.path.exists(args.o) is False: # Create output directory
    os.mkdir(args.o)
outdir = args.o

table_files = []
for fil in os.listdir(args.t):  # Get list of table files
    if fil.endswith(".table"):
        table_files.append(fil)


for file in table_files:  # Cycle over table files.  These are likely separated by chromosome

    pops = []
    oldpop = ""
    indkey = []
    popcount = []
    indcount = 0
    sitecount = []

    if os.path.exists(args.o + '/LDinput_' + file + '/') is False:
        os.mkdir(args.o + '/LDinput_' + file + '/')
    outdir = args.o + '/LDinput_' + file + '/'
    infofile = open(outdir + file + 'INFO.txt', 'w') # Create file to hold the summary information for each population

    with open(args.t + file, "rU") as table:  # open the current table file
        prefix = args.t[:-6]
        numsites = 0
        for i, line in enumerate(table):
            if i == 0:  # Get list of populations and number of individuals from the header line of the table file.  

                line = line.strip("\n")
                line = line.split("\t")
                numind = len(line[2:])
                scaf = line[0]

                for k, j in enumerate(line[2:]):
                    curpop = j.split("_")[0]  # This expects population name to be first part of individual name followed by underscore (e.g. POP1_19)
                    if k == 0:
                        oldpop = curpop
                        pops.append(curpop)
                        indkey.append(curpop)
                        indcount = 1
                        exec "f%dgeno = open('%s%s.genotypes.txt','w')" % (len(pops), outdir, curpop) # creates genotypes file using ascertained population name
                        exec "f%dpos = open('%s%s.positions.txt','w')" % (len(pops), outdir, curpop)
                    elif curpop != oldpop: # Found all individuals from previous pop.  Write out info and move on to next population 
                        popcount.append(indcount)
                        sitecount.append(0)
                        pops.append(curpop)
                        indkey.append(curpop)
                        indcount = 1
                        oldpop = curpop
                        exec "f%dgeno = open('%s%s.genotypes.txt','w')" % (len(pops), outdir, curpop)
                        exec "f%dpos = open('%s%s.positions.txt','w')" % (len(pops), outdir, curpop)
                    elif curpop == oldpop:
                        indcount += 1
                        indkey.append(curpop)

                popcount.append(indcount)
                sitecount.append(0)
                exec "f%dgeno = open('%s%s.genotypes.txt','w')" % (len(pops), outdir, curpop)
                exec "f%dpos = open('%s%s.positions.txt','w')" % (len(pops), outdir, curpop)

            if i > 0: # Start parsing the actual data in the table file
                line = line.strip("\n")
                line = line.split("\t")
                pos = line[1]

                for num, popp in enumerate(pops):
                    ref = '0'
                    GT = ''
                    alt = False
                    numobs = 0
                    if popcount[num] >= int(args.mip): # Only report data for populations that meet criteria regarding minimum number of individuals in the population
                        for ind in range(0, popcount[num]):
                            gt = line[2 + sum(popcount[:num]) + ind].split("/")
                            if gt[0] == '.':
                                GT += '9\t'
                            else:
                                gtc = int(args.p)
                                if ref == '0':
                                    ref = gt[0]
                                for g in gt:
                                    if g == ref:
                                        gtc -= 1
                                    else:
                                        alt = True
                                GT += str(gtc) + "\t"
                                numobs += 1

                        GT.strip("\t")

                        if alt is True and float(numobs) / float(popcount[num]) > args.mf:
                            exec "f%dgeno.write('%s')" % (num + 1, GT)
                            exec 'f%dgeno.write("""\n""")' % (num + 1)
                            exec "f%dpos.write('%s')" % (num + 1, str(pos))
                            exec 'f%dpos.write("""\n""")' % (num + 1)
                            sitecount[num] += 1


                if i % 100000 == 0:
                    print(i)

    for h, g in enumerate(pops):
        st = g + "\t" + str(popcount[h]) + "\t" + str(sitecount[h]) + "\n"
        infofile.write(st)
