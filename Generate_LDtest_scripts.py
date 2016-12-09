import os

inds=[5,10,50]
sites=[100000,200000]
ploidy=[2,4]

for ind in inds:
	for site in sites:
		for p in ploidy:
			shfile1 = open('LD.sh','w')
			shfile1.write('#!/bin/bash\n'+
						'#SBATCH -J LD.sh'+'\n'+
						'#SBATCH -e LD.'+str(ind)+str(site)+str(p)+'.err\n'+
						'#SBATCH -o GS.'+str(ind)+str(site)+str(p)+'.out\n'+
						'#SBATCH -p nbi-short\n'+
						'#SBATCH -n 1\n'+
						'#SBATCH -t 2-5:00\n'+
						'#SBATCH --mem=32000\n'+
						'LD_LIBRARY_PATH=/nbi/software/testing/bin/core/../..//gcc/5.1.0/x86_64/lib64:$LD_LIBRARY_PATH\n'+
						'export LD_LIBRARY_PATH\n'+
						'/nbi/software/testing/ldr2/0.1/x86_64/LD /nbi/Research-Groups/JIC/Levi-Yant/Patrick/LDtestdata/LDTest_'+str(ind)+'ind_'+str(site/1000)+'Ksites_ploidy'+str(p)+'_randomfreqs.genotypes.txt /nbi/Research-Groups/JIC/Levi-Yant/Patrick/LDtestdata/LDTest_'+str(ind)+'ind_'+str(site/1000)+'Ksites_ploidy'+str(p)+'_randomfreqs.positions.txt ' + str(ind)+' '+str(site)+' 0.01 0.05 /nbi/Research-Groups/JIC/Levi-Yant/Patrick/LDtestdata/LDTest_'+str(ind)+'ind_'+str(site/1000)+'Ksites_ploidy'+str(p)+'_randomfreqs '+str(p)+'\n')
			shfile1.close()

		    # cmd1 = ('sbatch LD.sh')
		    # p1 = subprocess.Popen(cmd1, shell=True)
		    # sts1 = os.waitpid(p1.pid, 0)[1]

			file1 = open('LD.sh','r')
			data1 = file1.read()
			print(data1)
			os.remove('LD.sh')
