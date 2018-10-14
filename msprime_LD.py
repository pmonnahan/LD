
import math
import subprocess
import msprime

def simulate(replicates_object, rep_id, num_ind):
    for j, tree_sequence in enumerate(replicates_object):
        dip_pos_out = open(f"Pos.Dip.MSP.{rep_id}.{j}.txt", 'w')
        dip_gen_out = open(f"Geno.Dip.MSP.{rep_id}.{j}.txt", 'w')
        tet_pos_out = open(f"Pos.Tet.MSP.{rep_id}.{j}.txt", 'w')
        tet_gen_out = open(f"Geno.Tet.MSP.{rep_id}.{j}.txt", 'w')
        for variant in tree_sequence.variants():
            dips = [sum(x) for x in chunks(variant.genotypes, 2)][0:num_ind+1]
            tets = [sum(x) for x in chunks(variant.genotypes, 4)][0:num_ind+1]
            if sum(dips) == 0 or sum(tets) == 0 or sum(dips) == (2 * num_ind) or sum(tets) == (4 * num_ind):
                pass
            else:
                site = int(math.ceil(variant.site.position))
                dip_pos_out.write(str(site) + "\n")
                tet_pos_out.write(str(site) + "\n")
                dip_gen_out.write("\t".join([str(x) for x in dips]) + "\n")
                tet_gen_out.write("\t".join([str(x) for x in tets]) + "\n")
        dip_pos_out.close()
        dip_gen_out.close()
        tet_pos_out.close()
        tet_gen_out.close()
    return()

def runLD(num_reps,rep_id, num_ind, mffg=0.1, minAF=0.05, LD_Path = "/Users/pmonnahan/Documents/Research/code/LD/LD"):
    for j in range(0,num_reps):    
        line_num = 0
        with open(f"Pos.Tet.MSP.{rep_id}.{j}.txt", 'r') as pos_file:
            for line in pos_file:
                line_num+=1
        cmd_t = f"{LD_Path} Geno.Tet.MSP.{rep_id}.{j}.txt Pos.Tet.MSP.{rep_id}.{j}.txt {num_ind} {line_num} {mffg} {minAF} LD.Tet.MSP.{rep_id}.{j}.af{minAF} 4"
        pp = subprocess.call(cmd_t.split())
        print(cmd_t)
        line_num = 0
        with open(f"Pos.Dip.MSP.{rep_id}.{j}.txt", 'r') as pos_file:
            for line in pos_file:
                line_num+=1
        cmd_d = f"{LD_Path} Geno.Dip.MSP.{rep_id}.{j}.txt Pos.Dip.MSP.{rep_id}.{j}.txt {num_ind} {line_num} {mffg} {minAF} LD.Dip.MSP.{rep_id}.{j}.af{minAF} 2"
        print(cmd_d)
        pp = subprocess.call(cmd_d.split())
    return()

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


num_ind = 10
num_reps = 100
Ne = 100000
hap_sampSize = 40
seed = 30
length = 15e5

rep_a = msprime.simulate(sample_size=hap_sampSize, Ne=Ne, length=length, mutation_rate=1e-8, random_seed=seed, recombination_rate=1e-8, num_replicates=num_reps)
rep_b = msprime.simulate(sample_size=hap_sampSize, Ne=Ne, length=length, mutation_rate=2e-8, random_seed=seed, recombination_rate=1e-8, num_replicates=num_reps)
rep_c = msprime.simulate(sample_size=hap_sampSize, Ne=Ne, length=length, mutation_rate=1e-8, random_seed=seed, recombination_rate=2e-8, num_replicates=num_reps)
rep_d = msprime.simulate(sample_size=hap_sampSize, Ne=Ne, length=length, mutation_rate=1e-8, random_seed=seed, recombination_rate=4e-8, num_replicates=num_reps)
rep_e = msprime.simulate(sample_size=hap_sampSize, Ne=Ne, length=length, mutation_rate=4e-8, random_seed=seed, recombination_rate=1e-8, num_replicates=num_reps)

prefix_ids = ["a","b","c","d","e"]
for i, rep in enumerate([rep_a, rep_b, rep_c, rep_d, rep_e]):
    simulate(rep, prefix_ids[i], num_ind)
    runLD(num_reps, prefix_ids[i], num_ind)





