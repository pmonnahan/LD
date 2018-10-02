
import math

replicates = msprime.simulate(sample_size=20, Ne=100000, length=20e6, mutation_rate=5e-7, random_seed=30, recombination_rate=5e-7, num_replicates=1000)



def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
num_ind = 10

replicates = msprime.simulate(sample_size=40, Ne=100000, length=15e5, mutation_rate=1e-8, random_seed=30, recombination_rate=1e-8, num_replicates=1000)


for j, tree_sequence in enumerate(replicates):
    dip_pos_out = open("Pos.Dip.MSP" + str(j) + ".txt", 'w')
    dip_gen_out = open("Geno.Dip.MSP" + str(j) + ".txt", 'w')
    tet_pos_out = open("Pos.Tet.MSP" + str(j) + ".txt", 'w')
    tet_gen_out = open("Geno.Tet.MSP" + str(j) + ".txt", 'w')
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




