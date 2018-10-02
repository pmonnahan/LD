import subprocess

for i in range(0,999):
    LD_Path = "/Users/pmonnahan/Documents/Research/code/LD/LD"    
    line_num = 0
    with open(f"Pos.Tet.MSP{i}.txt", 'r') as pos_file:
        for line in pos_file:
            line_num+=1

    cmd_t = f"{LD_Path} Geno.Tet.MSP{i}.txt Pos.Tet.MSP{i}.txt 10 {line_num} 0.1 0.05 LD.Tet.MSP{i}.af05 4"
    pp = subprocess.call(cmd_t.split())
    print(cmd_t)
    line_num = 0
    with open(f"Pos.Dip.MSP{i}.txt", 'r') as pos_file:
        for line in pos_file:
            line_num+=1

    cmd_d = f"{LD_Path} Geno.Dip.MSP{i}.txt Pos.Dip.MSP{i}.txt 10 {line_num} 0.1 0.05 LD.Dip.MSP{i}.af05 2"
    print(cmd_d)
    pp = subprocess.call(cmd_d.split())
