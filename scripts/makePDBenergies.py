import matplotlib.pyplot as plt
import math

def make_ene_distribution(PDBdata,SIMdata,outfile):
    PDBvalues=[]
    eneVALUES=[]
    DIHvalues={}
    with open(PDBdata,'r') as valuesFILE:
        lines = valuesFILE.readlines()
        for i in lines:
            PDBvalues.append(round(float(i.split()[0])))
    #values = valuesFILE.read()

    with open(SIMdata,'r') as dihFILE:
        lines = dihFILE.readlines()
        for i in lines:
            #print(i.split()[0],i.split()[1])
            DIHvalues[int(i.split()[0])]=float(i.split()[1])
            #dihedrals = dihFILE.read()

        for i in PDBvalues:
            if (i == 359):
                i = 358
            if (i == 360):
                i = 0
            if DIHvalues[i] > 0:
                eneVALUES.append((-1)*math.log(DIHvalues[i]))
            else:
                eneVALUES.append(10)

    eneMIN=min(eneVALUES)
    for i in range(len(eneVALUES)):
        eneVALUES[i]=eneVALUES[i]-eneMIN

    #print eneVALUES
    
    dist = plt.hist(eneVALUES, bins=range(0,12,1), density=True)
    with open(outfile,'w') as outfile:
        for i in range(len(dist[0])):
            outfile.write(str(dist[1][i]) + ' ' + str(dist[0][i]) + '\n')


files = [
    ["../Data/PC0values.dat","../Data/dihedral/b4f/866/b4f866c0dabffa6cd891e91841591d46590f34aa/a231c358dd8c0b7bb4cf558ccdbf022373354a9a/M_G3O3_MM_G3C4_MM_G3C5_MM_G3N6_M.dat","../Data/PC0eneDIST.dat"],
    ["../Data/PE0values.dat","../Data/dihedral/2d4/61b/2d461bc9d828af155146162ef42438974e4cbeaf/c8ce4cc36ef6ec7a085cf011176a92d9d746a792/M_G3O3_MM_G3C4_MM_G3C5_MM_G3N6_M.dat","../Data/PE0eneDIST.dat"],
    ["../Data/PG0values.dat","../Data/dihedral/0d5/d1d/0d5d1dcb43e775faf4e53c4f9ff255a67481bd38/9b487701b24d3fad83991e311188b08d3d5ea768/M_G3O3_MM_G3C4_MM_G3C5_MM_G3C6_M.dat","../Data/PG0eneDIST.dat"],
    ["../Data/PS0values.dat","../Data/dihedral/f40/bb6/f40bb6ab5d44402be07059e8df74b5a8200f031e/6774168dfec0a5a7377c8a46341eba603f320cf7/M_G3O3_MM_G3C4_MM_G3C5_MM_G3N6_M.dat","../Data/PS0eneDIST.dat"],

    ["../Data/PC1values.dat","../Data/dihedral/b4f/866/b4f866c0dabffa6cd891e91841591d46590f34aa/a231c358dd8c0b7bb4cf558ccdbf022373354a9a/M_G3P2_MM_G3O3_MM_G3C4_MM_G3C5_M.dat","../Data/PC1eneDIST.dat"],
    ["../Data/PE1values.dat","../Data/dihedral/2d4/61b/2d461bc9d828af155146162ef42438974e4cbeaf/c8ce4cc36ef6ec7a085cf011176a92d9d746a792/M_G3P2_MM_G3O3_MM_G3C4_MM_G3C5_M.dat","../Data/PE1eneDIST.dat"],
    ["../Data/PG1values.dat","../Data/dihedral/0d5/d1d/0d5d1dcb43e775faf4e53c4f9ff255a67481bd38/9b487701b24d3fad83991e311188b08d3d5ea768/M_G3P2_MM_G3O3_MM_G3C4_MM_G3C5_M.dat","../Data/PG1eneDIST.dat"],
    ["../Data/PS1values.dat","../Data/dihedral/f40/bb6/f40bb6ab5d44402be07059e8df74b5a8200f031e/6774168dfec0a5a7377c8a46341eba603f320cf7/M_G3P2_MM_G3O3_MM_G3C4_MM_G3C5_M.dat","../Data/PS1eneDIST.dat"],

    ["../Data/PC2values.dat","../Data/dihedral/b4f/866/b4f866c0dabffa6cd891e91841591d46590f34aa/a231c358dd8c0b7bb4cf558ccdbf022373354a9a/M_G3O1_MM_G3P2_MM_G3O3_MM_G3C4_M.dat","../Data/PC2eneDIST.dat"],
    ["../Data/PE2values.dat","../Data/dihedral/2d4/61b/2d461bc9d828af155146162ef42438974e4cbeaf/c8ce4cc36ef6ec7a085cf011176a92d9d746a792/M_G3O1_MM_G3P2_MM_G3O3_MM_G3C4_M.dat","../Data/PE2eneDIST.dat"],
    ["../Data/PG2values.dat","../Data/dihedral/0d5/d1d/0d5d1dcb43e775faf4e53c4f9ff255a67481bd38/9b487701b24d3fad83991e311188b08d3d5ea768/M_G3O1_MM_G3P2_MM_G3O3_MM_G3C4_M.dat","../Data/PG2eneDIST.dat"],
    ["../Data/PS3values.dat","../Data/dihedral/f40/bb6/f40bb6ab5d44402be07059e8df74b5a8200f031e/6774168dfec0a5a7377c8a46341eba603f320cf7/M_G3O1_MM_G3P2_MM_G3O3_MM_G3C4_M.dat","../Data/PS2eneDIST.dat"],

        ["../Data/PC3values.dat","../Data/dihedral/b4f/866/b4f866c0dabffa6cd891e91841591d46590f34aa/a231c358dd8c0b7bb4cf558ccdbf022373354a9a/M_G3_MM_G3O1_MM_G3P2_MM_G3O3_M.dat","../Data/PC3eneDIST.dat"],
    ["../Data/PE3values.dat","../Data/dihedral/2d4/61b/2d461bc9d828af155146162ef42438974e4cbeaf/c8ce4cc36ef6ec7a085cf011176a92d9d746a792/M_G3_MM_G3O1_MM_G3P2_MM_G3O3_M.dat","../Data/PE3eneDIST.dat"],
    ["../Data/PG3values.dat","../Data/dihedral/0d5/d1d/0d5d1dcb43e775faf4e53c4f9ff255a67481bd38/9b487701b24d3fad83991e311188b08d3d5ea768/M_G3_MM_G3O1_MM_G3P2_MM_G3O3_M.dat","../Data/PG3eneDIST.dat"],
    ["../Data/PS3values.dat","../Data/dihedral/f40/bb6/f40bb6ab5d44402be07059e8df74b5a8200f031e/6774168dfec0a5a7377c8a46341eba603f320cf7/M_G3_MM_G3O1_MM_G3P2_MM_G3O3_M.dat","../Data/PS3eneDIST.dat"],

        ["../Data/PC4values.dat","../Data/dihedral/b4f/866/b4f866c0dabffa6cd891e91841591d46590f34aa/a231c358dd8c0b7bb4cf558ccdbf022373354a9a/M_G2_MM_G3_MM_G3O1_MM_G3P2_M.dat","../Data/PC4eneDIST.dat"],
    ["../Data/PE4values.dat","../Data/dihedral/2d4/61b/2d461bc9d828af155146162ef42438974e4cbeaf/c8ce4cc36ef6ec7a085cf011176a92d9d746a792/M_G2_MM_G3_MM_G3O1_MM_G3P2_M.dat","../Data/PE4eneDIST.dat"],
    ["../Data/PG4values.dat","../Data/dihedral/0d5/d1d/0d5d1dcb43e775faf4e53c4f9ff255a67481bd38/9b487701b24d3fad83991e311188b08d3d5ea768/M_G2_MM_G3_MM_G3O1_MM_G3P2_M.dat","../Data/PG4eneDIST.dat"],
    ["../Data/PS4values.dat","../Data/dihedral/f40/bb6/f40bb6ab5d44402be07059e8df74b5a8200f031e/6774168dfec0a5a7377c8a46341eba603f320cf7/M_G2_MM_G3_MM_G3O1_MM_G3P2_M.dat","../Data/PS4eneDIST.dat"],

        ["../Data/PC5values.dat","../Data/dihedral/b4f/866/b4f866c0dabffa6cd891e91841591d46590f34aa/a231c358dd8c0b7bb4cf558ccdbf022373354a9a/M_G1_MM_G2_MM_G3_MM_G3O1_M.dat","../Data/PC5eneDIST.dat"],
    ["../Data/PE5values.dat","../Data/dihedral/2d4/61b/2d461bc9d828af155146162ef42438974e4cbeaf/c8ce4cc36ef6ec7a085cf011176a92d9d746a792/M_G1_MM_G2_MM_G3_MM_G3O1_M.dat","../Data/PE5eneDIST.dat"],
    ["../Data/PG5values.dat","../Data/dihedral/0d5/d1d/0d5d1dcb43e775faf4e53c4f9ff255a67481bd38/9b487701b24d3fad83991e311188b08d3d5ea768/M_G1_MM_G2_MM_G3_MM_G3O1_M.dat","../Data/PG5eneDIST.dat"],
    ["../Data/PS5values.dat","../Data/dihedral/f40/bb6/f40bb6ab5d44402be07059e8df74b5a8200f031e/6774168dfec0a5a7377c8a46341eba603f320cf7/M_G1_MM_G2_MM_G3_MM_G3O1_M.dat","../Data/PS5eneDIST.dat"],
 ]


for i in files:
    make_ene_distribution(i[0],i[1],i[2])


    #make_ene_distribution(files[i][0],files[0][1],files[0][2])
#make_ene_distribution("../Data/PC0values.dat","../Data/dihedral/b4f/866/b4f866c0dabffa6cd891e91841591d46590f34aa/a231c358dd8c0b7bb4cf558ccdbf022373354a9a/M_G3O3_MM_G3C4_MM_G3C5_MM_G3N6_M.dat","../Data/PC0eneDIST.dat")
        
