#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import mdtraj
import urllib.request
import seaborn as sns
import MDAnalysis as mda

from OrderParameter import *

# Download link
def download_link(doi, file):
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        return 'https://zenodo.org/record/' + zenodo_entry_number + '/files/' + file
    else:
        print ("DOI provided: {0}".format(doi))
        print ("Repository not validated. Please upload the data for example to zenodo.org")
        return ""
    
# read mapping file
def read_mapping_file(mapping_file, atom1):
    with open(mapping_file, 'rt') as mapping_file:
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
    return m_atom1

def read_mapping_fileRES(mapping_file, atom1):
    with open(mapping_file, 'rt') as mapping_file:
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
                    m_atom1res = line.split()[2]
    return m_atom1, m_atom1res

def read_mapping_filePAIR(mapping_file, atom1, atom2):
    with open(mapping_file, 'rt') as mapping_file:
            print(mapping_file)
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
#                    print(m_atom1)
                if atom2 in line: 
                    m_atom2 = line.split()[1]
#                    print(m_atom2)
    return m_atom1, m_atom2

def make_positive_angles(x):
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = np.degrees(x[i]) + 360
        else:
            x[i] = np.degrees(x[i])
    return x


def calcDihedrals(lipids,DIHatoms):
    colors = {'POPC' :'black','POPS':'red','POPE':'blue','POPG':'green'}
    for subdir, dirs, files in os.walk(r'../Data/Simulations/'):
        for filename in files:
            filepath = subdir + os.sep + filename
            if filepath.endswith("README.yaml"):
                READMEfilepath = subdir + '/README.yaml'
                with open(READMEfilepath) as yaml_file:
                    readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                    for molname in lipids:
                        doi = readme.get('DOI')
                        trj = readme.get('TRJ')
                        tpr = readme.get('TPR')
                        trj_name = subdir + '/' + readme.get('TRJ')[0][0]
                        tpr_name = subdir + '/' + readme.get('TPR')[0][0]
                        gro_name = subdir + '/conf.gro'
                        trj_url = download_link(doi, trj[0][0])
                        tpr_url = download_link(doi, tpr[0][0])
                        start_time = readme.get('TIMELEFTOUT')
    
                        #Download tpr and xtc files to same directory where dictionary and data are located
                        if (not os.path.isfile(tpr_name)):
                            response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                        if (not os.path.isfile(trj_name)):
                            response = urllib.request.urlretrieve(trj_url, trj_name)

                        dihedralFOLDERS = subdir.replace("Simulations","dihedral")
                        outfilename=str(dihedralFOLDERS) + '/' + molname + '_' + DIHatoms[0] + '_' + DIHatoms[1] + '_' + DIHatoms[2] + '_' + DIHatoms[3] +'.xvg'

                        print(outfilename)

                        get_ipython().system('mkdir -p {dihedralFOLDERS}')

                         
                        if not os.path.isfile(outfilename) and sum(readme['N' + molname]) > 0:
                            print('Analyzing '+molname+' in '+filepath)
                            #fig= plt.figure(figsize=(12,9))
                            if (not os.path.isfile(gro_name)):
                                get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name}  -dump 0 -o {gro_name}')
                        
                            xtcwhole=subdir + '/whole.xtc'
                            if (not os.path.isfile(xtcwhole)):
                                get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name} -o {xtcwhole} -b {start_time} -pbc mol ')
                        
                            #try:
                            #    traj = mdtraj.iterload(xtcwhole, chunk=1010, top = gro_name)
                            #    #print(traj)
                            #except FileNotFoundError or OSError:
                            #    continue

                            mapping_file = './mapping_files/'+readme['MAPPING_DICT'][molname] # readme.get('MAPPING')[0][0]
                            try:
                                atom1 = read_mapping_file(mapping_file, DIHatoms[0])
                                atom2 = read_mapping_file(mapping_file, DIHatoms[1])
                                atom3 = read_mapping_file(mapping_file, DIHatoms[2])
                                atom4 = read_mapping_file(mapping_file, DIHatoms[3])
                                print(atom1, atom2, atom3, atom4)
                            except:
                                #print(atom1 + " and " + atom2 + " not found in the mapping file.")
                                print("Some atom not found in the mapping file.")
                                continue
                            dihRESULT = []

                            ndx_name = dihedralFOLDERS + '/' + molname + '_' + atom1 +  '_' + atom2 + '_' + atom3 + '_' + atom4 + '.ndx'

                            try:
                                res1 = read_mapping_fileRES(mapping_file, DIHatoms[0])[1]
                                res2 = read_mapping_fileRES(mapping_file, DIHatoms[1])[1]
                                res3 = read_mapping_fileRES(mapping_file, DIHatoms[2])[1]
                                res4 = read_mapping_fileRES(mapping_file, DIHatoms[3])[1]
                            except:
                                res1 = readme.get(molname)
                                res2 = readme.get(molname)
                                res3 = readme.get(molname)
                                res4 = readme.get(molname)

                            #get_ipython().system('echo "r {res} & a {atom1} {atom2} {atom3} {atom4}  \n q \n " |  gmx make_ndx -f {tpr_name} -o {ndx_name}')

                            print(res1,atom1,res2,atom2,res3,atom3,res4,atom4)
                            
                            mol = mda.Universe(tpr_name, gro_name)
                            selection1 = mol.select_atoms("resname " + res1 + " and name " + atom1)
                            selection2 = mol.select_atoms("resname " + res2 + " and name " + atom2)
                            selection3 = mol.select_atoms("resname " + res3 + " and name " + atom3)
                            selection4 = mol.select_atoms("resname " + res4 + " and name " + atom4)
                            selection = []

                            if ((len(selection1) != len(selection2)) or (len(selection2) != len(selection3)) != (len(selection3) != len(selection4))):
                                continue
                            
                            for i in range(len(selection1)):
                                selection.append(selection1[i].id+1)
                                selection.append(selection2[i].id+1)
                                selection.append(selection3[i].id+1)
                                selection.append(selection4[i].id+1)

                            with open(ndx_name, 'w') as ndxfile:
                                ndxfile.write('[ '+ atom1 + '_' + atom2 + '_' + atom3 + '_' + atom4 + ' ] \n')
                                for item in selection:
                                    ndxfile.write(str(item) + " ")
                                    ndxfile.write(' \n')
                            
                            #get_ipython().system('echo "{res}_&_{atom1}_{atom2}_{atom3}_{atom4} \n" | gmx angle -f {xtcwhole} -od {outfilename} -type dihedral -n {ndx_name} -xvg none ')
                            get_ipython().system('gmx angle -f {xtcwhole} -od {outfilename} -type dihedral -n {ndx_name} -xvg none ')               
                            

                            get_ipython().system('cp {READMEfilepath} {dihedralFOLDERS}')
                            #outfile=open(str(dihedralFOLDERS) + '/' + DIHatoms[0] + DIHatoms[1] + DIHatoms[2] + DIHatoms[3] +'.dat','w')
                            #outfile=open(outfilename,'w') 
                            #for i in range(len(xaxis)):
                            #    outfile.write(str(xaxis[i]) + " " + str(dihRESULTav[i])+'\n')
                            #outfile.close()
                            #plt.close()
                            #del traj


dihedrals=[['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3O3_M'],
           ['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3P2O1_M'],
           ['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3P2O2_M'],
           ['M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M'],
           ['M_G3O1_M', 'M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M'],
           ['M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M'],
           ['M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M', 'M_G3N6_M'],
           ['M_G3O3_M','M_G3C4_M','M_G3C5_M','M_G3C6_M'],
           ['M_G3O3_M','M_G3C4_M','M_G3C5_M','M_G3C5O1_M'],
           ['M_G1C17_M', 'M_G1C16_M', 'M_G1C15_M', 'M_G1C14_M'],
           ['M_G1C16_M', 'M_G1C15_M', 'M_G1C14_M', 'M_G1C13_M'],
           ['M_G1C15_M', 'M_G1C14_M', 'M_G1C13_M', 'M_G1C12_M'],
           ['M_G1C14_M', 'M_G1C13_M', 'M_G1C12_M', 'M_G1C11_M'],
           ['M_G1C13_M', 'M_G1C12_M', 'M_G1C11_M', 'M_G1C10_M'],
           ['M_G1C12_M', 'M_G1C11_M', 'M_G1C10_M', 'M_G1C9_M'],
           ['M_G1C11_M', 'M_G1C10_M', 'M_G1C9_M', 'M_G1C8_M'],
           ['M_G1C10_M', 'M_G1C9_M', 'M_G1C8_M', 'M_G1C7_M'],
           ['M_G1C9_M', 'M_G1C8_M', 'M_G1C7_M', 'M_G1C6_M'],
           ['M_G1C8_M', 'M_G1C7_M', 'M_G1C6_M', 'M_G1C5_M'],
           ['M_G1C7_M', 'M_G1C6_M', 'M_G1C5_M', 'M_G1C4_M'],
           ['M_G1C6_M', 'M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M'],
           ['M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M', 'M_G1C2O1_M'],
           ['M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M', 'M_G1C2_M'],
           ['M_G1C4_M', 'M_G1C3_M', 'M_G1C2_M', 'M_G1O1_M'],
           ['M_G1C3_M', 'M_G1C2_M', 'M_G1O1_M', 'M_G1_M'],
           ['M_G2C19_M', 'M_G2C18_M', 'M_G2C17_M', 'M_G2C16_M'],
           ['M_G2C18_M', 'M_G2C17_M', 'M_G2C16_M', 'M_G2C15_M'],
           ['M_G2C17_M', 'M_G2C16_M', 'M_G2C15_M', 'M_G2C14_M'],
           ['M_G2C16_M', 'M_G2C15_M', 'M_G2C14_M', 'M_G2C13_M'],
           ['M_G2C15_M', 'M_G2C14_M', 'M_G2C13_M', 'M_G2C12_M'],
           ['M_G2C14_M', 'M_G2C13_M', 'M_G2C12_M', 'M_G2C11_M'],
           ['M_G2C13_M', 'M_G2C12_M', 'M_G2C11_M', 'M_G2C10_M'],
           ['M_G2C12_M', 'M_G2C11_M', 'M_G2C10_M', 'M_G2C9_M'],
           ['M_G2C11_M', 'M_G2C10_M', 'M_G2C9_M', 'M_G2C8_M'],
           ['M_G2C10_M', 'M_G2C9_M', 'M_G2C8_M', 'M_G2C7_M'],
           ['M_G2C9_M', 'M_G2C8_M', 'M_G2C7_M', 'M_G2C6_M'],
           ['M_G2C8_M', 'M_G2C7_M', 'M_G2C6_M', 'M_G2C5_M'],
           ['M_G2C7_M', 'M_G2C6_M', 'M_G2C5_M', 'M_G2C4_M'],
           ['M_G2C6_M', 'M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M'],
           ['M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M', 'M_G2C2O1_M'],
           ['M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M', 'M_G2C2_M'],
           ['M_G2C4_M', 'M_G2C3_M', 'M_G2C2_M', 'M_G2O1_M'],
           ['M_G2C3_M', 'M_G2C2_M', 'M_G2O1_M', 'M_G2_M'],
           ['M_G1O1_M', 'M_G1_M', 'M_G2_M', 'M_G3_M'],
           ['M_G1O1_M', 'M_G1_M', 'M_G2_M', 'M_G2O1_M'],
           ['M_G2C2_M', 'M_G2O1_M', 'M_G2_M', 'M_G1_M'],
           ['M_G1_M', 'M_G2_M', 'M_G3_M', 'M_G3O1_M'],
           ['M_G2O1_M', 'M_G2_M', 'M_G3_M', 'M_G3O1_M'],
           ['M_G2_M', 'M_G3_M', 'M_G3O1_M', 'M_G3P2_M'],
           ['M_G1C2O1_M', 'M_G1C2_M', 'M_G1O1_M', 'M_G1_M'],           
           ['M_G2C2O1_M', 'M_G2C2_M', 'M_G2O1_M', 'M_G2_M']
]

lipids = {'POPC','POPS','POPE','POPG'}
#lipids = {'POPC'}
#lipids = {'POPG','POPE','POPC'}
for i in dihedrals:
    DIHatoms = i
    calcDihedrals(lipids,i)

#mapping_file = "./mapping_files/mappingPOPEcharmm.txt"
#dihedrals = parseDihedralInput(mapping_file)
#parseDihedralInput(mapping_file)






