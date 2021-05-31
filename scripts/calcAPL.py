#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import sys
import re

import yaml
import json
import matplotlib.pyplot as plt
import mdtraj
import urllib.request
import seaborn as sns

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

def read_mapping_filePAIR(mapping_file, atom1, atom2, molname):
    with open(mapping_file, 'rt') as mapping_file:
            print(mapping_file)
            for line in mapping_file:
                if atom1 in line:
                    m_atom1 = line.split()[1]
                    try:
                        res = line.split()[2]
                    except:
                        res = molname
#                    print(m_atom1)
                if atom2 in line: 
                    m_atom2 = line.split()[1]
#                    print(m_atom2)
    return m_atom1, m_atom2, res

def make_positive_angles(x):
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = np.degrees(x[i]) + 360
        else:
            x[i] = np.degrees(x[i])
    return x


# ## Calculate the the angle between PN vector and z-axis  and plot the distribution of angles for each simulation

# In[3]:


#Calculate the angle between PN vector and z-axis for each lipid residues

lipids = {'POPS','POPE','POPG','POPC'}
#lipids = {'POPE'}
atom1 = 'M_G3P2_M'
atom2 = 'M_G3C5O1_M'   #for POPG
#atom2 = 'M_G3N6_M'

colors = {'POPC' :'black','POPS':'red','POPE':'blue','POPG':'green'}

h = []


for subdir, dirs, files in os.walk(r'../Data/Simulations/'):
    for filename in files:
        filepath = subdir + os.sep + filename
        if filepath.endswith("README.yaml"):
            READMEfilepath = subdir + '/README.yaml'
            with open(READMEfilepath) as yaml_file:
                readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                doi = readme.get('DOI')
                trj = readme.get('TRJ')
                tpr = readme.get('TPR')
                trj_name = subdir + '/' + readme.get('TRJ')[0][0]
                tpr_name = subdir + '/' + readme.get('TPR')[0][0]
                gro_name = subdir + '/conf.gro'
                trj_url = download_link(doi, trj[0][0])
                tpr_url = download_link(doi, tpr[0][0])
                EQtime=float(readme.get('TIMELEFTOUT'))*1000

                aplFOLDERS = subdir.replace("Simulations","APL")    
                outfilename = str(aplFOLDERS) + '/apl.dat' 

                #print(molname, readme['NPOPC'][0], outfilename)                    
                if not os.path.isfile(outfilename):
                                            
                    print('Analyzing '+filepath)
                        
                    #Download tpr and xtc files to same directory where dictionary and data are located
                    if (not os.path.isfile(tpr_name)):
                        response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                    if (not os.path.isfile(trj_name)):
                        response = urllib.request.urlretrieve(trj_url, trj_name)
                        
                    #fig= plt.figure(figsize=(12,9))
                    #if (not os.path.isfile(gro_name)):
                    #get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name}  -dump 0 -o {gro_name}')

                    get_ipython().system('mkdir -p {aplFOLDERS}')
                    get_ipython().system('cp {READMEfilepath} {aplFOLDERS}')

                    
                    box=aplFOLDERS + '/box.xvg'
                    if (not os.path.isfile(box)):
                        try:
                            get_ipython().system('echo System | gmx traj -f {trj_name} -s {tpr_name} -ob {box} -b {EQtime} -xvg none')
                        except:
                            print('Box calculation failed.')

                    try:
                        Nlipids = readme.get('NPOPC')[0] +  readme.get('NPOPC')[1] + readme.get('NPOPE')[0] +  readme.get('NPOPE')[1] + readme.get('NPOPG')[0] +  readme.get('NPOPG')[1]  +readme.get('NPOPS')[0] +  readme.get('NPOPS')[1] +readme.get('NCHOL')[0] +  readme.get('NCHOL')[1]
                    except:
                        Nlipids = readme.get('NPOPC')[0] +  readme.get('NPOPC')[1] + readme.get('NPOPE')[0] +  readme.get('NPOPE')[1] + readme.get('NPOPG')[0] +  readme.get('NPOPG')[1]  +readme.get('NPOPS')[0] +  readme.get('NPOPS')[1]

                    #os.system(awk -v Nlipids={Nlipids} '{areaSUM = areaSUM + ($2 * $3)}END{print areaSUM/Nlipids}' {box} > {aplFOLDERS}/apl.dat)
                    #os.system(awk '{areaSUM = areaSUM + ($2 * $3)}END{print areaSUM/100}' {box} > {aplFOLDERS}/apl.dat)

                    areaSUM = 0
                    lines = 0
                    with open(box) as boxDATA:
                        for line in boxDATA:
                            areaSUM = areaSUM + float(line.split()[1]) * float(line.split()[1])
                            lines = lines + 1
                        area = areaSUM / lines

                    outfile=open(outfilename,'w')
                        
                    #    for i in range(len(xaxis)):
                    #        #print(xaxis[i],distSUM[i])
                    outfile.write(str(2*area/Nlipids))
                    outfile.close()








