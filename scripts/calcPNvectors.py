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
                for molname in lipids:
                    doi = readme.get('DOI')
                    trj = readme.get('TRJ')
                    tpr = readme.get('TPR')
                    trj_name = subdir + '/' + readme.get('TRJ')[0][0]
                    tpr_name = subdir + '/' + readme.get('TPR')[0][0]
                    gro_name = subdir + '/conf.gro'
                    trj_url = download_link(doi, trj[0][0])
                    tpr_url = download_link(doi, tpr[0][0])
                    EQtime=float(readme.get('TIMELEFTOUT'))*1000

                    HGorientationFOLDERS = subdir.replace("Simulations","HGorientation")    
                    outfilename = str(HGorientationFOLDERS) + '/' + molname + 'PNvectorDIST.dat' 

                    print(molname, readme['NPOPC'][0], outfilename)                    
                    if not os.path.isfile(outfilename) and int(readme['N' + molname][0]) > 0:
                                            
                        print('Analyzing '+molname+' in '+filepath)
                        
                        #Download tpr and xtc files to same directory where dictionary and data are located
                        if (not os.path.isfile(tpr_name)):
                            response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                        if (not os.path.isfile(trj_name)):
                            response = urllib.request.urlretrieve(trj_url, trj_name)
                        
                        fig= plt.figure(figsize=(12,9))
                        if (not os.path.isfile(gro_name)):
                            get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name}  -dump 0 -o {gro_name}')
                        
                        xtcwhole=subdir + '/whole.xtc'
                        if (not os.path.isfile(xtcwhole)):
                            get_ipython().system('echo System | gmx trjconv -f {trj_name} -s {tpr_name} -o {xtcwhole} -pbc mol -b {EQtime}')
                        
                        mapping_file = './mapping_files/'+readme['MAPPING_DICT'][molname] # readme.get('MAPPING')[0][0]
                                                
                        try:
                            atoms = read_mapping_filePAIR(mapping_file, atom1 , atom2, readme[molname])
                        except:
                            print(atom1 + " and " + atom2 + " not found in the mapping file.")
                            continue

                        print(atoms, [atoms[0],atoms[1]])
                        
                        try:
                            anglesMeanError = read_trj_PN_angles(atoms[2], [atoms[0],atoms[1]] , tpr_name, xtcwhole, gro_name)
                        except OSError:
                            print("Could not calculate angles for " + molname + " from files " + tpr_name + " and " + trj_name)
                            continue
                        angles = anglesMeanError[0]
                        resAverageAngles = anglesMeanError[1]
                        totalAverage = anglesMeanError[2]
                        totalAverageError = anglesMeanError[3]
                        
                        #print(resAverageAngles)
                        
                        #fig= plt.figure(figsize=(12,9))
                        print(len(angles))
                        dist = [ 0 for i in range(len(angles))]
                        distSUM = [ 0 for i in range(0,178,2)]
                        # print(distSUM)
                        for i in range(len(angles)):
                            dist[i] = plt.hist(angles[i], bins=range(0,180,2), density=True)
                            distSUM = np.add(distSUM,dist[i][0])
                            #print(np.sum(dist[i][0]))
                            #print(np.sum(distSUM))
                            
                        # print(distSUM)
                        distSUM = [x / len(angles) for x in distSUM]
                        print(np.sum(distSUM))
                        
                        #print(len(dist))
                        xaxis = [ 0 for i in range(len(dist[0][1])-1)]
                        for i in range(len(dist[0][1])-1):
                            xaxis[i]=(dist[0][1][i])
                        
                        plt.legend()
                        plt.xlabel("Angle (°)")
                        #print(molname)
                        #print(colors[molname])
                        
                        
                        get_ipython().system('mkdir -p {HGorientationFOLDERS}')
                        get_ipython().system('cp {READMEfilepath} {HGorientationFOLDERS}')
                        outfile=open(outfilename,'w')
                        
                        for i in range(len(xaxis)):
                            #print(xaxis[i],distSUM[i])
                            outfile.write(str(xaxis[i]) + " " + str(distSUM[i])+'\n')
                        outfile.close()
                        outfile=open(str(HGorientationFOLDERS) + '/' + molname + 'PNvectorAVERAGE.dat','w')
                        outfile.write(str(totalAverage) + " " + str(totalAverageError)+'\n')
                        outfile.close()
            
#labels = [l.get_label() for l in h]
#plt.legend(h, labels, loc='best',ncol=2)
#plt.xlabel("Angle (°)")


# In[45]:


#print(totalAverage,totalAverageError)


# In[ ]:





# In[31]:







