# Script for analyzing the conformations of lipid headgroups
# Uses python3, mdtraj and pymol
# Developed by Buslaev Pavel (pbuslaev@gmail.com)

import mdtraj as md
import numpy as np
import sys
import pymol

def load_traj(traj_file, traj_top, resname, sf, ef, stride):
	# Load topology file
	topol=md.load(traj_top).topology
	# Get atom indices to load
	ailist = topol.select('not water and not type W H Hs WT4 NaW KW CLW MgW and resname %s' \
	% resname)
	# load trajectory
	return md.load(traj_file,top = traj_top,atom_indices=ailist,stride=stride)[sf:ef]

def concat_traj(traj,show_atoms):
	# Get number of frames and lipids
	nframes = traj.n_frames
	nlip = len(list(traj.topology.residues))
	
	# Get topologe
	top = traj.topology

	# Get ids of atoms to show
	show_ai = top.select(show_atoms)

	# Reshape the trajectory
	trajxyz=traj.atom_slice(show_ai).xyz
	trajxyz=trajxyz.reshape((nframes,nlip,len(show_ai)//nlip,3))
	trajxyz = np.swapaxes(trajxyz,0,1)
	trajxyz = trajxyz.reshape((nframes*nlip,len(show_ai)//nlip,3))

	return md.Trajectory(trajxyz,traj[0].atom_slice(show_ai[0:len(show_ai)//nlip]).topology)

def align_traj(traj, ref_atoms):
	trajxyz = traj.xyz
	top = traj.topology
	ref_ai = top.select("name " + ref_atoms)

	trajxyz = centerTraj(trajxyz, ref_ai)
	trajxyz = rotateTraj(trajxyz, trajxyz[0], ref_ai)
	avg_xyz = trajxyz.mean(axis=0)
	trajxyz = centerTraj(trajxyz, ref_ai)
	trajxyz = rotateTraj(trajxyz,avg_xyz, ref_ai)

	return md.Trajectory(trajxyz,top)

def save_traj(traj, n_out_f,outf):
	nf = traj.n_frames
	st = (nf // n_out_f) + 1
	traj[::st].save_pdb(outf)

# Kabsch algorithm for alignment
def kabsch(p,q): # p - moving vector; q - reference 
    # Computation of covariance matrix
    C = np.dot(np.transpose(p),q)
    # Singular value decomposition
    V,S,W = np.linalg.svd(C)
    d = (np.linalg.det(V)*np.linalg.det(W)) < 0.
    if d:
        S[-1]= -S[-1]
        V[:, -1] = -V[:, -1]

    # Create rotation matrix
    U = np.dot(V,W)
    return U

# Moving selections to the same center
def centerTraj(trajxyz,align_i):
    xyz = trajxyz
    # Get number of frames, atoms and coordinate
    (nf,na,nc) = xyz.shape
    # Find the center of mass for each frame
    list_of_cm = np.mean(xyz[:,align_i,:], axis=1)
    # Center the trajectory
    xyz = np.array(list(map(lambda x,y: x - y, xyz,list_of_cm)))
    return xyz

# Rotating selections to minimize rmsd
def rotateTraj(trajxyz,refxyz,align_i):
    xyz = trajxyz
    # Find center of mass for the reference
    rcm = np.mean(refxyz,axis=0)
    # Center referencd
    ref = refxyz - rcm
    # Get number of frames, atoms and coordinate
    (nf,na,nc) = xyz.shape
    xyz = np.array(list(map(lambda x: np.dot(x,kabsch(x[align_i],ref[align_i])),xyz)))
    return xyz

# Calculate dihedral angle distributions
def parseDihedFile(file,traj):
	# Get topology
	top = traj.topology
	diheds = []
	with open(file) as f:
		for line in f:
			# for each line from the input file select atoms and 
			# calculate dihed for the trajectory
			ai = top.select("name " + line)
			ai = np.array(np.split(ai,ai.shape[0]//4))
			dihed = np.concatenate(md.compute_dihedrals(traj,ai),axis=None)
			# Append calculated dihedrals to the diheds
			diheds.append(dihed)
	diheds=np.array(diheds)
	np.savetxt("diheds.txt",diheds)

# Show dihedrals
def showDihed(ref,structure,a_atoms,pngf):
	topol =md.load(ref).topology
	satoms = list(i.split('-')[1] for i in str(list(topol.atoms))[1:-1].split(","))
	
	refatoms = {\
	"hg":["C11","N","P","C12","O12"],\
	"gb":["C1","C2","C3","P","O11","C21","O21","C31","O31"]
	}
	
	# views for different sets of atoms
	# hg - P to N atoms for different lipids
	# gb - glycerol backbone
	views={\
	"hg":"\
    -0.816361189,   -0.477001756,   -0.325615883,\
    -0.414168328,    0.876450181,   -0.245559067,\
     0.402518690,   -0.065604433,   -0.913057089,\
    -0.000000283,    0.000000288,  -14.178191185,\
     0.205995172,   -0.187574118,    0.333097398,\
    11.678203583,   16.678203583,   20.000000000 ",
    "gb":"\
     0.972944081,    0.089319512,   -0.213074163,\
    -0.202513024,   -0.114210837,   -0.972596884,\
    -0.111205898,    0.989433169,   -0.093032047,\
     0.000000296,   -0.000000288,  -18.411260605,\
    -0.077071443,    0.284420639,    0.109496683,\
    14.169631958,   21.652887344,   20.000000000 "}

	viewk = ""

	for k,v in refatoms.items():
		if all(elem in satoms for elem in v):
			viewk = k
			break

	pymol.cmd.reinitialize()

	#pymol.finish_launching()

	pymol.cmd.bg_color("white")
	pymol.cmd.set("orthoscopic","true")
	pymol.cmd.set("ray_shadow","off")

	pymol.cmd.viewport(600,400)

	pymol.cmd.load(ref,"ref")
	pymol.cmd.load(structure,"st")

	pymol.cmd.align("st and "+a_atoms,"ref "+a_atoms)

	pymol.cmd.split_states("st")

	pymol.cmd.delete("st")
	pymol.cmd.delete("ref")

	pymol.cmd.set_view (views[viewk])

	pymol.cmd.hide("everything")
	pymol.cmd.show("sticks")
	pymol.cmd.set("stick_radius",0.05)
	pymol.util.cbag()

	pymol.cmd.ray()
	pymol.cmd.png(pngf)
	pymol.cmd.quit()

class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        return self.value != None
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        elif isinstance(v,str):
        	self.value = self.func(v)
        else:
        	self.value = [ self.func(i) for i in v ]

def main(args):

	options = [
	# options for concat feature
	"Input/output options",
	("-f", Option(str, 1, None, "Input trajectory file (.xtc, .trr, ...)")),
	("-t", Option(str, 1, None, "Input topology file (.pdb, .gro, ...)")),
	("-st", Option(int, 1, 1, "Only read every Nth frame")),
	("-sf", Option(int, 1, 0, "First frame (ps) to read from trajectory")),
	("-ef", Option(int, 1, -1, "Last frame (ps) to read from trajectory")),
	("-sn", Option(str, -1, None, "Atoms to show")),
	("-an", Option(str, -1, None, "Atoms for alignment")),
	("-nf", Option(int, 1, 100, "Number of frames in the output")),
	("-l", Option(str, 1, None, "Lipid type = resname")),
	("-o", Option(str, 1, "out.pdb", "Output frames (.pdb)")),
	("-di", Option(str, 1, None, "List of dihedrals")),
	("-r", Option(str, 1, None, "Reference structure (.pdb)")),
	("-so", Option(str, 1, "out.png", "Dihedral structures output figure"))		
	]

	# if the user asks for help: pcalipids.py concat -h
	if (len(args)>0 and (args[0] == '-h' or args[0] == '--help')) or len(args)==0:
		print("\n",__file__[__file__.rfind('/')+1:])
		for thing in options: # print all options for selected feature
			print(type(thing) != str and "%10s: %s"%(thing[0],thing[1].description) or thing)
		print()
		sys.exit()

	options = dict([i for i in options if not type(i) == str])

	print(args)
	while args:
		ar = args.pop(0) # choose argument
		if options[ar].num == -1:
			listOfInputs = ""
			while args:
				ar1 = args.pop(0)
				if ar1 in list(options.keys()):
					options[ar].setvalue(listOfInputs)
					args.insert(0,ar1)
					break
				else:
					listOfInputs += (ar1+" ")
			options[ar].setvalue(listOfInputs)
		else:
			options[ar].setvalue([args.pop(0) for i in range(options[ar].num)]) # set value

	if not options["-f"].value or not options["-t"].value or not options["-l"].value or not options["-r"].value:
		print("Trajectory, structure, reference structure and the lipid resname have to be provided")
		return
	
	print(options["-sn"].value)

	if not options["-sn"].value:
		show_atoms = "all"
	else:
		show_atoms = "name "+options["-sn"].value

	if not options["-an"].value:
		ref_atoms = "all"
		ref_atoms_pymol = "name *"
	else:
		ref_atoms = "name "+options["-an"].value
		ref_atoms_pymol = "name "+(options["-an"].value).replace(" ","+")

	traj = load_traj(options["-f"].value,\
					 options["-t"].value,\
					 options["-l"].value,\
					 options["-sf"].value,\
					 options["-ef"].value,\
					 options["-st"].value)

	if options["-di"].value:
		parseDihedFile(options["-di"].value, traj)

	ctraj = concat_traj(traj, show_atoms)
	atraj = align_traj(ctraj,ref_atoms)

	save_traj(atraj,options["-nf"].value,options["-o"].value)

	print("pymol input")

	showDihed(options["-r"].value,options["-o"].value,ref_atoms_pymol,options["-so"].value)

if __name__ == '__main__':
	args = sys.argv[1:]
	main(args)
