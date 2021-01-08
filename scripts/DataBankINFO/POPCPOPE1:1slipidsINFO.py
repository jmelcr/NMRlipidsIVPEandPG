DOI="10.5281/zenodo.3605386"
user_information = """
POPC:POPE 1:1 slipids from M. Javanainen
#NMRLIPIDS BEGIN

@SIM
@SYSTEM=POPC:POPE_1:1_T298K
@MAPPING=POPC,mappingPOPCslipids.txt,POPE,mappingPOPEslipids.txt
@SOFTWARE=gromacs
@FF=Slipids
@FF_SOURCE=http://www.fos.su.se/~sasha/SLipids/
@FF_DATE=??
@TRJ=pcpe_slipids.xtc
@TPR=pcpe_slipids.tpr
@PREEQTIME=200
@TIMELEFTOUT=0

@POPC=POPC
@POPG=POPG
@POPS=POPS
@POPE=POPE

@POT=K
@SOD=NA
@CLA=CL
@CAL=CA
@SOL=SOL

@NPOPC=[0,0]
@NPOPG=[0,0]
@NPOPS=[0,0]
@NPOPE=[0,0]

@NPOT=0
@NSOD=0
@NCLA=0
@NCAL=0
@NSOL=0

@TEMPERATURE=0
@TRJLENGTH=0

#NMRLIPIDS END

"""
dir_wrk = "/media/osollila/Data/tmp/DATABANK/"
