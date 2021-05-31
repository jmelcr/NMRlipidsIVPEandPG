DOI="10.5281/zenodo.1402449"

user_information = """
POPC:POPE 1:1 300K berger
#NMRLIPIDS BEGIN

@SIM
@MAPPING=POPC,mappingPOPCberger.txt,POPE,mappingPOPEberger.txt
@SYSTEM=POPE_T300K
@SOFTWARE=gromacs
@FF=berger
@FF_SOURCE=??
@FF_DATE=?/?/????
@TRJ=md_dt100_OK_centered.xtc
@TPR=md.tpr
@PREEQTIME=0
@TIMELEFTOUT=100
@UNITEDATOM=POPC,Berger_POPC,POPE,Berger_POPE

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


# Working directory
dir_wrk  = "/media/osollila/Data/tmp/DATABANK/"
