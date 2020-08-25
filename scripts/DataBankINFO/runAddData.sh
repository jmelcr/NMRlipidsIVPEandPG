#!/bin/bash
#Loop AddData.py over many simulations 

for f in POPCcharmmINFO.py  POPScharmm2INFO.py  POPScharmmINFO.py POPEcharmmINFO.py
do
    ipython3 AddData.py $f
done
