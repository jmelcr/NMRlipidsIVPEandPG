#!/bin/bash
#Loop AddData.py over many simulations 

for f in *INFO.py
do
    ipython AddData.py $f
done
