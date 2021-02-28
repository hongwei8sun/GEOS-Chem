#!/bin/sh 
    git add .
    git status
    git commit -m " add one more inert tracer (for un-dissolved plume) in GEOS-Chem, not in plume model "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
