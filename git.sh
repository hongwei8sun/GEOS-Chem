#!/bin/sh 
    git add .
    git status
    git commit -m " plot gif for plume 1D 1mon "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
