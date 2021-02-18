#!/bin/sh 
    git add .
    git status
    git commit -m " lagrange for GC_Class v13.0 "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
