#!/bin/sh 
    git add .
    git status
#    git commit -m " Binary file for Lagrange output "
    git commit -m "plot figure for [mol/cm3]; lon-Z; lat-Z"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
