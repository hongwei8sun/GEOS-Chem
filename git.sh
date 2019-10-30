#!/bin/sh 
    git add .
    git status
    git commit -m " Binary file for Lagrange output "
#    git commit -m "AGU: aircraft release aerosol boxes between (-30N,30N) for 1 month"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
