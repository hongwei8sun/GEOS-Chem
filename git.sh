#!/bin/sh 
    git add .
    git status
    git commit -m "success! pay attention to CFL condition"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
