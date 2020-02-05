#!/bin/sh 
    git add .
    git status
    git commit -m " seperate lagrangian and plume files before creating connections between them"
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
