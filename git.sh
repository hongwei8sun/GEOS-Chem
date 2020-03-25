#!/bin/sh 
    git add .
    git status
    git commit -m " last version of combining rings without adding new rings "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
