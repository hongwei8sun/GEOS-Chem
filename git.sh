#!/bin/sh 
    git add .
    git status
    git commit -m " Correct: put lagrange_write in front of lagrange_run "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
