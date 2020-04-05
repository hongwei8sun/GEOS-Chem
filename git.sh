#!/bin/sh 
    git add .
    git status
    git commit -m " first version: dissolve PIG by comparing 2-order reaction rates "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
