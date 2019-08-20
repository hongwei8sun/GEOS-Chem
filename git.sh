#!/bin/sh 
    git add .
    git status
    git commit -m " just delete all the gif before this comment "
#    git remote add origin https://github.com/hongwei8sun/GEOS-Chem.git
    git pull origin master
    git push origin master
