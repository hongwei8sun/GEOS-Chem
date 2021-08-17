!           $Id: aerosol.i,v 1.5 2012/12/14 21:04:09 dkweis Exp $
      PARAMETER (NSIZE=150,RZERO=3.9376E-4,Vrat=1.2,NHA=43,ICN=14)
!      PARAMETER (NSIZE=40,RZERO=3.9376E-4,Vrat=2.0,NHA=43,ICN=14)
!  AL2O3 number of particles, smallest size bin, density (gm/cc), 
!  fractal dimension for volume (Df) for area (Dh)
!  for sphere, Df=3, Dfh=3
!      PARAMETER (NSZAL=24,RZAL=0.01,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
!      PARAMETER (NSZAL=24,RZAL=0.01,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=2.8,Dfhwet=2.8) 
!      PARAMETER (NSZAL=18,RZAL=0.04,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
!      PARAMETER (NSZAL=18,RZAL=0.04,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=2.8,Dfhwet=2.8) 
!!      PARAMETER (NSZAL=18,RZAL=0.08,rhoal=2.7,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
!      PARAMETER (NSZAL=18,RZAL=0.08,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
!      PARAMETER (NSZAL=18,RZAL=0.08,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=2.8,Dfhwet=2.8) 
!      PARAMETER (NSZAL=10,RZAL=0.16,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
!      PARAMETER (NSZAL=10,RZAL=0.16,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=2.8,Dfhwet=2.8) 
!!      PARAMETER (NSZAL=8,RZAL=0.24,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
!      PARAMETER (NSZAL=8,RZAL=0.24,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=2.8,Dfhwet=2.8) 
!      PARAMETER (NSZAL=5,RZAL=0.32,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
!      PARAMETER (NSZAL=5,RZAL=0.32,rhoal=3.8,Df=1.6,Dfh=2.0,Dfwet=2.8,Dfhwet=2.8) 
!   for CaCO3 solid particles
!!     PARAMETER (NSZAL=8,RZAL=0.275,rhoal=2.71,almwt=100.09,Df=1.6,Dfh=2.0,Dfwet=1.6,Dfhwet=2.0) 
  PARAMETER (NSZAL=8,RZAL=0.275,rhoal=2.71,almwt=100.09,Df=2.1,Dfh=2.1,Dfwet=2.1,Dfhwet=2.1) 
      PARAMETER (rhoCaCO3=2.71,rhoCaCl2=2.15,rhoCaBr2=3.35,rhoCaN2O6=2.50,rhoCaSO4=2.32)
      PARAMETER (wmolCaCO3=100.09,wmolCaCl2=110.99,wmolCaBr2=199.89,wmolCaN2O6=164.09,wmolCaSO4=136.14)
