SIMMHT + MLTA - Simulation of Magnetic Hyperthermia and Machine Learning Trainning Analysis

This repository includes three codes:

1 - ##SIMMHT## - a FORTRAN simulator that solves the time dependent behavior of the temperature field of a circular tumour subjected to magnetic hyperthermia;

2 - ##MLTA - Matlab version## - a matlab script that apply different machine learning regressor methods to predict the final temperature at the center of a spherical prototypical tumor and the time it takes to achive a steady-state condition based on the following inputs:
- nominal radius $a$ of the magnetic particles;
- the volume fraction $\phi$ of particles in the magnetic fluid;
- the amplitude $H_0$ of the applied field;
- the frequency $\omega$ of the applied field;
- 
3 - ##MLTA - Python version## - an equivalent python script that performs the same actions of the matlab version but using an open-source based philosophy that does not depend on proprietary software. This script uses de `sklearn` library.

Bellow we provide a description of each one of these codes.
  
