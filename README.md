# SIMMHT + MLTA 

## Simulation of Magnetic Hyperthermia and Machine Learning Trainning Analysis

This repository includes three codes:

**1 - SIMMHT** - a FORTRAN simulator that solves the time dependent behavior of the temperature field of a circular tumour subjected to magnetic hyperthermia;

**2 - MLTA - Matlab version** - a matlab script that apply different machine learning regressor methods to predict the final temperature at the center of a spherical prototypical tumor and the time it takes to achive a steady-state condition based on the following inputs:
- nominal radius $a$ of the magnetic particles;
- the volume fraction $\phi$ of particles in the magnetic fluid;
- the amplitude $H_0$ of the applied field;
- the frequency $\omega$ of the applied field;
 
**3 - MLTA - Python version** - an equivalent python script that performs the same actions of the matlab version but using an open-source based philosophy that does not depend on proprietary software. This script uses de `sklearn` library.

Bellow we provide a description of each one of these codes.
  

## SIMMHT

## MLTA - Matlab version

## MLTA - Python version

To install the scikit-learn (sklearn) package in the Spyder software, use the following steps:

- First install the miniconda anaconda distribution in your machine;
- Open the anaconda prompt from miniconda and create an environment with the command `conda create -n spyder-env -y`
- Activate the environment with `conda activate spyder-env`
- Then install the scikit-learn package: `pip install scikit-learn`
- Then, open the Spyder software and click in "Tools" => "Preferences" => "Python interpreter";
- Click on "use the following Python interpreter" and choose *spyder-env/python.ex*
- Click on apply. Then, click in *restart kernel* in the options button on the right corner;

The package can now be used.