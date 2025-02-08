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

This code was writen in FORTRAN and solves the time-dependent temperature field of a circular tumor subjected to magnetic hyperthermia (MHT) in a 2D finite differences mesh. The time evolution is solved explicitly. The governing equation solved is a non-linear second order PDE that incorporates the following physical mechanisms:

- Heat diffusion;
- Metabolic generation;
- Heat tranfer through the blood perfusion;
- Magnetic production term due to MHT;

This code has been used before in other works of the group **[1]** and has been validated using an experimental in-vivo work **[2]**. The source code is organized in several files that are compiled using `make`. Therefore, in order to produce an executable file (program) based on the source code you can simply run the `make` command on a Linux terminal. The code has been tested based on compilations performed by the intel fortran compiler `ifort`, but it should also work with `gfortran`. In order to change the compilation instructions you should alter the `makefile` included in this repository inside the folder **SIMMHT**.

The source code files are:

- *config.dat*: configuration file where the user can alter the simulation data. This file is not necessary for compilation purposes, but the code was designed to read it in this specific format;
  
- *input.f90*: FORTRAN file that reads the configuration file and stores all information in terms of variables used through the calculations performed;
  
- *variables.f90*: FORTRAN module with the declaration of all global variables used in the code;
  
- *simmht.f90*: FORTRAN file that calls the necessary functions in a given order to perform the simulations that produce a database;
  
- *functions.f90*: FORTRAN module with all the subroutines used in the calculations;

In order to run the simulations, the user needs to follow these steps:

1 - Compile the solver by running the command: `make`. This shall produce an executable file named **simmht.ex**;

2 - Run the simulations by typing: `./simmht.ex`. This procedure will generate two new files and at the same time will print on the terminal the time and the temperature at the center of the each realization. A given realization ends when the tumor reachs a steady state condition. When this occurs the program goes to the next realization. A large database with hundreds of realizations can take weeks of CPU time. The generated files are:

- *generated_inputs.dat*: this is a large table with $N_{rea}$ lines, where $N_{rea}$ denotes the number of realizations selected by the user in the *entrada.dat* file. In each line the file the program prints the random values of $\phi$, $H_0$, $\omega$ and $a$ that will be simulated;
    
- *outputs.dat*: after the convergence of each realization, the same input variables are printed here with two new information: the steady state temperature at the center of the tumor and the time the system takes to reach this condition;

## MLTA - Matlab version

This code simply reads the database in a txt file named DS.txt, learn the patterns with these data and apply the regressor algorithms to predict values that are not in the trainning data. It also check for the performance of the regressor schemes. The hyperparameters used in the machine learning algorithms are estimated through another code that uses Bayesian optimization.

## MLTA - Python version

This code is an extension of the MATLAB version of the MLTA programm, but focused on the regressor that provided the best performance, the ANN. It also comes with a graphic user interface developed using the customtkinter library where the user can define a single set of input parameters and after clicking on a button the algorithm prints the output parameters predicted by the ANN.

In order to use this version the user must install the following packages:

- matplotlib;
- customtkinter;
- sklearn;

In order to install matplotlib you could type the following command in a terminal (for Debian/Ubuntu users):

`sudo apt-get install python3-matplotlib`

To install customtkinter please type (for Debian/Ubuntu users):

`pip3 install customtkinter`  

To install the scikit-learn (sklearn) package in the Spyder software, use the following steps:

- First install the miniconda anaconda distribution in your machine;
- Open the anaconda prompt from miniconda and create an environment with the command `conda create -n spyder-env -y`
- Activate the environment with `conda activate spyder-env`
- Then install the scikit-learn package: `pip install scikit-learn`
- Then, open the Spyder software and click in "Tools" => "Preferences" => "Python interpreter";
- Click on "use the following Python interpreter" and choose *spyder-env/python.ex*
- Click on apply. Then, click in *restart kernel* in the options button on the right corner;

The package can now be used.



## References

[1] Gontijo, Rafael Gabler, and Andrey Barbosa Guimarães. "Langevin dynamic simulations of magnetic hyperthermia in rotating fields." Journal of Magnetism and Magnetic Materials 565 (2023): 170171. [DOI: 10.1016/j.jmmm.2022.170171](https://doi.org/10.1016/j.jmmm.2022.170171).

[2] Salloum, Maher, Ronghui Ma, and Liang Zhu. "An in-vivo experimental study of temperature elevations in animal tissue during magnetic nanoparticle hyperthermia." International Journal of Hyperthermia 24.7 (2008): 589-601. [DOI: 10.1080/02656730802203377](https://doi.org/10.1080/02656730802203377).