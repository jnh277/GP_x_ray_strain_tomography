# A Bayesian Approach to Triaxial Strain Tomography from High-energy X-ray Diffraction

A Gaussian process based approach for the reconstruction of strain from high-energy X-ray strain tomography. This repository contains matlab code to run the simulation example given in the paper available at https://arxiv.org/abs/1903.02158. The code simulates the generation of high-energy X-ray strain images where the sample is rotated about only a single axis. 

The matlab script run_example.m will simulate the measurements and then run the reconstruction. The number of sample rotation angles can be changed by a parameter at the top of the script.

Note: the added measurement noise is random and so the results will vary slightly with each run, particularly if the number of rotations is set low (i.e. a low number of measurements is used).
