FFT_MA
======

Gaussian random fields by fast fourier transform using the moving average approach

Based on the journal paper: The FFT Moving Average (FFT-MA) Generator: An Efficient Numerical Method for Generating and Conditioning Gaussian Simulations,
http://dx.doi.org/10.1023/A%3A1007542406333 by Mickaële Le Ravalec, Benoît Noetinger and Lin Y. Hu.

This implementation supports full 3D Gaussian random fields and cosimulation. See the main files

The prerequests for the this project is fftw3 which is avaliable in every linux distribution. For the C++ implementation there is also a second prerequest called fftw++ - http://fftwpp.sourceforge.net/ which should be downloaded into a superfolder compared to this folder. If it is installed in another place - modify the Makefile.

Windows and Mac users you are on your own.

The fortran version uses the fairly new feature of classes in fortran 2003 - the gfortran compiler supports this from 4.6 mayby earlier versions will also compile and run, however, I have not tried it out.
