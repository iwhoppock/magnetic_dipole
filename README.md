# Magnetic Dipole/Magnetic Mirror/Three Adiabatic Invariants

This models a charged particle moving in a magnetic dipole, and shows all three adiabatic invariants (magnetic moment, longitudnal invariant/magnetic mirror, and the toroidal drift/magnetic flux). A visualisation of the output can be seen in dipole_bw_1.eps. 

All three algorithms make use of the Boris algorithm for non-relativistic particle pushing.

Three versions: (1) C; (2) Python; and (3) Matlab.

For it to work properly, the initial magnetic field at the starting location of the particle must be 1. Changes to the initial position will likely change the parameters of the magnetic field (all one must do is print out the initial field, and then multiply/divide the components of said field to make it one initially.) 

The units are normalised by the: speed is in terms of the Alfven speed, frequencies the proton cyclotron frequency, lengths the proton inertial length, etc. To change to an alpha, replace mass with 4 and charge with 2. Change the velocities gradually. This code implements the Boris algorithm for non-relativistic particle tracing. The quantity vAc is the ratio of the alfven speed by the speed of light (c). 

Obviously the C code will be much faster. 

## The C code

I make use of autotools. Compling instructions are: (1) autoreconf -1 (2) ./configure (3) make (4) ./dipole

The output is a csv file, and I have included a matlab plotting script (easily adopted for matplotlib).

![alt text](https://github.com/iwhoppock/magnetic_dipole/blob/master/dipole_bw_1.eps)

