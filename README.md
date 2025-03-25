# Fitting scaling functions to the data

In polymer physics scaling law is observed for a number of properties. In particular size characteristics like: gyration radius, hydrodynamic radius, end-to-end distance, spanning radius. There are also a number of other characteristics that can be fitted

## Functions awaliable 
- Scaling exponent and amplitude: the function to calculate scaling exponent and amplitude for a function like a gyration radius for a relatively small data set.  This function however does not account for corrections to scaling. It takes a file with two columns in lin-lin scale or log-log scale.
- Scaling form RG chain: In general gyration radius can be presented as $R^2_g=A(1+bN^{-\Delta})N^{2\nu}$. Values of the exponents are taken from the paper by  N.Clisby, PRL 104, 055702 (2010). 
This function takes an input file with two columns one containing number of monomers (molecular mass) and second ether $R_g^2$ or $R_g$. With order of gyration radius specified by the second parameter.
- Rg chain prediction: This function can predict a value of the gyrtion radius for linear polymer with a  given number of monomers in good solution ($\nu=0.5882$). There are build in approximation functions for Langeven dynamic and lattice Monte-Carlo on cubic lattice. It can also give a prediction based on a new set of data. The function takes a number of monomers as parameter. And in case of new set of data it is also required.
- Dynamic functions: Dynamic function $g_1$ and $g_3$ have a dependence $\sim t$ on large time scales. This function is designed to fined a region where this dependence is true and to calculate an amplitude in this region. This function takes in a file with two columns one is simulation time and second eather $g_1$ or $g_3$ in linear scale. It also takes step size (number of points that are added to lower cut off on each iteration step), a minimum size of the fiting regeon and an accuracy that is how close the exponent has to be to 1. 
- Is there scaling: this function takes a file with two column (A,B) and tests wether there is a scaling dependance between columns, or $B \sim \exp(\alpha A)$
- Finit size scaling: This function is desighned to calculate finit size scaling for universal size ratios. Takes in a file with two columns and returns a value of the size ration for an infinitely long chain

Some functions have an option to plot the fit.

## Files
- Functions.py contains functions described above plus few support functions
- Scaling_calc.py is a main file
- re_MC.txt is a sample file that contains number of monomers in the first column and the squared end to end distance in the second. This data was received for self avoiding walks on cubic lattice using pivot algorithm for chains of up to 1000 steps.
- re_MC.txt is a sample file that contains number of monomers in the first column and the squared gyration radius in the second. This data was received for self avoiding walks on cubic lattice using pivot algorithm for chains of up to 1000 steps.
- rg_MD_chain.txt is a sample file that contains number of monomers in the first column and the squared gyration radius in the second. This data was received for linear chains with up to 900 beads in continuous space using Langeven dynamics with a repulsive potential. Simulations were conducted using LAMMPS.
- re_vs_rg.txt is a sample file that contains number of monomers in the first column and $R_e^2/R_g^2$ in the second. This data was received for self avoiding walks on cubic lattice using pivot algorithm for chains of up to 1000 steps.

## Requirements

To run the main python file following libraries are required: 
- tkinter
- matplotlib
- numpy
- scikit-learn
- pandas
- statsmodels

## Note!
When using the program for scientific research please site the papers if the citation is provided as an output of the function
