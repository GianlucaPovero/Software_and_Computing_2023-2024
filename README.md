# Software and Computing project 2023/2024 #

## Pourpose ##
This code was made be me in order to fit the distribution of the D<sup>*</sup> into KK.
Trying to do a simultaneous fit of the two possible distributions D<sup>+</sup> and D<sup>-</sup> and then computing the chi<sup>2</sup>/n<sub>DoF</sub> of the combined fit.

The code uses the RooFit package to define the functon and plot the distributions, but uses RooMinuit when it comes to computing the actual Chi<sup>2</sup> by defining a 
RooChi2Var. This gives the fit a better control over the fitting action.

In order to make the code user freandly an example of a histogram.root can be found in the repository. Also, an example of output result for both plots and logbook has been uploaded
and can be found in the repository.

## Requirements and Usage ##
In order to run the code it is only necessary to have the latest version of ROOT, which is esily found on line.
Onece dowloaded root, it is suggested to have the script.cpp and the histogram.root file in the same directory but most importatly to open root inside the same
directory of the script.cpp.

In order to run the code it is suggested to go in the directory where the script.cpp is and write `root -l script_name.cpp`, this will make it so that root will open without 
printing the infromations about its current version and will automatically launch the script.

## Notes ##
The code is not perfect and does not work as intended, meaning that the fit computed by the code does not converge perfectly. The rason for this behavior can be found 
in the large statistics that this code is meant to analyze, reaching the order of the 10e6 events. 

In the code you will find that the model pdf for D<sup>-</sup> mostly uses parameters defined for D<sup>+</sup>, the reason for this is that the fit is non optimized to perfection and
RooFit itself has a hard time fitting with a high number of free parameters. In order to reduce this effect, most of the variables are "shared" between the two models.

However, all parameters for the D<sup>-</sup> model are defined beacuse the final version of the code is yet to be found, and surely the best configuration will be a mix 
of both shared and independent variables between the two models.

Although unstable, the code reaches an accetable result as MIGRAD = 1 and HESSE = 1, meaning that although the values of the parameters are non perfect and the minimum is not reached precisely, the errors computed by the hessian matrix are good. 

Latsly, it is woth mentioning that the example result in the repository has been found after two runs of the code with only the means of the Johnson functions not shared. 
