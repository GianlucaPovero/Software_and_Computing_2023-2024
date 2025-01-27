# Software and Computing project 2023/2024 #

## Pourpose ##
This code was made to fit the distribution of the D<sup>*</sup> into KK.
Trying to do a simultaneous fit of the two possible dostributions D<sup>+</sup> and D<sup>-</sup> and then computing the chi<sup>2</sup>/n<dub>DoF</sub> of the combined fit.

An example of histogram.root can be found in the repository and also an example of the results for both plots and logbook.

## Requirements and Usage ##
In order to run the code it is necessary to have the latest version of ROOT, which is esily found on line.
Onece dowloaded root, it is suggested to have the script.cpp and the histogram.root file in the same directory but most importatly to open root inside the same
directory of the script.cpp.

In order to run the code it is suggested to go in the directory where the script.cpp is and write "root -l script_name.cpp", this will make it so that root will open without 
printing the infromations about its current version and will outomatically lounch the script.

## Notes ##
The code is not perfect and does not work s intended, meaning that the fit computed by the code does not converge perfectly. The rason for this behavior cabn be found 
in the large statistics that this code is meant to analyze, reaching the order of the 10e6. 

However unstable, the code reaches an accetable result in MIGRAD = 1 and HESSE = 1, meaning that although the values of the paraemters are non perfect and the minimum is not reaced precisely, the errors computed by the hessian matrix are good. 

The example result in the repository has been found after two runs of the code.
