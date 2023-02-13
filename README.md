# multigaussfit
Python code to perform line fitting on a single spectrum
based on lmfit  https://lmfit.github.io/lmfit-py/ 
you can find many other models to use on the lmfit page (this code uses only gaussian curves)

To use: >python gaussfit.py file.dat
file.dat is an ascii file with 3 columns: WL, flux, flux_error

It is important to define the model with the amount of single gaussian components
and to update the 'prefixes' for each single gaussian, also in the parameter definition stage

the printing is done in LaTeX format
