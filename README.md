optimizator
===========
- First run :

reduceOpttree.C 

in order to flatten the ntuple. Need to properly set the weights if merging different MC samples.

- Before doing optimization update Variables.hh file in order to select the variables to be used. For each variable the type of cuts (e.g. upper or lower bound type) has to be set with the boolean. Also the initial limits can be set.
- In optimizator_v3.C defines the category to use and the number of iterations to perform

Usage:

root -l

.L VarCat.cc++

.L optimizator_v3.C++

optimizator_v3(cat_number)

Once the optimization is completed the working points are stored in a ROOT file in the form of TCut and can be checked
by using check_WP.C macro.
