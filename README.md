# SEP-pATTPC-analysis
Code used for the analysis of the 8Li + 40Ar experiment @NSL, ND, SEP/09-16/2019

firstMacro.C is a simple routine to be executed or compiled as a root macro. Root version6+. To use the deconvolution option requires compilation.

The routine works by activating/deactivating a series of options in the form of boolean variables at the beggining of the main function. A series of warnings are in place in case conflicting options are selected.
A set of extra variables containing constants for the runtime is found right underneath the "options". The variables names' should help with reading.

It does require an extra file containing the list of runs to be read.
Inside the 'requiredFiles' folder one can find a series of files required for some of the options in firstMacro.C (unsurprisingly). 
The 'makeList.C' and 'mapping.C' were used pre-run to create the 'list-analysis.txt' and 'TPC-SEP2019-fullMap.dat' files.

The 'braggCalculations' folder contains a 'braggPeakCalculator.C' and all the files required to run it. It was used to predict the fusion signal's characteristics based on SRIM and PACE4 calculations. It is run as an invoked root macro, with a few variables to select at the beginning of the file, but the energy at which the fusion is happening has to be given while running (prompted).

