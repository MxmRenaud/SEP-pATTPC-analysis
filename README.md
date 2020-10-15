# SEP-pATTPC-analysis
Code used for the analysis of the 8Li + 40Ar experiment @NSL, ND, SEP/09-16/2019

firstMacro.C is a simple routine to be executed or compiled as a root macro. Root version6+. To use the deconvolution option requires compilation

The routine works by activating/deactivating a series of options in the form of boolean variables at the beggining of the main function. A series of warnings are in place in case conflicting options are selected.
A set of extra variables containing constants for the runtime is found right underneath the "options". The variables names' should help with reading.

It does require an extra file containing the list of runs to be read.

DISCLAIMER/WARNING: some of the error messages are needlessly rude/aggressive. Will be fixed once calmer waters have been reached.


