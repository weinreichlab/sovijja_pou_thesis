README (Sovijja Pou Thesis):

The Excel file “initialWorkPyrimethamineProbabilitiesandTrajectory” holds the initial work I did on the topic including preliminary modeling of the PK equation. The Excel file “OgbunugaforAlleleGR” is the growth rate data I used from Ogbunugafor et al. 2016.

All the functions and code used to generate the results found in the thesis are in the MATLAB and Mathematica Files. 

For the Matlab results, the two major scripts used are “finalGenerateResults” and “probTrajectory.” “finalGenerateResults” is how I generated the fitness landscapes with various functions, saved them, and looked at average/standard deviation of fitness results as well as fitness peak statistics. “probTrajectory” is our incomplete/imperfect method of generating probabilities of trajectories

Matlab:

arithavgfit: a function used to calculate the fitness of an allele under the IV drip/no variance regime.

fit_func: function used to calculate the fitness of an allele under the standard PK regime.

fit_landscape: Generates various fitness landscapes

p_fix: Calculates probability of fixation

peak_analysis: gives peak statistics of fitness landscapes

sel_coeff, selmatrix: calculates selection coefficients

8000RunsGeoFit: stores the landscapes used to generate the Results section of the thesis

PeakInfo: stores the peak statistics of the landscapes used in the Results section.



Mathematica:

Holds the functions and plots used to study Jensen’s Inequality in the thesis.


