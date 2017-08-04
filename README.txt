Codes are organised as follows:



ABC-rej scripts: 
	
	We simulate our 2D lattice based model using ABC_mex.cpp (compiled using mex ABC_mex.cpp, binary works for ubuntu)
	
We perform ABC and save out using 'ABC_COMBINED.m'
	
We loop through the 10 replicates, excluding each one to compute distances and then the corresponding information gain for posterior due to those distances. 
This is done in 'ComputeJackknifing.m'. We save out the distances and information gains into the "ABC-rej data" folder

 


Plotting scripts:
	
	'Plots_For_Manuscript.m' is the file containing the different plotting scripts for plotting single posteriors/bar plots/prettying the plots for the manuscript.



ABC-DC:
	
	Contains the ABC-DC algorithm, the data from the algorithm for statistic CXY, and a means of plotting the results shown in the report.





ABC-rej data:
	Note: This folder is available from http://datadryad.org/submit?journalID=RSOS&manu=RSOS-171045
	Contains 'infoRE_SD#DF#.mat' resulting from 'ComputeJackknifing.m' for varying Scratch Densities and Density Fractions

