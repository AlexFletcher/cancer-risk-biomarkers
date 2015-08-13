# cancer-risk-biomarkers
Code associated with the manuscript "A computational modelling approach for deriving biomarkers to predict cancer risk in premalignant disease"

The provided code is in a number of files, which should all be stored in the same directory.

A description of the files is provided below:

runner.m — highest level code file, from which multiple simulations with different parameters can be run. For each set of parameters, calls on the summary_stats function to aggregate data over N_runs simulations of the model.

summary_stats.m — second-highest level code file, which, when given a set of parameters, runs the Gillespie algorithm code for the Moran model N_runs times, and aggregates the summary statistic and endpoint data of each run into cell arrays. In this code, the parallelization of the code is done, in that the N_runs runs of the simulation are run in parallel when possible.

To use this code, run the following: 

summary_stats(fName, s_pos, s_del, mu, N_muts_necessary, bx_size, weights_matrix, total_weight_sum, weights_matrix_bx, total_weight_sum_bx)

Where fName - filename of the file for all of the saved data, s_pos, s_del are the fitness advantage and disadvantage respectively of positive and negative mutations, mu is the mutation rate, N_muts_necessary - cutoff number of positive mutations for cancer, bx_size - radius of biopsy, weights_matrix - weight matrix, total_weight_sum - weight matrix sum, weights_matrix_bx - biopsy weight matrix, total_weight_sum_bx - biopsy weight matrix sum. Note that the last 4 variables must be loaded into memory before use of this function by loading the appropriate .mat files.

new_model2d.m - the primary code file, containing the Gillespie algorithm details for the spatial simulation, returning a number of the variables such as the endpoint or biopsy statistic measures. Lines 53-58 in the code represent an important variable that can be changed to decide whether multiple biopsies should be included at one time point or another. Otherwise, biopsies occur just once per time point (default setting).

diversityMeasures.m - a helper function returning the Shannon and Simpson indices, given the list of cell genotypes and numbers of each type of genotype

calcMoranI.m - a helper function computing Moran’s I, requiring the input from the appropriate weight matrix variable

gearyC.m - a helper function computing Geary’s C, requiring the input from the appropriate weight matrix variable, and cellular lattice, analogous to that for Moran’s I

fitness_prolif_index.m - a helper function computing the proliferative indices.

The .mat files called weight_matrix.mat, weight_matrix_bx.mat, weight_matrix_bx_small.mat, and weight_matrix_bx_large.mat are all data files containing the necessary weight matrices, and total sum of weights as defined in the manuscript. These are necessary for the computations of Moran’s I and Geary’s C. Note that the weight_matrix file contains the weight matrix for the full lattice of size 100x100, weight_matrix_bx contains the weight matrix for the biopsy of radius 20, weight_matrix_bx_small contains the weight matrix for the biopsy of radius 5, and weight_matrix_bx_large contains the weight matrix for the biopsy of radius 40.

The make_weight_matrix_bx.m function creates the weight matrix for a circular biopsy with radius specified by the variable N_bx_radius. This must be run before the main program, so that runner.m can load the file containing the weight matrix (or so that summary_stats.m can use the weight matrix, that is loaded by the user).

The make_weight_matrix.m function creates the weight matrix for a the lattice for calculation of Moran’s I and Geary’s C, with lattice size set by N_size. This must be run before the main program, so that runner.m can load the file containing the weight matrix (or so that summary_stats.m can use the weight matrix, that is loaded by the user).

The .mat file called initial_locs.mat is a matrix containing all of the initial cell locations, which is initially every grid point, to speed up computations by having the variable pre-computed.

GRAPHING FUNCTIONS

km_plotting.m - The main graph/table generating function. To use this, set the variables in the first 2 lines manually, (i.e. file to be loaded containing the simulation output, and also the number of positive mutations necessary for cancer, if plotting. Lines 72-87 MUST be uncommented if the data generated includes data for multiple biopsies at the same time point and plotting anything other than the data for that time point. That is, this code should be uncommented only when using the data set with extra biopsies at a single time point, and generating a Hazards Ratio table, KM Curves, or correlation coefficient-time curves. 

For the generation of the hazards ratio tables, the times at which the data is sampled must be set in lines 101-102.

To make Kaplan-Meier curves, lines 163 - 237 should be uncommented.

To make Correlation coefficient vs Time graphs, for lattice, biopsy, and scrapings data, lines 241-444 should be uncommented.

To make 2D heatmaps of the difference index, using biopsy information from a prior time point, lines 447-523 should be uncommented.

To make the correlation coefficient vs number of biopsies plot, lines 527-583 should be uncommented.

makeKMgraph.m - a helper function for creating the Kaplan-Meier curves, as well as the input data for the R script performing the log rank test. Note that in order for R to use this csv, the file must be manually modified to have the header "end_time" for column 1 and "quartile" for column 2.

logrank_r_script.txt - contains the script that can be used with R Commander, when set to the appropriate directory containing the CSV file, when referenced to the correct file. Note that the CSV must have the appropriate headers of “end_time” for column 1 and “quartile” for column 2, that must be manually added, for this script to work.
