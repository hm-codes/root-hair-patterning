'selectModel.m' runs the simulation for any chosen model. The user input parameters are as follows:
parameter set, j - any of the 20,000 parameter sets tested, successful parameter set numbers are given in, Figure 2B, Figure 3A and 3B, Table 9, and Supplementary file S4.
root - "wildtype" or "mutant" for the scm mutant
nonlinearity - can be none; "linear", a single nonlinearity; "myb"; "complex"; or "cpc", a double nonlinearity; "mybcpc"' "mybcomplex"; or "cpccomplex", or all three; "all".
movement - (of CPC), eiter by both diffusion and directed movement, "DM", or diffusion only "Diff".
werR - the mechanism of wer repression which can be via the CPC complex, "complex", CPC protein, "protein", or competitive inhibition, "CI".
DEGL3 - can EGL3 diffuse "Y" or "N".
signal - this can be found for any of the successful parameter sets in column 6 of 'SuccessData.mat'
noise - this can be found for any of the successful parameter sets in column 7 of 'SuccessData.mat'

The script calls other MATLAB files:
'model_parameters.m' - assigns the parameter values according to both the parameter set and specified model
'ftcs.m' - a forward time cenetered space finite difference scheme for the diffusion of any diffusive components
'rangetest_count.m' and 'rangetest_count_WT.m' - test the simulaton results against experimental data. The output of these files is in the form:
HH HN NH NN
giving a binary 1x4 matrix, where 1 indicates the particular cell type matches experimental data and 0 otherwise.
The command line output of the programme is two 3x4 matrices. The first one is binary as described above with the results of the test against experimental data for GL2 (row 1), cpc (row 2), and WER (row 3). The second matrix is the average of the 20 simulation repeats, giving numbers comparable to those in Table 3, with the rows matching the order of components described above.
'crossSection_plot.m' - plots a cross section of a simulated root (set to the last of 20 random initial conditions), with the cortical cell positions shown as input into the model. It also plots a timeseries of GL2 concentration in each cell, by proxy of the total activator complex.
