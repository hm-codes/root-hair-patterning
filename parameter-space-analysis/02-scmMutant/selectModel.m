%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ALL
% cell type according to cortical celft position
corticalClefts = ["H","N","H","N","H","N","H","N","N","H","N","N","H","H","N","N","H","N","N"];

% where are nonlinearities [linear, myb, cpc, complex, cpccomplex, mybcomplex, mybcpc, all]
nonlinearity = "all";

% movement of CPC [DM, Diff]
movement = "DM";

% competitive inhibition or protein/complex repressing wer [CI, complex, protein]
werR = "complex";

% does EGL3 diffuse or not [Y, N]
DEGL3 = "N";

%% INPUT DISTRIBUTIONS
% cortical signal distribution [A (binary), B (diffuse)]
distC = "A";

% SCM distribution [A (binary), B (diffuse)]
distS = "A";

%% MUTANTS

% noisy initial conditions
noise = 1.0;
signal = [0.25, 0.125, 0.0625, 0.0312, 0.0156, 0.0078, 0.0039, 0.0019, 0.0010, 0.0005];

% run
run_M(corticalClefts,nonlinearity,movement,werR,DEGL3,noise,signal,distC,distS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%