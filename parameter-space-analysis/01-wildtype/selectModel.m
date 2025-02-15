%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% noisy initial conditions
noise = 0;

% cortical signal strength
CSs = 1;

% scrambled strength
Ss = 1;

% cortical signal distribution [A (binary), B (diffuse), Z (none)]
distC = "A";

% SCM distribution [A (binary), B (diffuse), Z (none)]
distS = "A";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_wt(corticalClefts,nonlinearity,movement,werR,DEGL3,noise,CSs,Ss,distC,distS)