function run_wt(corticalClefts,nonlinearity,movement,werR,DEGL3,noise,CSs,Ss,distC,distS)

%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corticalClefts - vector of cell type's according to cortical celft position ["H", "N",..,"H"]
% nonlinearity   - where are nonlinearities                         [linear, myb, cpc, complex]
% movement       - how CPC moves by directed movement or diffusion                   [DM, Diff]
% werR           - protein or complex repressing wer                         [complex, protein]
% DEGL3          - does EGL3 diffuse or not                                              [Y, N]
% noise          - level of noise in initial conditions
% CSs            - cortical signal strength
% Ss             - scrambled strength
% distC          - cortical signal distribution                       [A (binary), B (diffuse)]
% distS          - SCM distribution                                   [A (binary), B (diffuse)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

load parspace.mat

%% coordinates

Time = 10000;                    % total simulation time
dt = 0.001;                      % timestep
T = 0:dt:Time;                   % time vector
NT = length(T);                  % number of timesteps

L = length(corticalClefts);      % number of cells
dx = 1;                          % spatial step
x = dx:dx:L;                     % space vector
Nx = length(x);                  % number of space steps

H = find(corticalClefts=="H");   % Hair cell positions
NH = find(corticalClefts=="N");  % Nonhair cell positions

%% run from start or continue if stopped

numParsets = 20000;
numVars = 10;

filename = sprintf('full_results_wt_%s_%s_%s_%s.mat',nonlinearity,movement,werR,DEGL3);
if isfile(filename)
    load(filename)
    ps1 = min(find(fullResults(:,1)==0));
else
    fullResults = zeros(numParsets,4); % [result category; pass threshold test; if reached steady state; if simulation broke]
    fullCellSpec = zeros(numParsets,numVars,L); % cellspec for all parsets, all variables and all cells
    sumCellSpec = zeros(numParsets,L); % cellspec sum for all parsets and all cells
    ps1 = 1;
end

%% asign cortical signal

if distC == "A"
    CS = zeros(Nx,1)';
    CS(:,H) = CSs;
elseif distC == "B"
    CS = CS_input(L,H);
end

if distS == "A"
    S = zeros(Nx,1)';
    S(:,H) = Ss;
elseif distS == "B"
    S = CS_input(L,H);
end

%% loop

for j = ps1:numParsets

    simulate_wt
    j
    
end

%% save
filename = sprintf('full_results_wt_%s_%s_%s_%s.mat',nonlinearity,movement,werR,DEGL3);
save(filename,'fullResults','fullCellSpec','sumCellSpec','-v7.3')

% find parameter sets that are yeses
wt_yesResult = find(fullResults(:,1)==1);

% find parameter sets that passed the threshold test, reached steady state and did not break
wt_passThresholdTest = find(fullResults(:,2)==1);
wt_reachedSteadyState = find(fullResults(:,3)==1);
wt_didNotBreak = find(fullResults(:,4)==0);

% find intersect to get successful parameter sets
wt_passes = intersect(intersect(wt_passThresholdTest,wt_reachedSteadyState),wt_didNotBreak);
wt_yes = intersect(wt_yesResult,wt_passes);

% save wild type results
filename = sprintf('yesParsets_wt_%s_%s_%s_%s.mat',nonlinearity,movement,werR,DEGL3);
save(filename,'wt_yes');

end