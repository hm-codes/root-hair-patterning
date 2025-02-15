function run_M(corticalClefts,nonlinearity,movement,werR,DEGL3,noise,signal,distC,distS)

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

%% find parameter sets to test

% load wild type yeses
filename = sprintf('WildtypeResults/yesParsets_wt_%s_%s_%s_%s.mat',nonlinearity,movement,werR,DEGL3);
load(filename)

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

numParsets = length(wt_yes);

filename = sprintf('full_results_M_%s_%s_%s_%s_noise%.2f.mat',nonlinearity,movement,werR,DEGL3,noise);
if isfile(filename)
    load(filename)
    tracker = max(find(fullResults(:,1)~=0))+1;
    ps1 = wt_yes(tracker);
    load(sprintf('tables_M_%s_%s_%s_%s_noise%.2f.mat',nonlinearity,movement,werR,DEGL3,noise))
else
    fullResults = zeros(numParsets,length(signal));
    tables = zeros(numParsets,length(signal),3,4);
    tracker = 1;
    ps1 = 1;
end

%% loop

for j = ps1:20000

    if ismember(j,wt_yes)

        for jj = 1:length(signal)

            if movement == "Diff"
                % scrambled strength varies
                Ss = signal(jj);
                CSs = 0;
            elseif movement == "DM"
                % cortical signal strength varies
                CSs = signal(jj);
                Ss = 0;
            end

            % asign cortical signal
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
    
            % simulate
            simulate_M

                % test for mutant pattern
                scm_results(1,:) = mean(results_GL2);
                scm_results(2,:) = mean(results_cpc);
                scm_results(3,:) = mean(results_WER);
                %scm_results
                tables(tracker,jj,:,:) = scm_results;

                % if homogeneous, stop
                if any(any(scm_results>100))
                    fullResults(tracker,jj) = 999;           %"homogeneous"
                end
                if any(any(scm_results>100))
                    break
                end
    
                % if pass then finish
                if sum(sum(ttest_count_results)) == 12 || sum(sum(rangetest_count_results)) == 12
                    fullResults(tracker,jj) = 1;                    %"pass"
                end
                if sum(sum(ttest_count_results)) == 12 || sum(sum(rangetest_count_results)) == 12
                    break
                end

                % not failed but not passed so lower cortical signal
                fullResults(tracker,jj) = 2;                    %"continue"
           
        end

        % save results every 100 parsets
        if mod(tracker,100) == 0
            filename = sprintf('full_results_M_%s_%s_%s_%s_noise%.2f.mat',nonlinearity,movement,werR,DEGL3,noise);
            save(filename,'fullResults','-v7.3')
            filename = sprintf('tables_M_%s_%s_%s_%s_noise%.2f.mat',nonlinearity,movement,werR,DEGL3,noise);
            save(filename,'tables','-v7.3')
        end

        tracker = tracker+1;

    end
    
    j
    
end

%% save
filename = sprintf('full_results_M_%s_%s_%s_%s_noise%.2f.mat',nonlinearity,movement,werR,DEGL3,noise);
save(filename,'fullResults','-v7.3')
filename = sprintf('tables_M_%s_%s_%s_%s_noise%.2f.mat',nonlinearity,movement,werR,DEGL3,noise);
            save(filename,'tables','-v7.3')

% find parameter sets that are yeses
M_yes = [];
for i = 1:length(signal)
    M_yes = [M_yes; find(fullResults(:,i)==1)];
end
    
% save mutant results
filename = sprintf('yesParsets_M_%s_%s_%s_%s_noise%.2f.mat',nonlinearity,movement,werR,DEGL3,noise);
save(filename,'M_yes');

end