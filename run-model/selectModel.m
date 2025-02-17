%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter set number
j = 15122;

% wildtype or mutant
root = "mutant";

%% SPECIFY MODEL

% where are nonlinearities [linear, complex, cpc, myb, mybcpc, cpccomplex, mybcomplex, all]
nonlinearity = "all";

% movement of CPC [DM, Diff]
movement = "DM";

% wer downregulation [complex, protein, CI]
werR = "complex";

% does EGL3 diffuse or not [Y, N]
DEGL3 = "N";

%% MUTANTS

% signal for chosen parameter set
signal = 0.25;

% noise for chosen parameter set
noise = 0.5;

%% INPUT DISTRIBUTIONS

% cell type according to cortical celft position
% experiment
corticalClefts = ["H","N","H","N","H","N","H","N","N","H","N","N","H","H","N","N","H","N","N"];
% wt1
% corticalClefts = ["H","N","N","H","N","H","N","H","N","N","H","N","H","N","N","H","N","H","N"];
%scm1
% corticalClefts = ["H","N","H","N","H","N","H","N","N","H","N","N","H","N","H","N","N"];

% cortical signal distribution [A (binary), B (diffuse)]
distC = "A";

% SCM distribution [A (binary), B (diffuse)]
distS = "A";

%% run

%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corticalClefts - vector of cell type's according to cortical celft position ["H", "N",..,"H"]
% nonlinearity   - where are nonlinearities                         [linear, complex, cpc, myb,
%                                                          mybcpc, cpccomplex, mybcomplex, all]
% movement       - how CPC moves by directed movement or diffusion                   [DM, Diff]
% werR           - downregulation of wer                                 [complex, protein, CI]
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

%%  assign signal strength

if movement == "DM"
    CSs = signal;
    % scrambled strength
    if root == "wildtype"
        Ss = 1;
    elseif root == "mutant"
        Ss = 0;
    end
elseif movement == "Diff"
    CSs = 0;
    % scrambled strength
    if root == "wildtype"
        Ss = 1;
    elseif root == "mutant"
        Ss = signal;
    end
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

%% simulate

% generate parameters
model_parameters

% number of experimental repeats
numRepeats = 20;

% set up storage matrices
results_GL2 = zeros(numRepeats,4);
results_cpc = zeros(numRepeats,4);
results_WER = zeros(numRepeats,4);

concentrations = zeros(numRepeats,16,L);
timecourse = zeros(numRepeats,L,NT-1);

%% experiment

for seed = 1:numRepeats

    % initial conditions

    rng(seed)
    % proteins
    G = rand(Nx,1)*noise;    % GL3
    E = rand(Nx,1)*noise;    % EGL3
    C = rand(Nx,1)*noise;    % CPC
    W = rand(Nx,1)*noise;    % WER
    M = rand(Nx,1)*noise;    % MYB23
    % mRNA
    g = rand(Nx,1)*noise;    % gl3 mRNA
    e = rand(Nx,1)*noise;    % egl3 mRNA
    c = rand(Nx,1)*noise;    % cpc mRNA
    w = rand(Nx,1)*noise;    % wer mRNA
    m = rand(Nx,1)*noise;    % myb23 mRNA
    % complexes
    A_GW = rand(Nx,1)*noise; % activator complex GL3-WER
    A_EW = rand(Nx,1)*noise; % activator complex EGL3-WER
    A_GM = rand(Nx,1)*noise; % activator complex GL3-MYB23
    A_EM = rand(Nx,1)*noise; % activator complex EGL3-MYB23
    I_GC = rand(Nx,1)*noise; % activator complex GL3-CPC
    I_EC = rand(Nx,1)*noise; % activator complex EGL3-CPC


    % directed movement function
    fDM = @(S,C) k3*circshift(C,1).*S' + k3*circshift(C,-1).*S'...
        - k3*C.*circshift(S,1)' - k3*C.*circshift(S,-1)';

    % ODE functions
    cf1 = @(G,W,A_GW) k11*(G.^n1).*(W.^n2) - k12*A_GW;          % formation of the activator complex
    cf2 = @(E,W,A_EW) k13*(E.^n1).*(W.^n2) - k14*A_EW;          % formation of the activator complex
    cf3 = @(G,M,A_GM) k15*(G.^n1).*(M.^n2) - k16*A_GM;          % formation of the activator complex
    cf4 = @(E,M,A_EM) k17*(E.^n1).*(M.^n2) - k18*A_EM;          % formation of the activator complex
    cf5 = @(G,C,I_GC) k21*G.*C - k22*I_GC;                      % formation of the inhibitor complex
    cf6 = @(E,C,I_EC) k23*E.*C - k24*I_EC;                      % formation of the inhibitor complex

    fG = @(g,G,W,A_GW,M,A_GM,C,I_GC) qG*g - n1*cf1(G,W,A_GW)...
        - n1*cf3(G,M,A_GM) - cf5(G,C,I_GC) - dG*G;              % reactions of GL3 protein
    fE = @(e,E,W,A_EW,M,A_EM,C,I_EC) qE*e - n1*cf2(E,W,A_EW)...
        - n1*cf4(E,M,A_EM) - cf6(E,C,I_EC) - dE*E;              % reactions of EGL3 protein
    fC = @(c,G,C,I_GC,E,I_EC) qC*c - cf5(G,C,I_GC)...
        - cf6(E,C,I_EC) - dC*C + fDM(S,C);                      % reactions of CPC protein
    fW = @(w,G,W,A_GW,E,A_EW) qW*w - n2*cf1(G,W,A_GW)...
        - n2*cf2(E,W,A_EW) - dW*W;                              % reactions of WER protein
    fM = @(m,G,M,A_GM,E,A_EM) qM*m - n2*cf3(G,M,A_GM)...
        - n2*cf4(E,M,A_EM) - dM*M;                              % reactions of MYB23 protein

    hg = @(A_T,g) bg - rg*A_T.*g - dg*g;                        % production and degradation of gl3 mRNA
    he = @(A_T,e) be - re*A_T.*e - de*e;                        % production and degradation of egl3 mRNA
    hc = @(A_T,c) pc*(A_T.^n4)./(kdC+A_T.^n4) - dc*c;                          % production and degradation of cpc mRNA
    hw = @(I_T,C,w) bw - rw11*C.*w - rw12*I_T.*w - rw21*S'.*w...
        - rw22*CS'.*w - dw*w;                                   % production and degradation of wer mRNA
    hm = @(A_T,m) pm*(A_T.^n3)./(kdM+A_T.^n3) - dm*m;                          % production and degradation of myb23 mRNA


    % loop over time
    broken = 0;
    steadyState = 0;
    for i = 2:NT

        % total inhibior and activator
        I_T = I_GC + I_EC;
        A_T = A_GW + A_EW + A_GM + A_EM;
        timecourse(seed,:,i-1) = A_T;

        % reactions
        G_new = G + dt*fG(g,G,W,A_GW,M,A_GM,C,I_GC);
        E_new = E + dt*fE(e,E,W,A_EW,M,A_EM,C,I_EC);
        C_new = C + dt*fC(c,G,C,I_GC,E,I_EC);
        W_new = W + dt*fW(w,G,W,A_GW,E,A_EW);
        M_new = M + dt*fM(m,G,M,A_GM,E,A_EM);

        g_new = g + dt*hg(A_T,g);
        e_new = e + dt*he(A_T,e);
        c_new = c + dt*hc(A_T,c);
        w_new = w + dt*hw(I_T,C,w);
        m_new = m + dt*hm(A_T,m);

        A_GW_new = A_GW + dt*cf1(G,W,A_GW);
        A_EW_new = A_EW + dt*cf2(E,W,A_EW);
        A_GM_new = A_GM + dt*cf3(G,M,A_GM);
        A_EM_new = A_EM + dt*cf4(E,M,A_EM);
        I_GC_new = I_GC + dt*cf5(G,C,I_GC);
        I_EC_new = I_EC + dt*cf6(E,C,I_EC);

        % diffusion
        G_new = ftcs(dt,dx,Nx,DG,G_new);
        E_new = ftcs(dt,dx,Nx,DE,E_new);
        C_new = ftcs(dt,dx,Nx,DC,C_new);

        % steady state check
        if all(abs([G_new,E_new,C_new,W_new,M_new,g_new,e_new,c_new,w_new,m_new,A_GW_new,A_EW_new,A_GM_new,A_EM_new,I_GC_new,I_EC_new]-[G,E,C,W,M,g,e,c,w,m,A_GW,A_EW,A_GM,A_EM,I_GC,I_EC])<=1e-8)
            steadyState=1;
        end
        if all(abs([G_new,E_new,C_new,W_new,M_new,g_new,e_new,c_new,w_new,m_new,A_GW_new,A_EW_new,A_GM_new,A_EM_new,I_GC_new,I_EC_new]-[G,E,C,W,M,g,e,c,w,m,A_GW,A_EW,A_GM,A_EM,I_GC,I_EC])<=1e-8)
            break
        end

        % broken simulation check
        % record broken sims
        if any(any(imag([G_new,E_new,C_new,W_new,M_new,g_new,e_new,c_new,w_new,m_new,A_GW_new,A_EW_new,A_GM_new,A_EM_new,I_GC_new,I_EC_new])~=0))...
                || any(any([G_new,E_new,C_new,W_new,M_new,g_new,e_new,c_new,w_new,m_new,A_GW_new,A_EW_new,A_GM_new,A_EM_new,I_GC_new,I_EC_new]<0))
            broken = 1;
        end
        % break loop
        if any(any(imag([G_new,E_new,C_new,W_new,M_new,g_new,e_new,c_new,w_new,m_new,A_GW_new,A_EW_new,A_GM_new,A_EM_new,I_GC_new,I_EC_new])~=0))...
                || any(any([G_new,E_new,C_new,W_new,M_new,g_new,e_new,c_new,w_new,m_new,A_GW_new,A_EW_new,A_GM_new,A_EM_new,I_GC_new,I_EC_new]<0))
            break
        end

        % update variables
        G = G_new;
        E = E_new;
        C = C_new;
        W = W_new;
        M = M_new;

        g = g_new;
        e = e_new;
        c = c_new;
        w = w_new;
        m = m_new;

        A_GW = A_GW_new;
        A_EW = A_EW_new;
        A_GM = A_GM_new;
        A_EM = A_EM_new;
        I_GC = I_GC_new;
        I_EC = I_EC_new;

    end

    concentrations(seed,:,:) = [G,E,C,W,M,g,e,c,w,m,A_GW,A_EW,A_GM,A_EM,I_GC,I_EC]';

    %% TESTS
    % make tables

    % total concentrations of markers
    GL2_T = A_GW + A_EW + A_GM + A_EM;    % total activator complex -> GL2
    cpc_T = c;                            % total cpc mRNA
    WER_T = A_GW + A_EW + W;              % total WER protein

    th = 1.1;                             % threshold

    % test GL2
    if max(GL2_T)/min(GL2_T) <= th
        results_GL2(seed,1:4) = [999,999,999,999];
    elseif broken == 1 || steadyState == 0
        results_GL2(seed,1:4) = [999,999,999,999];
    else
        % average concentration across all cells
        Ave_GL2 = mean(GL2_T);

        % decide if each cell has GL2 (1) or not (0)
        GL2 = zeros(1,L);
        for cell = 1:L
            if GL2_T(cell) >= Ave_GL2
                GL2(cell) = 1;
            else
                GL2(cell) = 0;
            end
        end

        GH = find(GL2==0);            % find indices of hair cells according to GL2 marker
        GN = find(GL2==1);            % find indices of non-hair cells according to GL2 marker

        results_GL2(seed,1) = length(intersect(H,GH))/length(H)*100;
        results_GL2(seed,2) = length(intersect(H,GN))/length(H)*100;
        results_GL2(seed,3) = length(intersect(NH,GH))/length(NH)*100;
        results_GL2(seed,4) = length(intersect(NH,GN))/length(NH)*100;
    end

    % test CPC
    if max(cpc_T)/min(cpc_T) <= th
        results_cpc(seed,1:4) = [999,999,999,999];
    elseif broken == 1 || steadyState == 0
        results_cpc(seed,1:4) = [999,999,999,999];
    else
        % average concentration across all cells
        Ave_cpc = mean(cpc_T);

        % decide if each cell has cpc (1) or not (0)
        cpc = zeros(1,L);
        for cell = 1:L
            if cpc_T(cell) >= Ave_cpc
                cpc(cell) = 1;
            else
                cpc(cell) = 0;
            end
        end

        cH = find(cpc==0);            % find indices of hair cells according to cpc marker
        cN = find(cpc==1);            % find indices of non-hair cells according to cpc marker

        results_cpc(seed,1) = length(intersect(H,cH))/length(H)*100;
        results_cpc(seed,2) = length(intersect(H,cN))/length(H)*100;
        results_cpc(seed,3) = length(intersect(NH,cH))/length(NH)*100;
        results_cpc(seed,4) = length(intersect(NH,cN))/length(NH)*100;
    end

    % test WER
    if max(WER_T)/min(WER_T) <= th
        results_WER(seed,1:4) = [999,999,999,999];
    elseif broken == 1 || steadyState == 0
        results_WER(seed,1:4) = [999,999,999,999];
    else
        % average concentration across all cells
        Ave_WER = mean(WER_T);

        % decide if each cell has WER (1) or not (0)
        WER = zeros(1,L);
        for cell = 1:L
            if WER_T(cell) >= Ave_WER
                WER(cell) = 1;
            else
                WER(cell) = 0;
            end
        end

        WH = find(WER==0);            % find indices of hair cells according to WER marker
        WN = find(WER==1);            % find indices of non-hair cells according to WER marker

        results_WER(seed,1) = length(intersect(H,WH))/length(H)*100;
        results_WER(seed,2) = length(intersect(H,WN))/length(H)*100;
        results_WER(seed,3) = length(intersect(NH,WH))/length(NH)*100;
        results_WER(seed,4) = length(intersect(NH,WN))/length(NH)*100;
    end

end

%% tests

if root == "wildtype"

    rangetest_count_results(1,:)=rangetest_count_WT(mean(results_GL2));
    rangetest_count_results(2,:)=rangetest_count_WT(mean(results_cpc));
    rangetest_count_results(3,:)=rangetest_count_WT(mean(results_WER))

    % test for wildtype pattern
    wt_results(1,:) = mean(results_GL2);
    wt_results(2,:) = mean(results_cpc);
    wt_results(3,:) = mean(results_WER)

elseif root == "mutant"

    rangetest_count_results(1,:)=rangetest_count(mean(results_GL2));
    rangetest_count_results(2,:)=rangetest_count(mean(results_cpc));
    rangetest_count_results(3,:)=rangetest_count(mean(results_WER))

    % test for mutant pattern
    scm_results(1,:) = mean(results_GL2);
    scm_results(2,:) = mean(results_cpc);
    scm_results(3,:) = mean(results_WER)

end

%% plot

crossSection_plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%