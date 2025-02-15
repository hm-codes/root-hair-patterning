%% parameters

model_parameters

%% check stability condition r < 0.5

alpha = max([DG,DC,DE]);     % maximum diffusion coefficient
r = alpha*dt/(dx^2);

if r <= 0.5
    %msgbox(sprintf('Stability condition satisfied, r = %f',r)) % print r
else
    while r >= 0.5                              % continue until stability condition is met
        NT = NT + 10;                           % increase number of timesteps within Time
        dt = Time/(NT-1);                       % new smaller timestep
        r = alpha*dt/(dx^2);                    % new r
    end
    %msgbox(sprintf('New r = %f, dt = %d',r,dt)) % print new r and new dt
end

%% initial conditions

rng(100)
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


%% simulation

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


%% TEST - ARE ALL COMPONENTS IN THE CORRECT PLACE

% find totals of each component
Totals = zeros(8,L);
Totals(1,:) = g;                                           % gl3
Totals(2,:) = e;                                           % egl3
Totals(3,:) = C + I_GC + I_EC;                             % total CPC
Totals(4,:) = E + A_EW + A_EM + I_EC;                      % total EGL3
Totals(5,:) = m;                                           % myb23
Totals(6,:) = w;                                           % wer
Totals(7,:) = W + A_GW + A_EW;                             % total WER
Totals(8,:) = M + A_GM + A_EM;                             % total MYB23
Totals(9,:) = c;                                           % cpc
Totals(10,:) = G + A_GW + A_GM + I_GC;                     % total GL3

% find average concentrations across all cells
Averages = mean(Totals,2);

numVars = 10;
cellSpec = zeros(numVars,L);

% decide if each cell is a H (1) or NH (0)
% according to each variable
for cell = 1:L                                      % each cell
    for v = 1:4                                     % varibles CPC, gl3c, egl3, and EGL3 which are expected to be predominately in H cells
        if Totals(v,cell) >= Averages(v)
            cellSpec(v,cell) = 1;                   % if concentration is greater than the average, we expect a H cell (1)
        else
            cellSpec(v,cell) = 0;
        end
    end
    for v = 5:10                                    % variables wer, WER, myb, MYB23, cpc and GL3c which are expected to be predominantly in NH cells
        if Totals(v,cell) > Averages(v)
            cellSpec(v,cell) = 0;                   % if concentration is greater than the average, we expect a NH cell (0)
        else
            cellSpec(v,cell) = 1;
        end
    end
end

Ha = find(sum(cellSpec) == numVars);                              % find indices of hair cells
NHa = find(sum(cellSpec) == 0);                                   % find indices of non-hair cells

patternW = cellSpec;
patternW(7,:) = [];                                               % remove WER row
Ha_mWER = find(sum(patternW(:,:)) == numVars-1);                  % find indices of hair cells w/o WER
NHa_mWER = find(sum(patternW(:,:)) == 0);                         % find indices of non-hair cells w/o WER

patternW(6,:) = [];                                               % remove wer row
Ha_mW = find(sum(patternW(:,:)) == numVars-2);                    % find indices of hair cells w/o WER and wer
NHa_mW = find(sum(patternW(:,:)) == 0);                           % find indices of non-hair cells w/o WER and wer

patternW = cellSpec;
patternW(6,:) = [];                                               % remove wer row
Ha_mwer = find(sum(patternW(:,:)) == numVars-1);                  % find indices of hair cells w/o wer
NHa_mwer = find(sum(patternW(:,:)) == 0);                         % find indices of non-hair cells w/o wer

patternC = cellSpec;
patternC(9,:) = [];                                               % remove cpc row
Ha_mcpc = find(sum(patternC(:,:)) == numVars-1);                  % find indices of hair cells w/o cpc
NHa_mcpc = find(sum(patternC(:,:)) == 0);                         % find indices of non-hair cells w/o cpc

patternC(3,:) = [];                                               % remove CPC row
Ha_mC = find(sum(patternC(:,:)) == numVars-2);                    % find indices of hair cells w/o CPC and cpc
NHa_mC = find(sum(patternC(:,:)) == 0);                           % find indices of non-hair cells w/o CPC and cpc

patternC = cellSpec;
patternC(3,:) = [];                                               % remove CPC row
Ha_mCPC = find(sum(patternC(:,:)) == numVars-1);                  % find indices of hair cells w/o CPC
NHa_mCPC = find(sum(patternC(:,:)) == 0);                         % find indices of non-hair cells w/o CPC

patternG = cellSpec;
patternG(10,:) = [];                                              % remove GL3c row
Ha_mGL3 = find(sum(patternG(:,:)) == numVars-1);                  % find indices of hair cells w/o GL3c
NHa_mGL3 = find(sum(patternG(:,:)) == 0);                         % find indices of non-hair cells w/o GL3c

patternG(1,:) = [];                                               % remove gl3c row
Ha_mG = find(sum(patternG(:,:)) == numVars-2);                    % find indices of hair cells w/o GL3c and gl3c
NHa_mG = find(sum(patternG(:,:)) == 0);                           % find indices of non-hair cells w/o GL3c and gl3c

patternG = cellSpec;
patternG(1,:) = [];                                               % remove gl3c row
Ha_mgl3c = find(sum(patternG(:,:)) == numVars-1);                 % find indices of hair cells w/o gl3c
NHa_mgl3c = find(sum(patternG(:,:)) == 0);                        % find indices of non-hair cells w/o gl3c

patternD = cellSpec;
patternD(10,:) = [];                                              % remove CPC as well as GL3c row
patternD(3,:) = [];                                               % remove CPC as well as GL3c row
Ha_mD = find(sum(patternD(:,:)) == numVars-2);                    % find indices of hair cells w/o CPC and GL3c
NHa_mD = find(sum(patternD(:,:)) == 0);                           % find indices of non-hair cells w/o CPC and GL3c

patternE = cellSpec;
patternE(4,:) = [];                                               % remove EGL3 row
Ha_mEGL3 = find(sum(patternE(:,:)) == numVars-1);                 % find indices of hair cells w/o EGL3
NHa_mEGL3 = find(sum(patternE(:,:)) == 0);                        % find indices of non-hair cells w/o EGL3

patternE(2,:) = [];                                               % remove egl3 row
Ha_mE = find(sum(patternE(:,:)) == numVars-2);                    % find indices of hair cells w/o EGL3 and egl3
NHa_mE = find(sum(patternE(:,:)) == 0);                           % find indices of non-hair cells w/o EGL3 and egl3

patternE = cellSpec;
patternE(2,:) = [];                                               % remove egl3 row
Ha_megl3 = find(sum(patternE(:,:)) == numVars-1);                 % find indices of hair cells w/o egl3
NHa_megl3 = find(sum(patternE(:,:)) == 0);                        % find indices of non-hair cells w/o egl3

patternM = cellSpec;
patternM(8,:) = [];                                               % remove MYB23 row
Ha_mMYB23 = find(sum(patternM(:,:)) == numVars-1);                % find indices of hair cells w/o MYB23
NHa_mMYB23 = find(sum(patternM(:,:)) == 0);                       % find indices of non-hair cells w/o MYB23

patternM(5,:) = [];                                               % remove myb23 row
Ha_mM = find(sum(patternM(:,:)) == numVars-2);                    % find indices of hair cells w/o MYB23 and myb23
NHa_mM = find(sum(patternM(:,:)) == 0);                           % find indices of non-hair cells w/o MYB23 and myb23

patternM = cellSpec;
patternM(5,:) = [];                                               % remove myb23 row
Ha_mmyb23 = find(sum(patternM(:,:)) == numVars-1);                % find indices of hair cells w/o myb23
NHa_mmyb23 = find(sum(patternM(:,:)) == 0);                       % find indices of non-hair cells w/o myb23


% record pattern type
if length(Ha)==length(H) && length(NHa)==length(NH) && all(H == Ha) && all(NH == NHa)
    "pattern_CS"
    result = 1;
elseif length(Ha) + length(NHa) == L && ~isempty(Ha) && ~isempty(NHa)
    "pattern_yes"
    result = 2;
elseif ~isempty(Ha) && ~isempty(NHa)
    "pattern_part"
    result = 2.1;
elseif length(Ha_mWER) + length(NHa_mWER) == L && ~isempty(Ha_mWER) && ~isempty(NHa_mWER)
    "pattern_mWER"
    result = 3.1;
elseif length(Ha_mwer) + length(NHa_mwer) == L && ~isempty(Ha_mwer) && ~isempty(NHa_mwer)
    "pattern_mwer"
    result = 3.2;
elseif length(Ha_mW) + length(NHa_mW) == L && ~isempty(Ha_mW) && ~isempty(NHa_mW)
    "pattern_mW"
    result = 3;
elseif length(Ha_mCPC) + length(NHa_mCPC) == L && ~isempty(Ha_mCPC) && ~isempty(NHa_mCPC)
    "pattern_mCPC"
    result = 4.1;
elseif length(Ha_mcpc) + length(NHa_mcpc) == L && ~isempty(Ha_mcpc) && ~isempty(NHa_mcpc)
    "pattern_mcpc"
    result = 4.2;
elseif length(Ha_mC) + length(NHa_mC) == L && ~isempty(Ha_mC) && ~isempty(NHa_mC)
    "pattern_mC"
    result = 4;
elseif length(Ha_mGL3) + length(NHa_mGL3) == L && ~isempty(Ha_mGL3) && ~isempty(NHa_mGL3)
    "pattern_mGL3c"
    result = 5.1;
elseif length(Ha_mgl3c) + length(NHa_mgl3c) == L && ~isempty(Ha_mgl3c) && ~isempty(NHa_mgl3c)
    "pattern_mgl3c"
    result = 5.2;
elseif length(Ha_mG) + length(NHa_mG) == L && ~isempty(Ha_mG) && ~isempty(NHa_mG)
    "pattern_mG"
    result = 5;
elseif length(Ha_mMYB23) + length(NHa_mMYB23) == L && ~isempty(Ha_mMYB23) && ~isempty(NHa_mMYB23)
    "pattern_mMYB23"
    result = 6.1;
elseif length(Ha_mmyb23) + length(NHa_mmyb23) == L && ~isempty(Ha_mmyb23) && ~isempty(NHa_mmyb23)
    "pattern_mmyb23"
    result = 6.2;
elseif length(Ha_mM) + length(NHa_mM) == L && ~isempty(Ha_mM) && ~isempty(NHa_mM)
    "pattern_mM"
    result = 6;
elseif length(Ha_mEGL3) + length(NHa_mEGL3) == L && ~isempty(Ha_mEGL3) && ~isempty(NHa_mEGL3)
    "pattern_mEGL3c"
    result = 7.1;
elseif length(Ha_megl3) + length(NHa_megl3) == L && ~isempty(Ha_megl3) && ~isempty(NHa_megl3)
    "pattern_megl3c"
    result = 7.2;
elseif length(Ha_mE) + length(NHa_mE) == L && ~isempty(Ha_mE) && ~isempty(NHa_mE)
    "pattern_mE"
    result = 7;
elseif length(Ha_mD) + length(NHa_mD) == L && ~isempty(Ha_mD) && ~isempty(NHa_mD)
    "pattern_mD"
    result = 8;
elseif length(Ha) == L || length(NHa) == L
    "homogeneous"
    result = 9;
else
    "no pattern"
    result = 10;
end

%% check thresholds

th = 1.1;
if any(max(Totals,[],2)./min(Totals,[],2) <= th)
    thresh = 0;     % does not pass threshold test
else
    thresh = 1;     % passes threshold test
end

%% collect results
fullResults(j,:) = [result,thresh,steadyState,broken];

% cellspec
fullCellSpec(j,:,:) = cellSpec;
sumCellSpec(j,:) = squeeze(sum(cellSpec));

% save results every 1000 parsets
if mod(j,1000) == 0
    filename = sprintf('full_results_wt_%s_%s_%s_%s.mat',nonlinearity,movement,werR,DEGL3);
    save(filename,'fullResults','fullCellSpec','sumCellSpec','-v7.3')
end

% save final concentrations
% filename = sprintf('./Data_wt_%s_%s_%s_%s/finalConc%d.mat',nonlinearity,movement,werR,DEGL3,j);
% save(filename,'G_new','E_new','C_new','W_new','M_new','g_new',...
%     'e_new','c_new','w_new','m_new','A_GW_new','A_EW_new','A_GM_new',...
%     'A_EM_new','I_GC_new','I_EC_new','G','E','C','W','M','g','e','c',...
%     'w','m','A_GW','A_EW','A_GM','A_EM','I_GC','I_EC','-v7.3')
