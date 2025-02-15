%% parameters

%% parameter which change between models

% nonlinearities
if nonlinearity == "linear"
    n1 = 1;                      % number of GL3/EGL3 proteins in the activator complex
    n2 = 1;                      % number of WER/MYB23 proteins in the activator complex
    n3 = 1;                      % number of WER binding sites on the myb23 promoter
    n4 = 1;                      % number of WER binding sites on the cpc promoter
elseif nonlinearity == "myb"
    n1 = 1;                      % number of GL3/EGL3 proteins in the activator complex
    n2 = 1;                      % number of WER/MYB23 proteins in the activator complex
    n3 = parspace(j,45);         % number of WER binding sites on the myb23 promoter
    n4 = 1;                      % number of WER binding sites on the cpc promoter
elseif nonlinearity == "cpc"
    n1 = 1;                      % number of GL3/EGL3 proteins in the activator complex
    n2 = 1;                      % number of WER/MYB23 proteins in the activator complex
    n3 = 1;                      % number of WER binding sites on the myb23 promoter
    n4 = parspace(j,46);         % number of WER binding sites on the cpc promoter
elseif nonlinearity == "complex"
    n1 = parspace(j,43);         % number of GL3/EGL3 proteins in the activator complex
    n2 = parspace(j,44);         % number of WER/MYB23 proteins in the activator complex
    n3 = 1;                      % number of WER binding sites on the myb23 promoter
    n4 = 1;                      % number of WER binding sites on the cpc promoter
elseif nonlinearity == "mybcpc"
    n1 = 1;                      % number of GL3/EGL3 proteins in the activator complex
    n2 = 1;                      % number of WER/MYB23 proteins in the activator complex
    n3 = parspace(j,45);         % number of WER binding sites on the myb23 promoter
    n4 = parspace(j,46);         % number of WER binding sites on the cpc promoter
elseif nonlinearity == "mybcomplex"
    n1 = parspace(j,43);         % number of GL3/EGL3 proteins in the activator complex
    n2 = parspace(j,44);         % number of WER/MYB23 proteins in the activator complex
    n3 = parspace(j,45);         % number of WER binding sites on the myb23 promoter
    n4 = 1;                      % number of WER binding sites on the cpc promoter
elseif nonlinearity == "cpccomplex"
    n1 = parspace(j,43);         % number of GL3/EGL3 proteins in the activator complex
    n2 = parspace(j,44);         % number of WER/MYB23 proteins in the activator complex
    n3 = 1;                      % number of WER binding sites on the myb23 promoter
    n4 = parspace(j,46);         % number of WER binding sites on the cpc promoter
elseif nonlinearity == "all"
    n1 = parspace(j,43);         % number of GL3/EGL3 proteins in the activator complex
    n2 = parspace(j,44);         % number of WER/MYB23 proteins in the activator complex
    n3 = parspace(j,45);         % number of WER binding sites on the myb23 promoter
    n4 = parspace(j,46);         % number of WER binding sites on the cpc promoter
end

% directed movement
if movement == "DM"
    k3 = parspace(j,16);         % rate of SCM directed CPC movement
    rw21 = 0;                    % repression of wer by SCM
    rw22 = parspace(j,42);       % repression of wer by cortical signal
elseif movement == "Diff"
    k3 = 0;                      % rate of SCM directed CPC movement
    rw21 = parspace(j,41);       % repression of wer by SCM
    rw22 = 0;                    % repression of wer by cortical signal
end

% repression of wer
if werR == "complex"
    rw11 = 0;                    % repression of wer by CPC
    rw12 = parspace(j,40);       % repression of wer by the inhibitor complex
elseif werR == "protein"
    rw11 = parspace(j,39);       % repression of wer by CPC
    rw12 = 0;                    % repression of wer by the inhibitor complex
elseif werR == "CI"
    rw11 = 0;                    % repression of wer by CPC
    rw12 = 0;                    % repression of wer by the inhibitor complex
end

% EGL3 diffusion
if DEGL3 == "Y"
    DE = parspace(j,14);         % rate of diffusion of EGL3
elseif DEGL3 == "N"
    DE = 0;                      % rate of diffusion of EGL3
end

%% static parameters

% reaction rates
k11 = parspace(j,1);         % A_GW association rate
k12 = parspace(j,2);         % A_GW disassociation rate
k13 = parspace(j,3);         % A_EW association rate
k14 = parspace(j,4);         % A_EW disassociation rate
k15 = parspace(j,5);         % A_GM association rate
k16 = parspace(j,6);         % A_GM disassociation rate
k17 = parspace(j,7);         % A_EM association rate
k18 = parspace(j,8);         % A_EM disassociation rate
k21 = parspace(j,9);         % I_GC association rate
k22 = parspace(j,10);        % I_GC disassociation rate
k23 = parspace(j,11);        % I_EC association rate
k24 = parspace(j,12);        % I_EC disassociation rate

% diffusion
DG = parspace(j,13);         % rate of diffusion of GL3
DC = parspace(j,15);         % rate of diffusion of CPC

% basal production
bg = parspace(j,17);         % basal production rate of gl3 mRNA
be = parspace(j,18);         % basal production rate of egl3 mRNA
bw = parspace(j,19);         % basal production rate of wer mRNA

% degradation
dG = parspace(j,20);         % degradation rate of GL3
dE = parspace(j,21);         % degradation rate of EGL3
dC = parspace(j,22);         % degradation rate of CPC
dW = parspace(j,23);         % degradation rate of WER
dM = parspace(j,24);         % degradation rate of MYB23
dg = parspace(j,25);         % degradation rate of gl3 mRNA
de = parspace(j,26);         % degradation rate of egl3 mRNA
dc = parspace(j,27);         % degradation rate of cpc mRNA
dw = parspace(j,28);         % degradation rate of wer mRNA
dm = parspace(j,29);         % degradation rate of myb23 mRNA

% translation
qG = parspace(j,30);         % rate of translation of GL3
qE = parspace(j,31);         % rate of translation of EGL3
qC = parspace(j,32);         % rate of translation of CPC
qW = parspace(j,33);         % rate of translation of WER
qM = parspace(j,34);         % rate of translation of MYB23

% promotion
pc = parspace(j,35);         % production rate of cpc mRNA
pm = parspace(j,36);         % production rate of myb23 mRNA

% promotion
rg = parspace(j,37);         % repression of gl3 by the activator complex
re = parspace(j,38);         % repression of egl3 by the activator complex

% michaelis menten terms
kdM = parspace(j,47);         % effective dissociation constant of WERc on myb23 promoter
kdC = parspace(j,48);         % effective dissociation constant of WERc on myb23 promoter