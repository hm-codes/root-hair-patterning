%% t-test lit

function t_results = ttest_count(sim_means,sim_stds,sim_Number_of_repeats)

% literature values

HH_mean = 61.8;
HN_mean = 38.2;
NH_mean = 22.4;
NN_mean = 77.6;

H_std = 12.9;
N_std = 15.2;
Number_of_repeats = 11;

% simulation values

sim_HH_mean = sim_means(1);
sim_HN_mean = sim_means(2);
sim_NH_mean = sim_means(3);
sim_NN_mean = sim_means(4);

sim_H_std = sim_stds(1);
sim_N_std = sim_stds(4);

%% HH

% t-statistic

deltaX_HH = sim_HH_mean - HH_mean;

sX_H = H_std/sqrt(Number_of_repeats);
sim_sX_H = sim_H_std/sqrt(sim_Number_of_repeats);

t_HH = (deltaX_HH)/sqrt(sim_sX_H^2+sX_H^2);

% degrees of freedom
v_H =  (sim_sX_H^2+sX_H^2)^2/(H_std^4/(Number_of_repeats^2*(Number_of_repeats-1)) +...
                        sim_H_std^4/(sim_Number_of_repeats^2*(sim_Number_of_repeats-1)));

% check tables
t_critical = [12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201,...
               2.179, 2.16, 2.145, 2.131, 2.12, 2.11, 2.101, 2.093, 2.086, 2.08, 2.074,...
               2.069, 2.064, 2.06, 2.056, 2.052, 2.048, 2.045, 2.042, 1.96];

if v_H >= 32
    v_H = 31;
end

if t_HH < t_critical(floor(v_H))
    t_results(1) = 1;
else
    t_results(1) = 0;
end

%% HN

% t-statistic

deltaX_HN = sim_HN_mean - HN_mean;

t_HN = (deltaX_HN)/sqrt(sim_sX_H^2+sX_H^2);

% check tables

if t_HN < t_critical(floor(v_H))
    t_results(2) = 1;
else
    t_results(2) = 0;
end

%% NH

% t-statistic

deltaX_NH = sim_NH_mean - NH_mean;

sX_N = N_std/sqrt(Number_of_repeats);
sim_sX_N = sim_N_std/sqrt(sim_Number_of_repeats);

t_NH = (deltaX_NH)/sqrt(sim_sX_N^2+sX_N^2);

% degrees of freedom
v_N =  (sim_sX_N^2+sX_N^2)^2/(N_std^4/(Number_of_repeats^2*(Number_of_repeats-1)) +...
                        sim_N_std^4/(sim_Number_of_repeats^2*(sim_Number_of_repeats-1)));

if v_N >= 32
    v_N = 31;
end

% check tables

if t_NH < t_critical(floor(v_N))
    t_results(3) = 1;
else
    t_results(3) = 0;
end

%% NN

% t-statistic

deltaX_NN = sim_NN_mean - NN_mean;

t_NN = (deltaX_NN)/sqrt(sim_sX_N^2+sX_N^2);

% check tables

if t_NN < t_critical(floor(v_N))
    t_results(4) = 1;
else
    t_results(4) = 0;
end


end
