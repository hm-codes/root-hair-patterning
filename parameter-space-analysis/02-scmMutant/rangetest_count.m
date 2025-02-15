%% t-test lit

function t_results = rangetest_count(sim_means)

% literature values

HH_upper = 74.7;
HH_lower = 48.9;

NN_upper = 92.8;
NN_lower = 62.4;


% simulation values

sim_HH_mean = sim_means(1);
sim_NN_mean = sim_means(4);


%% HH

if sim_HH_mean < HH_upper && sim_HH_mean > HH_lower
    t_results(1:2) = 1;
else
    t_results(1:2) = 0;
end

%% NN

if sim_NN_mean < NN_upper && sim_NN_mean > NN_lower
    t_results(3:4) = 1;
else
    t_results(3:4) = 0;
end


end
