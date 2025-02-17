%% t-test lit

function t_results = rangetest_count_WT(sim_means)

% literature values

HH_upper = 100;
HH_lower = 75;

NN_upper = 100;
NN_lower = 90.9;


% simulation values

sim_HH_mean = sim_means(1);
sim_NN_mean = sim_means(4);


%% HH

if sim_HH_mean <= HH_upper && sim_HH_mean >= HH_lower
    t_results(1:2) = 1;
else
    t_results(1:2) = 0;
end

%% NN

if sim_NN_mean <= NN_upper && sim_NN_mean >= NN_lower
    t_results(3:4) = 1;
else
    t_results(3:4) = 0;
end


end
