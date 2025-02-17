%%
%3B
load("sensitivity_interval.mat")

g1 = repmat({'SCM+'},1,1);
g2 = repmat({'SCM-'},1,1);
g = [g1; g2];

figure(1)
b = bar(data');
b(1).FaceColor = [0.1 0.6 0.8];
b(2).FaceColor = [0.2 0.9 0.5];
legend(g)
xlabel('parameter set')
ylabel('total insensitive region')

%%
%3C
load("pars_ordered_by_sensitivity.mat")
figure(2)
b = bar(total_sensitivity);
b.FaceColor = [0.1 0.6 0.8];
xlabel('parameter')
ylabel('total sensitivity difference')
