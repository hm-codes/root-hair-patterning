load("p_success.mat")
load("scm_data_HH.mat")
load("WT_data_HH.mat")
load("scm_data_NN.mat")
load("WT_data_NN.mat")

scm_data_HH = [scm_data_HH ; nan(9,1)];
scm_data_NN = [scm_data_NN' ; nan(9,1)];


x = 1:11;

%% WT 
figure(1)

subplot(2,2,1);hold

a = [ WT_data_HH [p3423(:,5) ; nan(2,1)]  [p3905(:,5); nan(2,1)] ...
    [p4209(:,5); nan(2,1)] [p4634(:,5); nan(2,1)] [p6365(:,5); nan(2,1)] ... 
    [p10330(:,5); nan(2,1)] [p15122(:,5); nan(2,1)] [p15348(:,5); nan(2,1)] ...
    [p18537(:,5); nan(2,1)] [p18878(:,5); nan(2,1)]];

% HH
boxplot(a);

colors = [[0.1 0.6 0.8].*ones(10,1) ; 0.2 0.9 0.5 ];
title('wt HH')

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:));
end
boxplot(a);

subplot(2,2,2);hold
% NN
boxplot([ WT_data_NN [p3423(:,8) ; nan(2,1)]  [p3905(:,8); nan(2,1)] ...
    [p4209(:,8); nan(2,1)] [p4634(:,8); nan(2,1)] [p6365(:,8); nan(2,1)] ... 
    [p10330(:,8); nan(2,1)] [p15122(:,8); nan(2,1)] [p15348(:,8); nan(2,1)] ...
    [p18537(:,8); nan(2,1)] [p18878(:,8); nan(2,1)]],...
    'Notch','off')

x = 1:11;
colors = [[0.1 0.6 0.8].*ones(10,1) ; 0.2 0.9 0.5 ];

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:));
end

title('wt NN')
boxplot([ WT_data_NN [p3423(:,8) ; nan(2,1)]  [p3905(:,8); nan(2,1)] ...
    [p4209(:,8); nan(2,1)] [p4634(:,8); nan(2,1)] [p6365(:,8); nan(2,1)] ... 
    [p10330(:,8); nan(2,1)] [p15122(:,8); nan(2,1)] [p15348(:,8); nan(2,1)] ...
    [p18537(:,8); nan(2,1)] [p18878(:,8); nan(2,1)]],...
    'Notch','off')

%% scm mutant
subplot(2,2,3);hold

% HH
boxplot([ scm_data_HH p3423(:,1) p3905(:,1) p4209(:,1) p4634(:,1) p6365(:,1) ... 
    p10330(:,1) p15122(:,1) p15348(:,1) p18537(:,1) p18878(:,1)] )

x = 1:11;
colors = [[0.1 0.6 0.8].*ones(10,1) ; 0.2 0.9 0.5 ];

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:));
end

boxplot([ scm_data_HH p3423(:,1) p3905(:,1) p4209(:,1) p4634(:,1) p6365(:,1) ... 
    p10330(:,1) p15122(:,1) p15348(:,1) p18537(:,1) p18878(:,1)] )

title('scm HH')

subplot(2,2,4);hold

% NN
boxplot([ scm_data_NN p3423(:,4) p3905(:,4) p4209(:,4) p4634(:,4) p6365(:,4) ... 
    p10330(:,4) p15122(:,4) p15348(:,4) p18537(:,4) p18878(:,4)],...
    'Notch','off')
title('scm NN')

x = 1:11;
colors = [[0.1 0.6 0.8].*ones(10,1) ; 0.2 0.9 0.5 ];

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:));
end
boxplot([ scm_data_NN p3423(:,4) p3905(:,4) p4209(:,4) p4634(:,4) p6365(:,4) ... 
    p10330(:,4) p15122(:,4) p15348(:,4) p18537(:,4) p18878(:,4)],...
    'Notch','off')
title('scm HH')
