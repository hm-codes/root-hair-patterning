%% loading data

load Dmax_actual.mat

% columns: 
% 1: A_GW	
% 2: A_EW	
% 3: A_GM	
% 4: A_EM	
% 5: I_GC	
% 6: I_EC


%% 220 data points murged

colors = flip([[150 190 231]./255; [134 13 44]/255; [236 129 132]./255; [106 114 74]/255; [0.1 0.6 0.8]; 0.2 0.9 0.5]);


figure
subplot(1,2,1); hold

boxplot(trich,'Colors','k');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:));
end
boxplot(trich,'Colors','k');
title('trichoblasts')
xticklabels({'WER/GL3','WER/EGL3','MYB23/GL3','MYB23/EGL3','CPC/GL3','CPC/EGL3'})
ylim([-1,50])

subplot(1,2,2);hold
boxplot(atrich,'Colors','k');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:));
end
boxplot(atrich,'Colors','k');
title('atrichoblasts')
xticklabels({'WER/GL3','WER/EGL3','MYB23/GL3','MYB23/EGL3','CPC/GL3','CPC/EGL3'})
ylim([-1,50])

