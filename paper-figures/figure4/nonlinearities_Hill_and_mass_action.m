%% 1D x^n1

x = [0:0.001:3];

colours = [0.1 0.6 0.8;
            0.2 0.9 0.5;
             [236 129 132]./255;
                [150 190 231]./255];

figure;
subplot(1,2,1);hold
for n = [1,2,3,4]

    plot(x,(x.^n)./(1+(x.^n)),'Color',colours(n,:),'LineWidth',4)
   
end
% plot(x,x)
ylim([0 1])
xlabel('transcription factor concentration')
ylabel('transcription rate')
legend('n=1','n=2','n=3','n=4',"Location","eastoutside")

x = [0:0.001:3];

subplot(1,2,2);hold
for n = [1,2,3,4]
    plot(x,x.^n,'Color',colours(n,:),'LineWidth',4)
end
xlabel('protein concentration')
ylabel('(protein concentration)^n')
ylim([0 3])
legend('n=1','n=2','n=3','n=4',"Location","eastoutside")