load("a.mat")

figure;hold

b1 = bar(1,a(1))
b2 = bar(3:5,a(2:4))
b3 = bar(7:9,a(5:7))
b4 = bar(11,a(8))

b1.FaceColor =  [0.1 0.6 0.8];
b2.FaceColor = [0.2 0.9 0.5];
b3.FaceColor = [236 129 132]/255;

plot([2 2],[0 14],'color', [.5 .5 .5],LineWidth=1)
plot([6 6],[0 14],'color', [.5 .5 .5],LineWidth=1)
plot([10 10],[0 14],'color', [.5 .5 .5],LineWidth=1)

plot([0 12],[10 10],'Color',[0.1 0.6 0.8],'LineStyle','--',LineWidth=2)

xticks([1 3 4 5 8 9])
ylabel('number of successful parameter sets')
xlabel('model')