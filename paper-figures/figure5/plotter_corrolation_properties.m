a = [5    4   3    1   1 4   2    4   5   4    5];    % number of conditions met
b = [0.05 1.1 0.95 0.7 3 0.8 2.05 1.6 0.3 0.05 0.65]; % D_E max
n3 = [4 3 2 4 2 3 2 3 4 4 3];

figure;hold
plot(a,-b,'o','MarkerSize',10,'MarkerFaceColor',[0.1 0.6 0.8])
xlabel('number of conditions met')
ylabel('-\alpha')
box on

[R,P] = corrcoef(-a,b);

R(2)
P(2)
