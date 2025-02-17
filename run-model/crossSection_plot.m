%% plot cross section

chosen_row = 20;
L = length(concentrations(1,1,:));

GL2_T = squeeze(sum(concentrations(chosen_row,11:14,:)));
GL2_T(end+1,:) = GL2_T(end);
plot_WER = ones(2,L+1)*max(GL2_T);
plot_WER(1,:) = GL2_T;

n = length(plot_WER(:,1))-1;
r = (0:n)'/n+3;
theta = flip(8*pi/8:2*pi/L:2*pi+8*pi/8);
X = r*cos(theta);
Y = r*sin(theta);

CS = corticalClefts=="H";
CS(end+1) = CS(end);
plot_CS(1,:) = CS;
plot_CS(2,:) = CS;

n2 = length(plot_CS(:,1))-1;
r2 = (0:n2)'/n2+2;
theta2 = (theta-pi/L).*CS;
X2 = r2*cos(theta2);
Y2 = r2*sin(theta2);
X2(:,find(corticalClefts=="N"))=[];
Y2(:,find(corticalClefts=="N"))=[];
X2(:,end) = X2(:,1);
Y2(:,end) = Y2(:,1);

figure(1)
hold on
map = [24 21 21
    30 46 46
    26 62 61
   24 77 75 
   29 96 94
30 109 108
35 119 118
32 140 140 
    34 148 146
    41 168 166
37 197 194]./255;
colormap(map);
h=pcolor(X,Y,plot_WER);
h.EdgeColor = [238 0 20]./255;
h.LineWidth = 5;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
axis equal tight 
grid off
set(gca,'Color','k')
set(gcf,'Color','k')
set(gcf, 'InvertHardCopy', 'off'); 
c = colorbar;
colorTitleHandle = get(c,'Title');
titleString = 'GL2';
set(colorTitleHandle ,'String',titleString,'Fontsize',10.5,'Color',"white");
c.Color="white";
c.FontSize = 10.5;

h2=pcolor(X2,Y2,ones(2,length(X2)));
h2.EdgeColor = [238 0 20]./255;
h2.LineWidth = 5;

%% plot timeseries

ind = find(squeeze(timecourse(chosen_row,1,:))==0);
tcplot = squeeze(timecourse(chosen_row,:,1:min(ind)-1))';
tcplot(:,L+1) = tcplot(:,1);


figure(2)
colormap(map)
pcolor(tcplot)
shading flat
set(gcf,'Color','k')
set(gcf, 'InvertHardCopy', 'off'); 
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

% c = colorbar;
% colorTitleHandle = get(c,'Title');
% titleString = 'GL2';
% set(colorTitleHandle ,'String',titleString,'Fontsize',10.5,'Color',"white");
% c.Color="white";
% c.FontSize = 10.5;
xlabel('cell',FontSize=10.5)
ylabel('time (ms)',FontSize=10.5)

hold on
x = (1:1:L);
y1 = zeros(length(x),1)';
y2 = ones(length(x),1)'*min(ind)-1;

plot([x;x],[y1;y2],'r',LineWidth=2)