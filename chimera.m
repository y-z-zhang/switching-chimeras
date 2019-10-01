function chimera
%% This code simulates the switching chimera system.

fonttype = 'Times';
fsize = 23;
txtattrib = {'FontName',fonttype,'FontSize',fsize,...
         'FontWeight','normal'};
txtattrib2 = {txtattrib{:},'Interpreter','Latex'};

tableau20 = [ 31, 119, 180; 174, 199, 232; 255, 127,  14; 255, 187, 120;...    
              44, 160,  44; 152, 223, 138; 214,  39,  40; 255, 152, 150;...    
             148, 103, 189; 197, 176, 213; 140,  86,  75; 196, 156, 148;...    
             227, 119, 194; 247, 182, 210;  23, 190, 207; 158, 218, 229;...    
             188, 189,  34; 219, 219, 141; 127, 127, 127; 199, 199, 199]/255;

%% coupling strength
sigma = 1.7;

%% parameter for the logistic maps
r = 3;

%% simulation time
t = 1e6;

%% length of simulation data shown in plots
T = 3000;

%% noise intensity
eta = 1e-10;

%% adjaency matrix for the two-ring network
n = 6; %% size of each ring
adj = zeros(1,n);
adj(1,2) = 1; adj(1,end) = 1;
ring = adj;
for ii = 1:n-1
  adj = [adj;circshift(ring,[0,ii])];
end

adj = [adj,.2*ones(n); .2*ones(n),adj];

%% Laplacian matrix for the two-ring network
rowsum = sum(adj,2);
lap = diag(rowsum) - adj;

%% initial conditions for the oscillators
x(1:2*n,1) = rand(2*n,1);

%% uncomment below if want to filter out short-wavelength components from noise
%swb1 = [1,-1,1,-1,1,-1,0,0,0,0,0,0]';
%swb1 = swb1/norm(swb1);
%swb2 = [0,0,0,0,0,0,1,-1,1,-1,1,-1]';
%swb2 = swb2/norm(swb2);

%% evolve the system for t iterations
for ih = 2:t
	noise = normrnd(0,eta,2*n,1);
	%% uncomment below if want to filter out short-wavelength components from noise
	%noise = noise - dot(swb1,noise)*swb1 - dot(swb2,noise)*swb2;
	for i = 1:2*n
  		x(i,ih) = f(i,ih);
  	end
end

%% dynamical equation describing coupled logistic maps under noise
function z = f(i,ih)
	Z = r*logistic(x(i,ih-1)) - sigma*lap(i,:)*logistic(x(:,ih-1)) + noise(i);
	z = mod(Z, 1);
end

function y = logistic(x)
	y = x.*(1-x);
end

%% sync error in the first and the second ring
syncerr1 = std(x(1:n,:));
syncerr2 = std(x(n+1:2*n,:));

%% identify switching by monitoring and comparing sync errors
interval = 50;
err1 = zeros(t/interval,1);
err2 = zeros(t/interval,1);
for ii = 1:t/interval
	err1(ii) = mean(syncerr1((ii-1)*interval+1:ii*interval));
	err2(ii) = mean(syncerr2((ii-1)*interval+1:ii*interval));
end
err = err1 > err2;
time = [];
for ii = 2:t/interval-1
	if err(ii) ~= err(ii+1)
		time = [time,ii];
	end
end

%% plotting the first and last T iterations of simulation data
set(0,'DefaultAxesFontSize',15)

figure
hAxis(1)=subplot(3,1,1);
pos = get( hAxis(1), 'Position');
pos(1)=.11;
pos(2)=.68;
pos(3)=.87;
pos(4)=.28;
set(hAxis(1), 'Position', pos);
hold on
for i = 1:n
	plot(1:t,x(i,:),'.','MarkerSize',2,'color',tableau20(13-2*i,:));
end
set(gca,'xtick',[]);
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
xlim([1 T])
box on;
hold off;

hAxis(2)=subplot(3,1,2);
pos = get( hAxis(2), 'Position');
pos(1)=.11;
pos(2)=.36;
pos(3)=.87;
pos(4)=.28;
set(hAxis(2), 'Position', pos);
semilogy(1:t,syncerr1,'Color',tableau20(1,:),'LineWidth',2);
hold on
semilogy(1:t,syncerr2,'Color',tableau20(3,:),'LineWidth',2);
hold off
box on
ylim([1e-10,1e1])
xlim([1 T])
set(gca,'YTick',[1e-9 1e-5 1e-1])
set(gca,'xtick',[]);

hAxis(3)=subplot(3,1,3);
pos = get( hAxis(3), 'Position');
pos(1)=.11;
pos(2)=.035;
pos(3)=.87;
pos(4)=.28;
set(hAxis(3), 'Position', pos);
hold on
for i = 1:n
	plot(1:t,x(n+i,:),'.','MarkerSize',2,'color',tableau20(15-2*i,:));
end
xlim([1 T])
set(gca,'xtick',[]);
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
box on;
hold off;

set(gcf, 'PaperPosition', [0 0 5 3]);
set(gcf, 'PaperSize', [5 3]);
saveas(gcf,'chimera_1.pdf');



figure
hAxis(1)=subplot(3,1,1);
pos = get( hAxis(1), 'Position');
pos(1)=.01;
pos(2)=.68;
pos(3)=.97;
pos(4)=.28;
set(hAxis(1), 'Position', pos);
hold on
for i = 1:n
	plot(1:t,x(i,:),'.','MarkerSize',2,'color',tableau20(13-2*i,:));
end
set(gca,'xtick',[]);
xlim([t-T t])
box on;
hold off;

hAxis(2)=subplot(3,1,2);
pos = get( hAxis(2), 'Position');
pos(1)=.01;
pos(2)=.36;
pos(3)=.97;
pos(4)=.28;
set(hAxis(2), 'Position', pos);
semilogy(1:t,syncerr1,'Color',tableau20(1,:),'LineWidth',2);
hold on
semilogy(1:t,syncerr2,'Color',tableau20(3,:),'LineWidth',2);
hold off
box on
ylim([1e-10,1e1])
xlim([t-T t])
set(gca,'YTick',[1e-9 1e-5 1e-1])
set(gca,'xtick',[]);

hAxis(3)=subplot(3,1,3);
pos = get( hAxis(3), 'Position');
pos(1)=.01;
pos(2)=.04;
pos(3)=.97;
pos(4)=.28;
set(hAxis(3), 'Position', pos);
hold on
for i = 1:n
	plot(1:t,x(n+i,:),'.','MarkerSize',2,'color',tableau20(15-2*i,:));
end
xlim([t-T t])
set(gca,'xtick',[]);
box on;
hold off;

set(gcf, 'PaperPosition', [0 0 5 3]);
set(gcf, 'PaperSize', [5 3]);
saveas(gcf,'chimera_2.pdf');

end