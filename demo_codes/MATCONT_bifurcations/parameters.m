function [x,v,s,h,f] = parameters

curdir = pwd;
init;
cd(curdir);

opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50000);
opt=contset(opt,'MinStepsize',0.1);
opt=contset(opt,'MaxStepsize',100);
opt=contset(opt,'Eigenvalues',1);


%% Paramaters
alpha = 10;
alpha0 = 10; 
beta = 0.5;  
n = 2;

ap = 1; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=1, thus bifurcation parameter is s (SNAIL levels)
handles = feval(@interactions);
tspan = 0:1:1000;

% initial condition
x_start = [6,5,7,8,9,10];

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,kmrgd)handles{2}(t, kmrgd, alpha, alpha0, beta, n),tspan,x_start);
x_init = x_time(end,:)';

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@interactions,x_init,[alpha, alpha0, beta, n],ap);
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);

end


