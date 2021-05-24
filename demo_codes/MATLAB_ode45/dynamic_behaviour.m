clear all

%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 1000];

%% Starting the loop for a inital conditions:
for i=1:1:10
    A = 10*rand(1);
    B = 10*rand(1);    
%% Calling ODE function:

    [t, x] = ode45(@(t,x) interactions(t,x),domain,[A;B]);
    %Ploting:
    plot(t,x(:,1)); 
    hold on
    
end
%% Adding feature to the plot:
gcf % get current figure:
xlabel('Time(a.u)')
ylabel('Concentration(a.u)')
