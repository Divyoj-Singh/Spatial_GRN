
%% Time Domain: 
% mention the time domain for which u want to run the ODE. 
domain = [0 500];

high=0;
low=0;
inter=0;

g=10;
A0B=2;
gamma_A = 0.5; 

% Transcription rate/ Production rate:
g_A = 10; 



g_divs=10;
g_array=linspace(0.7*g,1.3*g,g_divs);
A0B_divs=15;
A0B_array=linspace(0.7*A0B,1.4*A0B,A0B_divs);

%% Starting the loop for a inital conditions:
for i=1:1:200
    A = 10*rand(1); 
    B = 10*rand(1); 
    %% In general, you can generate N random numbers in the 
    %% interval (a,b) with the formula r = a + (b-a).*rand(N,1).
    
%% Calling ODE function:

    [t, x] = ode45(@(t,x) interactions(t,x),domain,[A;B]);
    %Ploting:
    plot(t,x(:,1)/(g_A/gamma_A), 'LineWidth', 2, 'DisplayName', 'A'); 
    %plot(x(:,1),x(:,2)); 
    hold on;
    %plot(t,x(:,2), 'LineWidth', 2, 'DisplayName', 'B');
    %hold on
    
    if x(end)>=150
        high=high+1;
    
    elseif x(end) <=60
        low=low+1;
    else 
        inter=inter+1;
    end

end
combine=low+high;
disp("High + Low= "+combine+', Hybrid= '+inter);
%% Adding feature to the plot:
% gcf % get current figure:
xlim([0 100]);
ylim([0 3.5]);
xlabel('Time (a.u)',"FontSize",12,"FontName",'Arial','FontWeight','bold');
ylabel('Concentration A(a.u)',"FontSize",12,"FontName",'Arial','FontWeight','bold');

str = {
'g_A = 5',   
'g_B = 5',
'\gamma_A = 0.1',
'\gamma_B = 0.1',  
'A_0B = 120',
'B_0A = 120', 
'A_0A = 80',   
'B_0B = 80',
'n_{AB} = 1',
'n_{BA} = 1',
'n_{AA} = 4',   
'n_{BB} = 4',
'\lambda_{AB} = 0.1',
'\lambda_{BA} = 0.1',
'\lambda_{AA} = 5',  
'\lambda_{BB} = 5'};

text(900,130,str,'FontSize',8);
% 
% saveas(gcf,'toggle_w_selfAct_4.png')

hold off 

for i=1:1:50
    A = 10*rand(1); 
    B = 10*rand(1); 
    %% In general, you can generate N random numbers in the 
    %% interval (a,b) with the formula r = a + (b-a).*rand(N,1).
    
%% Calling ODE function:

    [t, x] = ode45(@(t,x) interactions(t,x),domain,[A;B]);
    %Ploting:
    plot(t,x(:,2), 'LineWidth', 2, 'DisplayName', 'B'); 
    %plot(x(:,1),x(:,2)); 
    hold on;
    %plot(t,x(:,2), 'LineWidth', 2, 'DisplayName', 'B');
    %hold on
    
end
xlim([0 500]);
ylim([0 20]);
xlabel('Time (a.u)',"FontSize",12,"FontName",'Arial','FontWeight','bold');
ylabel('Concentration B (a.u)',"FontSize",12,"FontName",'Arial','FontWeight','bold');



Defining the Local functions used for the plotting
% hill function

function H = hill(X,X0,lambda,n)
% H+ 
H1 = (X^n)/(X0^n+X^n);
% H-
H2 = (X0^n)/(X0^n+X^n);

% Hill function =  H- plus (lambda)* H+
H = H2 + lambda*H1;

end


function dxdt = interactions(t,x)

% defining array of values for the time t: 
dxdt = zeros(2,1);

%% Paramaters
% Degradation rate:
gamma_A = 0.5;   gamma_B = 0.5;  

% Transcription rate/ Production rate:
g_A = 10; 
g_B = 10;   

% Hills function threshold :
A0B = 10;
B0A = 10;
A0A = 10;   B0B = 10;

% Cooperativity/ hill function coefficient:
nAtoB = 2;   nBtoA = 2;
nAtoA = 3;   nBtoB = 3;

% fold change/ lambda
lambda_AtoB = 0.1;   lambda_BtoA = 0.1;
lambda_AtoA = 3;   lambda_BtoB = 3;

% %% equations:
dxdt(1) = g_A*hill(x(1),A0A,lambda_AtoA,nAtoA)*hill(x(2),B0A,lambda_BtoA,nBtoA) - gamma_A*x(1);  %% equation A
dxdt(2) = g_B*hill(x(2),B0B,lambda_BtoB,nBtoB)*hill(x(1),A0B,lambda_AtoB,nAtoB) - gamma_B*x(2); % equation B
end
% % equations:
% dxdt(1) = g_A*hill(x(2),B0A,lambda_BtoA,nBtoA) - gamma_A*x(1) ;  %% equation A
% dxdt(2) = g_B*hill(x(1),A0B,lambda_AtoB,nAtoB) - gamma_B*x(2) ; % equation B
% end
