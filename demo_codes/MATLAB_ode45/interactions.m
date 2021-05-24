function dxdt = interactions(t,x)
% defining array of values for the time t: 
dxdt = zeros(2,1);

%% Paramaters
% Degradation rate:
gamma_A = 0.5;   gamma_B = 0.5;  

% Transcription rate/ Production rate:
g_A = 10;   g_B = 10;   

% Hills function threshold :
A0B = 1;   B0A = 1;   

% Cooperativity/ hill function coefficient:
nAtoB = 2;   nBtoA = 2;

% fold change/ lambda
lambda_AtoB =0.1;   lambda_BtoA = 0.1;


%% equations:
dxdt(1) = g_A*hill(x(2),B0A,lambda_BtoA,nBtoA) -gamma_A*x(1);  %% equation A
dxdt(2) = g_B*hill(x(1),A0B,lambda_AtoB,nAtoB)-gamma_B*x(2); % equation B
end
