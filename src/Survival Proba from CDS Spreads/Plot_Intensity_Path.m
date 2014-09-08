

load('CDX9_17_Dec_2007.mat');
NbPath = 20;

%% Parametre du CIR
c = 10;
lambda = 0.1;

freq_saut = c;
taille_saut = 1/lambda;
a=0.2;
sigma = 0.0001

% 
% eta = Index_data.CIR_Param1_a
% nu = Index_data.CIR_Param1_c
% Eta =  Index_data.Eta_CIR1(1,:);

% eta = Index_data.CIR_Param0_a
% nu = Index_data.CIR_Param0_c
% Eta =  Index_data.Eta_CIR0(5,:); %5=AIG

%piecewise-constant mean-reversion parameter
% mu1 = Eta(1,1)
% mu2 = Eta(1,2)
% mu3 = Eta(1,3)
%mu3 = Eta(1,3)


% mu1 = 0.0079138;
% mu2 =  0.0058239 ;
% mu3 =0.0062575 ;

% AIG
% 0.0201   -0.0060    0.0172    0.0009


b1 =  0.00027068;
b2 =  -0.025811 ;
b3 =  -0.0025663;

% CCR-HomeLoans
%0.4397   -0.5146    0.5210   -0.1760

% b1 =  0.4397;
% b2 = -0.5146 ;
% b3 =  0.5210;

x =  0;

T1 = 3;
T2 = 5;
T3 = 7;
%T4 = 7;
t = T3;


ColorSet = varycolor(NbPath);
set(gca, 'ColorOrder', ColorSet);

hold all;
for m = 1:NbPath


%% Simulation du CIR (Euler Scheme)
Nb_time_step = 1000;
Step = (t-0)/Nb_time_step;
T = 0+(0:Nb_time_step)'.*Step;
X = zeros(Nb_time_step+1, 1);
X(1) = x;
Integral_Approx = zeros(Nb_time_step+1, 1);

for i = 1:Nb_time_step
   if (T(i)<T1) 
       b = b1;
   else
       if (T(i)<T2)
            b = b2;
       else
            b = b3;
       end
   end
       
   X(i+1) = X(i) + a*(b - X(i))*Step + sigma*gamrnd(freq_saut*Step,taille_saut);
   Integral_Approx(i+1) = Integral_Approx(i) + X(i+1)*Step;
   
end

plot(T, X);
title('Intensity Path AIG')
xlabel('maturity')
ylabel('Intensity')
end

%plot(T, Integral_Approx, 'r');


