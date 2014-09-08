clear all
clc
% global nu;
% global eta;

%% Parametre du CIR
%  x = 0.1264; %starting point of the OU process
%  eta = 1;
%  nu = 0.05;

%starting point of the OU process
x =  0.0052;
%vitesse retour moyenne
a = 2;
sigma = 0.5;
 
% parametre du Levy
c = 1;
lambda = 1;
 
% x = 0.1264; %starting point of the OU process
% eta = 1;
% nu = 0.05;
% 
% x = 0.1264; %starting point of the OU process
% eta = 1;
% nu = 0.05;

%e=1; 
%piecewise-constant mean-reversion parameter

b1 = -0.052233;
b2 = -0.037872 ;
b3 = -0.048131 ;

% b1 = 0.1134;
% b2 = 0.0218;
% b3 = 0.0814;

T1 = 3;
T2 = 5;
T3 = 7;


%Horizon time
t = 6;


%% Siblation de Levy OU (Euler Scheme)
Nb_time_step = 1000;
Step = (t-0)/Nb_time_step;
T = 0+(0:Nb_time_step)'.*Step;
X = zeros(Nb_time_step+1, 1);
X(1) = x;
Y = zeros(Nb_time_step+1, 1);
Y(1) = x;
Integral_Approx = zeros(Nb_time_step+1, 1);
Integrale_Approx = zeros(Nb_time_step+1, 1);
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
       
   X(i+1) = X(i) + a*(b - X(i))*Step + sigma*gamrnd(c*Step,1/lambda); %gamma
   %Y(i+1) = Y(i) + a*(b - Y(i))*Step + sigma*randraw('ig',[c*Step,lambda]); % inverse gaussian
   Integral_Approx(i+1) = Integral_Approx(i) + X(i+1)*Step;
   %Integrale_Approx(i+1) = Integrale_Approx(i) + Y(i+1)*Step;
end

hold on;
plot(T, X);
%plot(T,Y,'r')
title('Siblation Levy OU piecewise-constant')
xlabel('temps')
legend('gamma','inverse gaussian')
plot(T, Integral_Approx, 'b');
%plot(T, Integrale_Approx, 'r');

%% Forble donnant la proba de survie comme dans l'article Jeanblanc-Crépey-Zargari (2009)

%Coeff_x = phi(T1, phi(T2-T1, phi(t-T2,0)));
%Coeff_b1 = xi(T1, phi(T2-T1, phi(t-T2, 0)));
%Coeff_b2 = xi(T2-T1, phi(t-T2,0));
%Coeff_b3 = xi(t-T2, 0);

%Proba_Survie_ver1 = exp(-(x*Coeff_x+b1*Coeff_b1+b2*Coeff_b2+b3*Coeff_b3))

%% Formule donnant la proba de survie dans le cas b = const par morceaux avec 3 paliers

Coeff_x = phii(t, a);
Coeff_b1 = xii(t,  a) -  xii(t-T1, a);
Coeff_b2 = xii(t-T1,  a) - xii(t-T2,  a);
Coeff_b3 = xii(t-T2,  a);
Coeff_c = zzeta(t, 0, a, sigma, lambda);
%Coeff_b7=0;
Proba_Survie_formule = exp(-(x*Coeff_x+b1*Coeff_b1+b2*Coeff_b2+b3*Coeff_b3 +c*Coeff_c))


%% Formule donnant la proba de survie dans le cas b = const
% Coeff_x = phii(t, 0, a)
% Coeff_b1 = xii(t, 0, a)
% Coeff_c = zzeta(t, 0, a, lambda)
% %Coeff_c_old = zzeta_old(t, a, lambda, 0);
% %Coeff_b7=0;
% Proba_Survie_ver2 = exp(-(x*Coeff_x + b1*Coeff_b1 + c*Coeff_c))
% %Proba_Survie_ver3 = exp(-(x*Coeff_x + b1*Coeff_b1 + c*Coeff_c_old))
% %Proba_Survie_ver4 = exp(-(x*Coeff_x + b1*Coeff_b1 ))



%% Test forble par Siblation de Monte Carlo pour 3 paliers
%  Nb_sibl = 1000;
%  Nb_Survie = 0;
%  for k=1:Nb_sibl
%  
%      Nb_time_step = 100;
%     Step = (t-0)/Nb_time_step;
%      T = 0+(0:Nb_time_step)'.*Step;
%      X = zeros(Nb_time_step+1, 1);
%      X(1) = x;
%      Integral_Approx = 0;
%      for i = 1:Nb_time_step
%         if (T(i)<T1) 
%             b = b1;
%         else
%            if (T(i)<T2)
%                  b = b2;
%             else
%                  b = b3;
%            end
%         end
%        X(i+1) = X(i) + a*(b - X(i))*Step + sigma*gamrnd(c*Step,1/lambda); %gamma
%        %X(i+1) = X(i) + a*(b - X(i))*Step + sigma*sqrt(Step)*sqrt(X(i))*normrnd(0,1); %gamma
%         Integral_Approx = Integral_Approx + X(i+1)*Step;
%      end
%      %Un default avant t ?
%      Survie_event = (Integral_Approx < exprnd(1));
%      Nb_Survie = Nb_Survie + Survie_event;
%      
%  end
%  Estimateur_MC = Nb_Survie/Nb_sibl 
%  Estimateur_ET = sqrt(Estimateur_MC*(1-Estimateur_MC));
%  Rayon_Int_Conf = 1.96*Estimateur_ET/sqrt(Nb_sibl);
%  Intervalle_Conf_95 = [Estimateur_MC - Rayon_Int_Conf  Estimateur_MC + Rayon_Int_Conf]

 %% Test forble par Siblation de Monte Carlo pour 1 palier (cas constant)
%  Nb_sibl = 50000;
%  Nb_Survie = 0;
%  for k=1:Nb_sibl
%  
%      Nb_time_step = 100;
%     Step = (t-0)/Nb_time_step;
%      T = 0+(0:Nb_time_step)'.*Step;
%      X = zeros(Nb_time_step+1, 1);
%      X(1) = x;
%      Integral_Approx = 0;
%      for i = 1:Nb_time_step
%        X(i+1) = X(i) + a*(b1 - X(i))*Step + gamrnd(c*Step, 1/lambda); %gamma
%        %X(i+1) = X(i) + eta*(b - X(i))*Step + nu*Step*sqrt(X(i))*normrnd(0,1); %gamma
%         Integral_Approx = Integral_Approx + X(i+1)*Step;
%      end
%      %Un default avant t ?
%      Survie_event = (Integral_Approx < exprnd(1));
%      Nb_Survie = Nb_Survie + Survie_event;
%      
%  end
%  Estimateur_MC = Nb_Survie/Nb_sibl 
%  Estimateur_ET = sqrt(Estimateur_MC*(1-Estimateur_MC));
%  Rayon_Int_Conf = 1.96*Estimateur_ET/sqrt(Nb_sibl);
%  Intervalle_Conf_95 = [Estimateur_MC - Rayon_Int_Conf  Estimateur_MC + Rayon_Int_Conf]
 
  
 
