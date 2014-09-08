function res  = Proba_survi_OU_CGMY(t_scheduled,  intensities, t_pillars, CGMY_OU_Param)
% --------------------------------------------------------------------------------------------------
% Compute survival probability in a reduced-from model where intensities
% driven by a IG-OU process with piecewise constant mean-reversion
% coefficient
% --------------------------------------------------------------------------------------------------
% t_scheduled                   ... [t1; ...; tn] : column vector of increasing dates
% CIR_Param                     ... Parameters [X0;a;b] X0 initial value, a speed of mean-reversion, b volaility coefficient
% intensities                   ... row column of piecewise-constant mean-reversion parameters 
% t_pillars                     ... row column of time-to-pillars
% --------------------------------------------------------------------------------------------------
% sample call: 
%%   t_scheduled = (0.25:0.25:5)'
% t_scheduled = [2.5; 4; 5]
% bp = 0.0001;
% x = 0.008; %starting point of the CIR process
% a = 1; 
% sigma = 1;
% c = 20; %Nb moyen de saut par an
% Taille_Saut = 3*bp;
% lambda = 1/Taille_Saut;
% Gamma_OU_Param = [x; a; sigma; c; lambda]  
% intensities = [0.5; 1; 0.1]  %les b_i
% t_pillars = [2;3;5]
% res  = Proba_survi_OU_IG(t_scheduled,intensities, t_pillars, Gamma_OU_Param)
% -------------------------------------------------------------------------
% -------------------------

x = CGMY_OU_Param(1);
a = CGMY_OU_Param(2);
sigma = CGMY_OU_Param(3);
c = CGMY_OU_Param(4);
M = CGMY_OU_Param(5);
Y = CGMY_OU_Param(6);

Nb_period = size(t_scheduled,1);
Nb_lambda  = size(intensities,1);
Nb_pillars = size(t_pillars,1);

Proba_Survie = zeros(Nb_period,1);

if (Nb_pillars ~= Nb_lambda) 
    exception = MException('VerifyInput:OutOfBounds', ...
       'Number of components in t_pillars and intensities must be the same');
    throw(exception);
end

for k = 1:Nb_period
    t = t_scheduled(k); %Current time t
    Pillars_inf_t = [0; t_pillars(t>t_pillars)]; %column vector containing (in increasing order) all pillar dates (strictly) smaller than t
    Nb_Pillars_t = size(Pillars_inf_t,1); %Nb of Pillars smaller than t (including T=0)
    Intensities_t = intensities(1:Nb_Pillars_t,1); %All relevant piecewise-constant param to compute survival proba at time t
    Diff_t = t - Pillars_inf_t;
    Xi_t = xi_OU(Diff_t, a); %does not depend on intensities and starting point
    Diff_Xi_t = - diff([Xi_t ; 0]); %does not depend on intensities and starting point
    %Proba_Survie(k) = exp(- phi(t,param)*x - sum(Diff_Xi_t.*Intensities_t));
    
    Coeff_c = zeta_CGMY(t, 0, a, sigma,M,Y);
    
    Proba_Survie(k) = exp(- phi_OU(t,a)*x - Diff_Xi_t'*Intensities_t - c*Coeff_c);
end

res = Proba_Survie;