function res  = Proba_survi_CIR(t_scheduled,  intensities, t_pillars, CIR_Param)
% --------------------------------------------------------------------------------------------------
% Compute survival probability in a reduced-from model where intensities
% driven by a CIR process with piecewise constant mean-reversion
% coefficient
% --------------------------------------------------------------------------------------------------
% t_scheduled                   ... [t1; ...; tn] : column vector of increasing dates
% CIR_Param                     ... Parameters [X0;a;b] X0 initial value, a speed of mean-reversion, b volaility coefficient
% intensities                   ... row column of intensities corresponding to the intensities
% t_pillars                     ... row column of time-to-pillars
% --------------------------------------------------------------------------------------------------
% sample call: 
%%   t_scheduled = (0.25:0.25:5)'
% t_scheduled = [2.5; 4; 5]
% x = 0.1; %starting point of the CIR process
% eta = 3;
% nu = 0.1;
% %  eta = 3;
% %  nu = 0.1;
% %  x = 3;
%   CIR_Param = [x; eta; nu]  
%   t_pillars = [2;3;5]
%   intensities = [0.5; 1; 0.1]
%   res  = Proba_survi_CIR(t_scheduled,  intensities, t_pillars, CIR_Param)
% -------------------------------------------------------------------------
% -------------------------

x = CIR_Param(1);
a = CIR_Param(2);
b = CIR_Param(3);
param = [a; b];

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
    t = t_scheduled(k);
    Pillars_inf_t = [0; t_pillars(t>t_pillars)];
    Nb_Pillars_t = size(Pillars_inf_t,1);
    Intensities_t = intensities(1:Nb_Pillars_t,1);
    Diff_t = t - Pillars_inf_t;
    Xi_t = xi_CIR(Diff_t, param);
    Diff_Xi_t = - diff([Xi_t ; 0]);
    Proba_Survie(k) = exp(- phi_CIR(t,param)*x - sum(Diff_Xi_t.*Intensities_t));
end

res = Proba_Survie;

