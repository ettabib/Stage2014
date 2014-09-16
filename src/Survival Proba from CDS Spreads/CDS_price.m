function res = CDS_price(pricing_date, next_coupon_date, pillars,  intensities, Recovery, Interest, Intensity_model, Param)
% --------------------------------------------------------------------------------------------------
% Compute Default Leg, Unitary Premium Leg and Fair Spread for a CDS with
% specified caracteristics
% ------------------------- INPUT -------------------------------------------------------------------------
% pricing_date                 
% next_coupon_date              ... coupon date just after the pricing date    
% pillars                       ... [T1; ...; Tp ] : row vector of increasing dates or pillars, Tp must be maturity date
% Intensities                   ... [lambda_1; ...; lambda_p ] : piecewise constant mean-reversion coefficients, lambda_1 associated with period [0,T1], ...,lambda_p with [Tp-1,Tp]
% Recovery                      ... Recovery rate
% Interest                      ... Interest rate
% Intensity_model               ... Identifier for the intensity model : 1 for CIR, 2 for Gamma-OU, 3 for IG-OU, 4 for CGMI-OU  
% Param                         ... Parameters [X0; Param]
% ------------------------- OUTPUT -------------------------------------------------------------------------
% res.DL                        ... Default leg
% res.PL                        ... Premium Leg for a coupon spread of 1 monetary unit
% res.Fair_spread               ... Fair Spread
% -------------------------------------------------------------------------
% -------------------------
% sample call 1: 
%   toDate = '17-Dec-2007';                                                             % today
%   coupDate = '20-Dec-2007';                                                           % first coupon payment date
%   Pillars = {'10-Dec-2010'; '12-Dec-2012'; '14-Dec-2014'; '14-Dec-2017'};             % pillars
%   Intensities = 0.1*ones(size(Pillars,1),1);   
% Intensities
%   Recovery = 0.40;
%   Interest = 0.03;
% x = 0.1; %starting point of the CIR process
% eta = 3;
% nu = 0.1;
% %  eta = 3;
% %  nu = 0.1;
% %  x = 3;
% Intensity_model 
%   Param = [x; eta; nu]  
%   Intensity_model = 1;
%   res = CDS_price(toDate, coupDate, Pillars, Intensities, Recovery, Interest, Intensity_model, Param)
% -------------------------------------------------------------------------



% Check of the number of pillars corresponds to the number of 'steps' for
% the piecewise constant intensities
Nb_pillars = size(pillars,1);
Nb_intensities = size(intensities,1);

if (Nb_pillars ~= Nb_intensities) 
    exception = MException('VerifyInput:OutOfBounds', ...
       'Number of pillars and intensities must be the same');
    throw(exception);
end

%Convert cell vector of dates pillars into vectors of real number t_pillars
%corresponding to lengths between today and pillars dates 
t_pillars = zeros(Nb_pillars, 1);
for i=1:Nb_pillars
    t_pillars(i) = (datenum(pillars(i),'dd-mmm-yyyy')-datenum(pricing_date,'dd-mmm-yyyy'))/365;
end


%Create time vector t_scheduled used for computation of DL, PL and Fair
%Spreads. In this computation, it corresponds to premium payment dates
dz = 0.25;  %length in years bewteen two premium payment dates
matDate = pillars(Nb_pillars,1);
coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(next_coupon_date,'dd-mmm-yyyy'))/365;
dz1 = (datenum(next_coupon_date,'dd-mmm-yyyy')-datenum(pricing_date,'dd-mmm-yyyy'))/365;
T = dz1 + dz*round(coupon2mat/dz);
t_scheduled = (dz1:dz:T)';
t_scheduled(end) = t_pillars(end);
%t_scheduled
numPeriod = size(t_scheduled,1);
%numPeriod = round(coupon2mat/dz) + 1;  

%Extend the vector of intensities (size = NbPillars) on a grid with time
%interval defined by length between two premium payment dates (size =
%numPeriod)
% intensities2 = zeros(numPeriod,1);
% Periods = zeros(numPeriod,1);  %number of period of each time t in t_scheduled
% intensities2(t_scheduled <= t_pillars(1))=intensities(1); %get intensities corresponding to times in t_scheduled smaller than t_pillars(1)
% Periods(t_scheduled <= t_pillars(1)) = 1;
% for i=2:Nb_pillars
%     %ind = find(t_scheduled>pillars(i-1) & t_scheduled <= pillars(i));
%     intensities2(t_scheduled > t_pillars(i-1) & t_scheduled <= t_pillars(i))=intensities(i);
%     Periods(t_scheduled > t_pillars(i-1) & t_scheduled <= t_pillars(i))=i;
% end

switch Intensity_model
    case 1 %CIR case
        Proba_survi = Proba_survi_CIR(t_scheduled, intensities, t_pillars, Param);
    case 2
        
        Proba_survi = Proba_survi_OU_Gamma(t_scheduled, intensities, t_pillars, Param);        
        %TEST GAMMA COMPENSATED
        %Compensator = Param(3)*Param(4)/(Param(2)*Param(5));
        %Proba_survi = Proba_survi_OU_Gamma(t_scheduled, intensities - Compensator, t_pillars, Param);
        
    case 3
        Proba_survi = Proba_survi_OU_IG(t_scheduled, intensities, t_pillars, Param);
    case 4
        Proba_survi = Proba_survi_OU_CGMY(t_scheduled, intensities, t_pillars, Param);
    otherwise
        exception = MException('VerifyInput:OutOfBounds', ...
       'Unknown model reference');
    throw(exception);     
end


price = CDS_price_from_survi(t_scheduled, Proba_survi, Recovery, Interest);

res.DL = price.DL;
res.PL = price.PL;
res.Fair_spread = price.Fair_spread;



function res = CDS_price_from_survi(t_scheduled, proba_survi, Recovery, Interest)
% --------------------------------------------------------------------------------------------------
% Compute Default Leg, Unitary Premium Leg and Fair Spread for a particular schedule of premium
% payment dates and related survival probabilities
% ------------------------- INPUT -------------------------------------------------------------------------
% t_scheduled                   ... [t1; ...; tn] : row vector of increasing premium payment dates
% proba_survi                   ... Survival probabilities at [t1; ...; tn]
% Recovery                      ... Recovery rate
% Interest                      ... Interest rate
% ------------------------- OUTPUT -------------------------------------------------------------------------
% res.DL                        ... Default leg
% res.PL                        ... Premium Leg for a coupon spread of 1 monetary unit
% res.Fair_spread                        ... Fair Spread
% -------------------------------------------------------------------------
% -------------------------
% sample call: 
%   t_scheduled = (0.25:0.25:5)';
%   lambda = 0.1;
%   proba_survi = exp(-lambda*t_scheduled);
%   Recovery = 0.40;
%   Interest = 0.03;
%   res = CDS_price_from_survi(t_scheduled, proba_survi, Recovery, Interest)
% -------------------------------------------------------------------------
% -------------------------


Nb_period = size(t_scheduled,1);
Nb_proba  = size(proba_survi,1);


if (Nb_period ~= Nb_proba) 
    exception = MException('VerifyInput:OutOfBounds', ...
       'Number of periods and probabilities must be the same');
    throw(exception);
end

Proba_default = 1 - proba_survi;
Discount = exp(-Interest*t_scheduled);
Discounted_proba = Discount.*Proba_default;
Period_length = diff([0;t_scheduled]);

% Default leg
%DL = (1-Recovery)*(Discounted_proba(Nb_period) + Interest*sum(Discounted_proba.*Period_length));

Discounted_survi = Discount.*proba_survi;
DL = (1-Recovery)*(1 - Discounted_survi(Nb_period) - Interest*sum(Discounted_survi.*Period_length));

% Premium Leg
PL = sum(Discount.*Period_length.*proba_survi);

res.DL = DL;
res.PL = PL;
res.Fair_spread = DL/PL;


% function res  = proba_survi(t_scheduled, CIR_Param,  intensities, t_pillars)
% % --------------------------------------------------------------------------------------------------
% % Compute survival probability in a reduced-from model where intensities
% % driven by a CIR process with piecewise constant mean-reversion
% % coefficient
% % --------------------------------------------------------------------------------------------------
% % t_scheduled                   ... [t1; ...; tn] : column vector of increasing dates
% % CIR_Param                     ... Parameters [X0;a;b] X0 initial value, a speed of mean-reversion, b volaility coefficient
% % intensities                   ... row column of piecewise-constant mean-reversion parameters 
% % t_pillars                     ... row column of time-to-pillars
% % --------------------------------------------------------------------------------------------------
% % sample call: 
% %%   t_scheduled = (0.25:0.25:5)'
% % t_scheduled = [2.5; 4; 5]
% % x = 0.1; %starting point of the CIR process
% % eta = 3;
% % nu = 0.1;
% % %  eta = 3;
% % %  nu = 0.1;
% % %  x = 3;
% %   CIR_Param = [x; eta; nu]  
% %   t_pillars = [2;3;5]
% %   intensities = [0.5; 1; 0.1]
% %   res  = proba_survi_test(t_scheduled, CIR_Param,  intensities, t_pillars)
% % -------------------------------------------------------------------------
% % -------------------------
% 
% x = CIR_Param(1);
% a = CIR_Param(2);
% sigma = CIR_Param(3);
% c = CIR_Param(4);
% lambda = CIR_Param(5);
% 
% Nb_period = size(t_scheduled,1);
% Nb_lambda  = size(intensities,1);
% Nb_pillars = size(t_pillars,1);
% 
% Proba_Survie = zeros(Nb_period,1);
% 
% if (Nb_pillars ~= Nb_lambda) 
%     exception = MException('VerifyInput:OutOfBounds', ...
%        'Number of components in t_pillars and intensities must be the same');
%     throw(exception);
% end
% 
% for k = 1:Nb_period
%     t = t_scheduled(k); %Current time t
%     Pillars_inf_t = [0; t_pillars(t>t_pillars)]; %column vector containing (in increasing order) all pillar dates (strictly) smaller than t
%     Nb_Pillars_t = size(Pillars_inf_t,1); %Nb of Pillars smaller than t (including T=0)
%     Intensities_t = intensities(1:Nb_Pillars_t,1); %All relevant piecewise-constant param to compute survival proba at time t
%     Diff_t = t - Pillars_inf_t;
%     Xi_t = xi_OU(Diff_t, a); %does not depend on intensities and starting point
%     Diff_Xi_t = - diff([Xi_t ; 0]); %does not depend on intensities and starting point
%     %Proba_Survie(k) = exp(- phi(t,param)*x - sum(Diff_Xi_t.*Intensities_t));
%     
%     Coeff_c = zzeta(t, 0, a, sigma, lambda);
%     
%     Proba_Survie(k) = exp(- phi_OU(t,a)*x - Diff_Xi_t'*Intensities_t - c*Coeff_c);
% end
% 
% res = Proba_Survie;

%Proba_Survie_formule = exp(-(x*Coeff_x+b1*Coeff_b1+b2*Coeff_b2+b3*Coeff_b3 +c*Coeff_c))

% function res  = proba_survi(t_scheduled, intensities)
% % --------------------------------------------------------------------------------------------------
% % Compute survival probability in a reduced-from model where intensities
% % are piecewise constant
% % --------------------------------------------------------------------------------------------------
% % t_scheduled                   ... [t1; ...; tn] : column vector of increasing dates
% % intensities                   ... row column of intensities corresponding to [t1; ...; tn]
% % --------------------------------------------------------------------------------------------------
% % sample call: 
% %   t_scheduled = (0.25:0.25:5)';
% %   intensities = 10*rand(size(t_scheduled,1),1);
% %   res = proba_survi(t_scheduled, intensities)
% % -------------------------------------------------------------------------
% % -------------------------
% 
% 
% Nb_period = size(t_scheduled,1);
% Nb_lambda  = size(intensities,1);
% 
% if (Nb_period ~= Nb_lambda) 
%     exception = MException('VerifyInput:OutOfBounds', ...
%        'Number of components in t_scedulded and intensities must be the same');
%     throw(exception);
% end
% 
% Period_length = diff([0;t_scheduled]);
% Hazard_rate = cumsum(intensities.*Period_length);
% res = exp(-Hazard_rate);


