function res = Bootstrap_intensities_ver_1(pricing_date, next_coupon_date, pillars, CDS_market_spreads, Recovery, Interest, Intensity_model, Param)
% --------------------------------------------------------------------------------------------------
% Find piecewise constant intensities in a CDS pricer that match CDS Spreads at different pillars
% ------------------------- INPUT -------------------------------------------------------------------------
% pricing_date                 
% next_coupon_date              ... coupon date just after the pricing date    
% pillars                       ... [T1; ...; Tp ] : row vector of increasing dates or pillars, Tp must be maturity date
% CIR_Param                     ... Parameters [a;b] a speed of mean-reversion, b volatility coefficient
% CDS_market_spreads            ... [S_1; ...; S_p ] : S_1 associated with protection period [0,T1], ..., S_p with [0,Tp]
% Recovery                      ... Recovery rate
% Interest                      ... Interest rate
% ------------------------- OUTPUT -------------------------------------------------------------------------
% intensities               ... [lambda_1; ...; lambda_p ] : piecewise constant intensities, lambda_1 associated with period [0,T1], ...,lambda_p with [Tp-1,Tp]
% -------------------------------------------------------------------------
% sample call 1: 
%   toDate = '17-Dec-2007';                                                             % today
%   coupDate = '20-Dec-2007';                                                           % first coupon payment date
%   Pillars = {'10-Dec-2010'; '12-Dec-2012'; '14-Dec-2014'; '14-Dec-2017'};             % pillars
%   CDS_market_spreads = [29; 47; 54; 63];   CDS_market_spreads =CDS_market_spreads/10000;
%   Recovery = 0.40;
%   Interest = 0.03;
%   res = Bootstrap_intensities(toDate, coupDate, Pillars, CDS_market_spreads, Recovery, Interest)
% -------------------------------------------------------------------------
% sample call 2: 
%   toDate = '24-Jul-2009';                                                             % today
%   coupDate = '20-Sep-2009';                                                           % first coupon payment date
%   Pillars = {'24-Jul-2011'; '24-Jul-2013'};                               % pillars 1 an et 2ans
%   CDS_market_spreads = [29; 54];   CDS_market_spreads =CDS_market_spreads/10000;
%   Recovery = 0.40;
%   Interest = 0.03;
%   res = Bootstrap_intensities(toDate, coupDate, Pillars, CDS_market_spreads, Recovery, Interest)
% -------------------------------------------------------------------------

Nb_pillars = size(pillars,1);
Nb_spreads = size(CDS_market_spreads,1);


if (Nb_pillars ~= Nb_spreads) 
    exception = MException('VerifyInput:OutOfBounds', ...
       'Number of pillars and market spreads must be the same');
    throw(exception);
end

intensities = [];
CDS_market_spread = CDS_market_spreads(1);
Int_pillars = pillars(1,1);
last_int = fzero(@(last_int) Diff_Market_Model(CDS_market_spread, Int_pillars, intensities,  last_int, pricing_date, next_coupon_date, Recovery, Interest, Intensity_model, Param),  0.01);
intensities = [intensities; last_int];


for i=2:Nb_pillars
    CDS_market_spread = CDS_market_spreads(i);
    Int_pillars = pillars(1:i,1);
    last_int = fzero(@(last_int) Diff_Market_Model(CDS_market_spread, Int_pillars, intensities,  last_int, pricing_date, next_coupon_date, Recovery, Interest, Intensity_model, Param),  0.01);
    intensities = [intensities; last_int];
end

%-----------------------------DEBUG----------------------------------------
% Tolerance = 1e-5;
% Model_spreads = zeros(Nb_pillars, 1);
% for i=1:Nb_pillars
%     price = CDS_price(pricing_date, next_coupon_date, pillars(1:i,1), intensities(1:i,1), Recovery, Interest);
%     Model_spreads(i) = price.Fair_spread;
% end
% Diff_Spreads = abs(CDS_market_spreads - Model_spreads)
% if sum(Diff_Spreads.^2)<Tolerance
%     disp('OK');
% end

res = intensities;

%  = fzero(fun,x0)
% 
%     options = optimset('MaxFunEvals',500,'MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'Display','final');    
%     para = lsqnonlin(@(para) findSprd(attpts,mktsprd,N0,para,1),para0,zeros(size(para0)),[],options);  




function res = Diff_Market_Model(CDS_market_spread, pillars, intensities_without_one, last_Intensity, pricing_date, next_coupon_date, Recovery, Interest, Intensity_model, Param)
% --------------------------------------------------------------------------------------------------
% Compute the difference between (a) market spread and (b) model spread in a reduced form model with piecewise constant
% intensities such that:
% ------------------------- INPUT ------------------------------------------------------------------------- 
% CDS_market_spread
% pillars                       ... [T1; ...; Tk ] : row vector of increasing dates or pillars, Tp must be maturity date
% Intensities_without_one       ... [lambda_1; ...; lambda_{k-1} ] : piecewise constant intensities, lambda_1 associated with period [0,T1], ...
% Last_Intensity                ... lambda_{k-1}
% CIR_Param                     ... Parameters [X0;a;b] X0 initial value, a speed of mean-reversion, b volatility coefficient
% pricing_date                 
% next_coupon_date              ... coupon date just after the pricing date
% Recovery                      ... Recovery rate
% Interest                      ... Interest rate
% ------------------------- OUTPUT -------------------------------------------------------------------------
% Fair Spread
% -------------------------------------------------------------------------
intensities = [intensities_without_one; last_Intensity];
Nb_pillars = size(pillars,1);
Nb_intensities = size(intensities,1);

if (Nb_pillars ~= Nb_intensities) 
    exception = MException('VerifyInput:OutOfBounds', ...
       'Number of pillars and intensities must be the same');
    throw(exception);
end

price = CDS_price(pricing_date, next_coupon_date, pillars, intensities, Recovery, Interest, Intensity_model, Param);
CDS_model_spread = price.Fair_spread;

res = CDS_market_spread - CDS_model_spread;