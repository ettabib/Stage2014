%% Plot Survival distr function for a Gamma-Levy OU intensity model
% mean-reversion parameter is a piecewise-constant function of time


bp = 0.0001;

%Lambda is fixed at 1/Taille_Saut
Taille_Saut = 1*bp;
lambda = 1/Taille_Saut;

% definition grille a (vitesse retour moyenne) et grille c
nb_a = 2;
nb_s = 4; %nb sigma param

min_a = 0.01;
max_a = 2;
pas_a = (max_a-min_a)/nb_a;

min_s = 0.2;
max_s = 2;
pas_s = (max_s-min_s)/nb_s;

a = min_a:pas_a:max_a;
s = min_s:pas_s:max_s;

% a = [1];
% s = [1];

Size_a = size(a,2);
Size_s = size(s,2);

%% Assumption on recovery rate and interest rate
Recovery = 0.40;
Interest = 0.03;
%Interest = 0;

Intensity_model = 1; %CIR model

%%  Bootstrap procedure
load('CDX9_17_Dec_2007.mat'); %Load Index_data containing all relevant information regarding this date
                              % and in particular CDS spreads of all consituents at different pillars

toDate = Index_data.Today;                                   % today
coupDate = Index_data.coupDate;                             % first coupon payment date
Pillars = Index_data.CDS_Pillars;

Index_AIG = 5;
CDS_spreads_AIG = bp*Index_data.CDS_Spreads(Index_AIG,:);

nb_x0 = 10;
First_Int_Approx = CDS_spreads_AIG(1)/(1-Recovery);
min_x0 = 0.01*First_Int_Approx;
max_x0 = 2.5*First_Int_Approx;
pas_x0 = (max_x0-min_x0)/nb_x0;
x0 = min_x0:pas_x0:max_x0;

%x0 = First_Int_Approx;
Size_x0 = size(x0,2);

Nb_pillars = size(Pillars, 1);

%Convert cell vector of dates pillars into vectors of real number t_pillars
%corresponding to lengths between today and pillars dates 
t_pillars = zeros(Nb_pillars, 1);
for i=1:Nb_pillars
    t_pillars(i) = (datenum(Pillars(i),'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
end


CDS_market_spreads = CDS_spreads_AIG';

Res = cell(Size_a, Size_s, Size_x0);

for i = 1:Size_a
    for j = 1:Size_s
        for k = 1:Size_x0
            %a_c = [a(i); c(j)];
            CIR_Param = [x0(k); a(i); s(j)]; %Param with starting point x0
                                             %Bootstraped_Int = Bootstrap_intensities(toDate, coupDate, Pillars, CDS_market_spreads, OU_Levy_Param, Recovery, Interest);

            Bootstraped_Int = Bootstrap_intensities_ver_1(toDate, coupDate, Pillars, CDS_market_spreads, Recovery, Interest, Intensity_model, CIR_Param);

            b = Bootstraped_Int;
            Res{i,j,k} = b;
        end        
    end
end

%% Compute Survival Proba
pas_t = 0.25;
t_scheduled=(0:pas_t:4)';
size_t_sched = size(t_scheduled);

P = ones(Size_x0,size(t_scheduled,1));
for i = 1:Size_a
    for j = 1:Size_s
        for k = 1:Size_x0
            b = Res{i,j,k};
            CIR_Param = [x0(k); a(i); s(j)];
            P(k,1:size(t_scheduled,1)) = Proba_survi_CIR(t_scheduled, b, t_pillars, CIR_Param);
           % plot(t_scheduled, P(k,1:size(t_scheduled,1)));           
        end
    end
end

Function_Prob =1 - P;

plot(t_scheduled,Function_Prob);


%% Generate the sum of all payoff terms between t and T, all terms discounted back at t, not subject to the counterparty default risk
%Distribution function of the underlying
Distribution_function_Tau1 = 0:0.001:1;
% Distribution function of the conterparty
DateQuotation = 1;
Distribution_function_Tau2 = Function_Prob(DateQuotation,:);
Spreads = 
T = max(t_scheduled);
delta = 1/4;
t = 1/4 * T;
Payoff_Without_CR = 0;
bt = t_scheduled(find(t_scheduled>t));
bt = find(t_scheduled==bt(1));

S = 0;
i = bt;
while (t_scheduled(i) < T)
   S = S + Discount(t,t_scheduled(bt),R) *  
    
end

