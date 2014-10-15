%% Plot Survival distr function for a Gamma-Levy OU intensity model
% mean-reversion parameter is a piecewise-constant function of time
% on a introduit une possibilite de rang de parametre afin de tester
% l'infulence de ces derniers sur l'evolution des courbes

bp = 0.0001;

%Lambda is fixed at 1/Taille_Saut
Taille_Saut = 1*bp;
lambda = 1/Taille_Saut;

% definition grille a (vitesse retour moyenne) et grille c

% a = (1);
% s = (1);

nb_a = 0;
nb_s = 0; %nb sigma param
nb_x0 = 10;
n_point = 10; %%%%% nombre de points : nb dates de quotation echantillonees


min_a = 1;
max_a = 2;
pas_a = (max_a-min_a)/nb_a;

min_s = 0.2;
max_s = 5;
pas_s = (max_s-min_s)/nb_s;

a = min_a:pas_a:max_a;
s = min_s:pas_s:max_s;


Size_a = size(a,2);
Size_s = size(s,2);

% Definition du nombre de points
debut = 1000; % premiere date de quotation
fin = 1338; % derniere date de quotation
delta = floor((fin - debut) / n_point); % le pas entre deux points


%% Assumption on recovery rate and interest rate
Recovery = 0.40;
Interest = 0.03;
Intensity_model = 1; %CIR model

%% Import the data
[~, ~, raw] = xlsread('/Users/mohammad/Documents/MyProject/Stage2014/data/France_1.xlsx','Sheet1');
raw = raw(1:end,:);


%%  Bootstrap procedure
%load('CDX9_17_Dec_2007.mat'); %Load Index_data containing all relevant information regarding this date
% and in particular CDS spreads of all consituents at different pillars

toDate = {'16-08-2005'};                          % today
coupDate = {'16-08-2005'};                        % first coupon payment date
Pillars = {'16-08-2006' ,  '16-08-2007' ,'16-08-2008', '16-08-2009', '16-08-2010' ,'16-08-2011', '16-08-2012', '16-08-2013' ,'16-08-2014' ,'16-08-2015'};
% Pillars = datestr(Pillars);
% Index_AIG = 5;
nb_cotation = size(raw(:,1),1) - 1;
CDS_spreads_AIG = bp*cell2mat(raw(1:nb_cotation,2:11));

% Ploting according to the diversity of X0

Nb_pillars = size(Pillars, 2);

% Convert cell vector of dates pillars into vectors of real number t_pillars
% corresponding to lengths between today and pillars dates
% t_pillars : the vector of maturities

t_pillars = zeros(Nb_pillars, 1);
for i=1:Nb_pillars
    t_pillars(i) = (datenum(Pillars(i),'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
end

% CDS_market_spreads corresponds to the matrix of (quotation,spreads)
% (1338 x 11)
CDS_market_spreads = CDS_spreads_AIG;

% Defining The parameter X0

First_Int_Approx = CDS_spreads_AIG(1)/(1-Recovery);
min_x0 = 0.01*First_Int_Approx;
max_x0 = 2.5*First_Int_Approx;
pas_x0 = (max_x0-min_x0)/nb_x0;
x0 = min_x0:pas_x0:max_x0;

%x0 = First_Int_Approx;
Size_x0 = size(x0,2);


L = 0;
%resultat du boot strap : tableau a x s x X0 x quotation
Res = cell(Size_a, Size_s, Size_x0,n_point);

%boucle du bootstrap
for i = 1:Size_a
    for j = 1:Size_s
        for k = 1:Size_x0
            for m = 1:n_point
                %a_c = [a(i); c(j)];
                CIR_Param = [x0(k); a(i); s(j)]; %Param with starting point x0                                                 
                L(m) = debut + delta * m; %date de quotation
                Res{i,j,k,m} = Bootstrap_intensities_ver_1(toDate, coupDate, Pillars', CDS_market_spreads(L(m),:)', Recovery, Interest, Intensity_model, CIR_Param);                                                 
            end
        end        
    end
end



%% Compute Survival Proba
% in this section we compute the distribution function
pas_t = 0.25;
t_scheduled=(0:pas_t:Nb_pillars)';
size_t_sched = size(t_scheduled);
size_t_sched = size_t_sched(1);

P = ones(n_point,size(t_scheduled,1));
Survival_P = cell(Size_a,Size_s,Size_x0,n_point,size_t_sched);
figure
hold all;
for i = 1:Size_a %parcours les a     
    for j = 1:Size_s %parcours les s
        for k = 1:Size_x0 % parcours les x0
            for m = 1:n_point % parcours le nombre de courbes des dates de quotation
                b = Res{i,j,k,m};
                CIR_Param = [x0(k); a(i); s(j)];
                Survival_P{i,j,k,m} = Proba_survi_CIR(t_scheduled, b, t_pillars, CIR_Param);                
                plot(t_scheduled, Survival_P{i,j,k,m});
            end
        end
    end
end



%% %% Generate the sum of all payoff terms between t and T, all terms discounted back at t, not subject to the counterparty default risk
Function_Prob = 1 - Survival_P{1,1,1,1};
%plot(t_scheduled,Function_Prob);

%Distribution function of the underlying
Distribution_function_Tau1 = 0:0.001:1;
% Distribution function of the conterparty
DateQuotation = 1;
Distribution_function_Tau2 = Function_Prob;
T = max(t_scheduled);
delta = 1/4;
t = 1/4 * T;
Payoff_Without_CR = 0;

% bt : the first coupon date over t
bt = t_scheduled(t_scheduled>t);
bt = find(t_scheduled==bt(1));

% computing sum S = si * D(t,t_i) * dt , for t_i > t
S = 0;
i = bt;
while (t_scheduled(i) < T)
    pi = find(t_pillars==max(t_pillars(t_pillars <= t_scheduled(i))));
    S = S + Discount(t,t_scheduled(i),Recovery) * CDS_market_spreads(DateQuotation,pi) * pas_t;
    i = i + 1 ;
end
