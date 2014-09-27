%% Plot Survival distr function for a Gamma-Levy OU intensity model
% mean-reversion parameter is a piecewise-constant function of time
% on a introduit une possibilite de rang de parametre afin de tester
% l'infulence de ces derniers sur l'evolution des courbes

bp = 0.0001;

%Lambda is fixed at 1/Taille_Saut
Taille_Saut = 1*bp;
lambda = 1/Taille_Saut;

% definition grille a (vitesse retour moyenne) et grille c

a = (1);
s = (1);

Size_a = size(a,2);
Size_s = size(s,2);


% nb_a = 2;
% nb_s = 10; %nb sigma param
% 
% min_a = 0.01;
% max_a = 2;
% pas_a = (max_a-min_a)/nb_a;
% 
% min_s = 0.2;
% max_s = 2;
% pas_s = (max_s-min_s)/nb_s;

% a = min_a:pas_a:max_a;
% s = min_s:pas_s:max_s;




%% Assumption on recovery rate and interest rate
Recovery = 0.40;
Interest = 0.03;
%Interest = 0;

Intensity_model = 1; %CIR model

%% Import the data
[~, ~, raw] = xlsread('/Users/mohammad/Documents/MyProject/Stage2014/data/France_1.xlsx','Sheet1');
raw = raw(1:end,:);


%%  Bootstrap procedure
%load('CDX9_17_Dec_2007.mat'); %Load Index_data containing all relevant information regarding this date
% and in particular CDS spreads of all consituents at different pillars

toDate = {'16-08-2005'};                                   % today
coupDate = {'16-08-2005'};                             % first coupon payment date
Pillars = {'16-08-2006' ,  '16-08-2007' ,'16-08-2008', '16-08-2009', '16-08-2010' ,'16-08-2011', '16-08-2012', '16-08-2013' ,'16-08-2014' ,'16-08-2015'};
% Pillars = datestr(Pillars);
% Index_AIG = 5;
nb_cotation = size(raw(:,1),1) - 1;
CDS_spreads_AIG = bp*cell2mat(raw(1:nb_cotation,2:11));

% Ploting according to the diversity of X0

Nb_pillars = size(Pillars, 2);

%Convert cell vector of dates pillars into vectors of real number t_pillars
%corresponding to lengths between today and pillars dates
% t_pillars : the vector of maturities
t_pillars = zeros(Nb_pillars, 1);
for i=1:Nb_pillars
    t_pillars(i) = (datenum(Pillars(i),'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
end

% CDS_market_spreads corresponds to the matrix of (quotation,spreads)
% (1338 x 11)
CDS_market_spreads = CDS_spreads_AIG;


% Definition du nombre de points
debut = 1000; %debut
fin = 1100; %fin
n_point = 11; % nombre de points : dates de quotation
delta = floor((fin - debut) / n_point); % le pas entre deux points

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
                L(m) = a + delta * m; %date de quotation
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

P = ones(n_point,size(t_scheduled,1));
for i = 1:Size_a %parcours les a     
    for j = 1:Size_s %parcours les s
        for k = 1:Size_x0 % parcours les x0
            for m = 1:n_point % parcours le nombre de courbes des dates de quotation
                b = Res{i,j,k,m};
                CIR_Param = [x0(k); a(i); s(j)];
                P(m,:) = Proba_survi_CIR(t_scheduled, b, t_pillars, CIR_Param);
                plot(t_scheduled, P(m,:));
            end
        end
    end
end

Function_Prob =1 - P;

plot(t_scheduled,Function_Prob);