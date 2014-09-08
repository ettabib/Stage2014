%% Plot Survival distr function for a Gamma-Levy OU intensity model
% mean-reversion parameter is a piecewise-constant function of time


bp = 0.0001;

%Lambda is fixed at 1/Taille_Saut
Taille_Saut = 1*bp;
lambda = 1/Taille_Saut;

% definition grille a (vitesse retour moyenne) et grille c
nb_a = 2;
nb_s = 10; %nb sigma param

min_a = 0.01;
max_a = 2;
pas_a = (max_a-min_a)/nb_a;

min_s = 0.2;
max_s = 2;
pas_s = (max_s-min_s)/nb_s;

a = min_a:pas_a:max_a
s = min_s:pas_s:max_s

a = [1];
s = [1];

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
CDS_spreads_AIG = bp*Index_data.CDS_Spreads(Index_AIG,:)

nb_x0 = 10;
First_Int_Approx = CDS_spreads_AIG(1)/(1-Recovery)
min_x0 = 0.01*First_Int_Approx
max_x0 = 2.5*First_Int_Approx;
pas_x0 = (max_x0-min_x0)/nb_x0;
x0 = min_x0:pas_x0:max_x0

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


%% Computation of Large AOA Bounds


% Construction of matrix p and t_sched

dz = 0.25;  %length in years bewteen two premium payment dates
p = zeros(Nb_pillars+1,1);
p(1,1) = 1;
t_sched = cell(Nb_pillars);
Period_length = cell(Nb_pillars);

for i=1:Nb_pillars
    matDate = Pillars(i,1);
    coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
    dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
    T = dz1 + dz*round(coupon2mat/dz);
    t_sched_temp = (dz1:dz:T)';
    t_sched_temp(end) = t_pillars(i);
    t_sched{i} = t_sched_temp;    
    Period_length{i} = diff([0;t_sched{i}]);
    p(i+1,1) = length(t_sched{i});
end



% Construction of matrix K
K = zeros(Nb_pillars,Nb_pillars);
for i = 1:Nb_pillars
    S_i = CDS_spreads_AIG(i) + (1-Recovery)*Interest
    t_sched_temp = t_sched{i};
    Period_length_temp = Period_length{i};
    for j = 1:i  
        Periods = Period_length_temp(p(j):(p(j+1)-1));
        Times = t_sched_temp(p(j):(p(j+1)-1));
        Discount = exp(-Interest*Times);
        K(i,j) = S_i*Periods'*Discount;      
    end
end
%K    

% Construction of matrix J
J = zeros(Nb_pillars,1);
for i = 1:Nb_pillars
   S_i_pi = CDS_spreads_AIG(i) + (1-Recovery)*Interest;
   t_sched_temp = t_sched{i};
   Period_length_temp = Period_length{i};
   T_i = t_sched_temp(p(i+1));
   %Period_length(p(i+1))
   J(i) = (Period_length_temp(p(i+1))*S_i_pi + (1-Recovery))*exp(-Interest*T_i);
end

% Construction des bornes AOA larges 
Bounds_Large = zeros(Nb_pillars,2);
Bounds_Sharp = zeros(Nb_pillars,2);

%P_min pour T1
Bounds_Large(1,1) = (1-Recovery-K(1,1))/J(1);
Bounds_Sharp(1,1) = (1-Recovery-K(1,1))/J(1);
%P_max pour T1
Bounds_Large(1,2) = (1-Recovery)/(J(1)+K(1,1));
Bounds_Sharp(1,2) = (1-Recovery)/(J(1)+K(1,1));

hold on;
eps = 0.01;
t_sched_temp = t_sched{1};
T_1 = t_sched_temp(p(1+1));
line([T_1 T_1], [Bounds_Large(1,1) Bounds_Large(1,2)],'LineWidth',2, 'Color', 'k')
line([T_1+eps T_1+eps], [Bounds_Sharp(1,1) Bounds_Sharp(1,2)],'LineWidth',2, 'Color', 'r')

Index_between_1 = p(1):(p(2)-1);
Size_P_extreme_1 = size(Index_between_1);  
P_extreme_1 = ones(Size_P_extreme_1)';
scatter(t_sched_temp(Index_between_1), P_extreme_1, '.r' );

T_i_minus_1 = 0;
line([T_i_minus_1 T_1], [1 1],'LineWidth',1, 'Color', 'r')  

Index_between_2 = (p(1)+1):p(2);
Size_P_extreme_2 = size(Index_between_2);  
P_extreme_2 = Bounds_Sharp(1,2)*ones(Size_P_extreme_2)';
scatter(t_sched_temp(Index_between_2), P_extreme_2, '.b' )

T_i_minus_1 = 0;
line([T_i_minus_1 T_1], [Bounds_Sharp(1,2) Bounds_Sharp(1,2)],'LineWidth',1, 'Color', 'b')


% TEST
% B1 = (1-Recovery-K(1,1))
 t_sched_temp = t_sched{1};
 Period_length_temp = Period_length{1};
 p1 = p(2);
 S_1_j = CDS_spreads_AIG(1) + (1-Recovery)*Interest;
% Times2 = Period_length_temp(1:(p1-1));
% Discount2 = exp(-Interest*t_sched_temp(1:(p1-1)));
% B2 = 1-Recovery - S_1_j*sum(Times2.*Discount2)
A1 = J(1)
A2 = (Period_length_temp(end)*S_1_j + (1-Recovery))*exp(-Interest*t_sched_temp(end))



for i=2:Nb_pillars
%for i=2:2        
    
    %Min : large bounds approach
    Bounds_Large(i,1) = (1-Recovery - K(i,1) - K(i,2:i)*Bounds_Large(1:(i-1),2))/J(i);
    
    %Max : large bounds approach
    Bounds_Large(i,2) = (1-Recovery - K(i,1:(i-1))*Bounds_Large(1:(i-1),1))/(J(i)+K(i,i));
  
    %Min : large bounds approach
    Bounds_Sharp(i,1) = (1-Recovery - K(i,1) - K(i,2:i)*Bounds_Sharp(1:(i-1),1))/J(i);
    %Max : large bounds approach
    Bounds_Sharp(i,2) = (1-Recovery - K(i,1:(i-1))*Bounds_Sharp(1:(i-1),2))/(J(i)+K(i,i));    
    
    t_sched_temp = t_sched{i};
    T_i = t_sched_temp(p(i+1));
    line([T_i T_i], [Bounds_Large(i,1) Bounds_Large(i,2)],'LineWidth',2, 'Color', 'k')
    line([T_i+eps T_i+eps], [Bounds_Sharp(i,1) Bounds_Sharp(i,2)],'LineWidth',2, 'Color', 'r')
    
    
    Index_between_1 = p(i):(p(i+1)-1);
    Size_P_extreme_1 = size(Index_between_1);  
    P_extreme_1 = Bounds_Sharp(i-1,1)*ones(Size_P_extreme_1)';
    scatter(t_sched_temp(Index_between_1), P_extreme_1, '.r' );
    
    T_i_minus_1 = t_pillars(i-1);
    line([T_i_minus_1 T_i], [Bounds_Sharp(i-1,1) Bounds_Sharp(i-1,1)],'LineWidth',1, 'Color', 'r')  
    

    Index_between_2 = (p(i)+1):p(i+1);
    Size_P_extreme_2 = size(Index_between_2);  
    P_extreme_2 = Bounds_Sharp(i,2)*ones(Size_P_extreme_2)';
    scatter(t_sched_temp(Index_between_2), P_extreme_2, '.b' )
    
    T_i_minus_1 = t_pillars(i-1);
    line([T_i_minus_1 T_i], [Bounds_Sharp(i,2) Bounds_Sharp(i,2)],'LineWidth',1, 'Color', 'b')
    
     
end
scatter(t_sched_temp(end), Bounds_Sharp(Nb_pillars,1), '.r' )


%
% hax=axes; 
% T1 = 3; %your point goes here 
% line([T1 T1],get(hax,'YLim'))
% T2 = 5;
% line([T2 T2],get(hax,'YLim'))






    

   