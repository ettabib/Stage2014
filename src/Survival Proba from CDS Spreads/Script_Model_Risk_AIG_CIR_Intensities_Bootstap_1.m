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

Res = cell(Size_a, Size_s, Size_x0);
%% 

for i = 1:Size_a
    for j = 1:Size_s
        for k = 1:Size_x0
            %a_c = [a(i); c(j)];
            CIR_Param = [x0(k); a(i); s(j)]; %Param with starting point x0
                                             %Bootstraped_Int = Bootstrap_intensities(toDate, coupDate, Pillars, CDS_market_spreads, OU_Levy_Param, Recovery, Interest);

            Bootstraped_Int = Bootstrap_intensities_ver_1(toDate, coupDate, Pillars, CDS_market_spreads, Recovery, Interest, Intensity_model, CIR_Param)

            b = Bootstraped_Int;
            Res{i,j,k} = b;

            %Eta(i, :) = Bootstraped_Int';

            %disp(Index_data.Names(Index_AIG,1));
            %disp(num2str(Eta(i, :)));
        end
        
    end
end



%%Plot Survival Proba
pas_t = 0.1;
t_scheduled=(0:pas_t:10)';
size_t_sched = size(t_scheduled);

fig=figure; 
ColorSet = varycolor(Size_a*Size_s*Size_x0);
set(gca, 'ColorOrder', ColorSet);
hold all;

for i = 1:Size_a
    for j = 1:Size_s
        for k = 1:Size_x0
            b = Res{i,j,k};
            CIR_Param = [x0(k); a(i); s(j)];

            P = Proba_survi_CIR(t_scheduled, b, t_pillars, CIR_Param);
            plot(t_scheduled, P);

            Price = CDS_price(toDate, coupDate, Pillars,  b, OU_Levy_Param, Recovery, Interest);
            %Price.Fair_spread;

            %% Plot Cum_sum
            %Cum_Sum_P = cumsum(pas_t*P)
            %plot(t_scheduled, Cum_Sum_P);

            %% Test de Coherence - CASE Interest = 0
            %         %CDS maturity 3y
            %             dz = 0.25;  %length in years bewteen two premium payment dates
            %             matDate = Pillars(1,1);
            %             coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
            %             dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
            %             T = dz1 + dz*round(coupon2mat/dz);
            %             t_sched = (dz1:dz:T)';
            %             t_sched(end) = t_pillars(1);
            %             %t_scheduled
            %             %numPeriod = size(t_scheduled,1);
            %             Period_length = diff([0;t_sched]);
            %             P2 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param)
            %             Spread_3y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2)
            % 
            %             S1 = CDS_spreads_AIG(1);
            %             P_T1 = P2(end)
            %             P_T1_Inf = (1-Recovery - S1*sum(Period_length(1:(end-1))))/ (S1*Period_length(end)+(1-Recovery))
            %             P_T1_Sup = (1-Recovery)/ (S1*sum(Period_length)+(1-Recovery))
            % 
            %             %plot(t_scheduled, ones(size_t_sched)*P_T1_Inf,'-b')
            %             %plot(t_scheduled, ones(size_t_sched)*P_T1_Sup,'-r')
            %             line([3 3], [P_T1_Inf P_T1_Sup],'LineWidth',2, 'Color', 'k')
            % 
            %     %         Sum_Coeff = Spread_3y*sum(Period_length)+ (1-Recovery);
            %     %         Coeff_T_3_y = Spread_3y*Period_length(end)+(1-Recovery);
            %     %         Rel_Coeff_T_3_y = Coeff_T_3_y/Sum_Coeff
            % 
            %              %CDS maturity 5y
            %             dz = 0.25;  %length in years bewteen two premium payment dates
            %             matDate = Pillars(2,1);
            %             coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
            %             dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
            %             T = dz1 + dz*round(coupon2mat/dz);
            %             t_sched = (dz1:dz:T)';
            %             t_sched(end) = t_pillars(2);
            %             %t_scheduled
            %             %numPeriod = size(t_scheduled,1);
            %             Period_length = diff([0;t_sched]);
            %             P2 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param);
            %             Spread_5y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2);
            % 
            % 
            %             S2 = CDS_spreads_AIG(2);
            %             p1=13;
            %             p2 = length(t_sched);
            %             Sum_up_T1 = sum(Period_length(1:p1));
            %             Sum_T1_T2 = sum(Period_length((p1+1):p2));
            %             Sum_T1_T2_minus_one = sum(Period_length((p1+1):(p2-1)));
            % 
            % 
            %             P_T2 = P2(end)
            %             P_T2_Inf = (1-Recovery - S2*(Sum_up_T1+P_T1_Sup*Sum_T1_T2_minus_one))/ (S2*Period_length(end)+(1-Recovery))
            %             P_T2_Sup = (1-Recovery - S2*Sum_up_T1*P_T1_Inf)/ (S2*Sum_T1_T2+(1-Recovery))
            %             
            %             %P_T2 = P2(end)
            %             P_T2_Inf_2 = (1-Recovery - S2*(Sum_up_T1+P_T1_Inf*Sum_T1_T2_minus_one))/ (S2*Period_length(end)+(1-Recovery))
            %             P_T2_Sup_2 = (1-Recovery - S2*Sum_up_T1*P_T1_Sup)/ (S2*Sum_T1_T2+(1-Recovery))            
            % 
            % 
            %             %plot(t_scheduled, ones(size_t_sched)*P_T2_Inf,'-b')
            %             %plot(t_scheduled, ones(size_t_sched)*P_T2_Sup,'-r')  
            %             line([5 5], [P_T2_Inf P_T2_Sup],'LineWidth',2, 'Color', 'k')
            %             line([5.05 5.05], [P_T2_Inf_2 P_T2_Sup_2],'LineWidth',2, 'Color', 'r')

            %         
            %         %CDS maturity 7y
            %         dz = 0.25;  %length in years bewteen two premium payment dates
            %         matDate = Pillars(3,1);
            %         coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
            %         dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
            %         T = dz1 + dz*round(coupon2mat/dz);
            %         t_sched = (dz1:dz:T)';
            %         t_sched(end) = t_pillars(3);
            %         %t_scheduled
            %         %numPeriod = size(t_scheduled,1);
            %         Period_length = diff([0;t_sched]);
            %         P2 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param);
            %         Spread_7y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2)
            %  
            %         %CDS maturity 10y
            %         dz = 0.25;  %length in years bewteen two premium payment dates
            %         matDate = Pillars(4,1);
            %         coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
            %         dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
            %         T = dz1 + dz*round(coupon2mat/dz);
            %         t_sched = (dz1:dz:T)';
            %         t_sched(end) = t_pillars(4);
            %         %t_scheduled
            %         %numPeriod = size(t_scheduled,1);
            %         Period_length = diff([0;t_sched]);
            %         P2 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param);
            %         Spread_10y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2);
            %         R = sum(P2)
            %         
            %         Sum_Coeff = Spread_10y*sum(Period_length)+ (1-Recovery);
            %         Coeff_T_10_y = Spread_10y*Period_length(end)+(1-Recovery);
            %         Rel_Coeff_T_10_y = Coeff_T_10_y/Sum_Coeff
            
            
            
        end
    end
end


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

eps = 0.01;
t_sched_temp = t_sched{1};
T_1 = t_sched_temp(p(1+1));
line([T_1 T_1], [Bounds_Large(1,1) Bounds_Large(1,2)],'LineWidth',2, 'Color', 'k')
line([T_1+eps T_1+eps], [Bounds_Sharp(1,1) Bounds_Sharp(1,2)],'LineWidth',2, 'Color', 'r')

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



% %TEST T1
% dz = 0.25;  %length in years bewteen two premium payment dates
% matDate = Pillars(1,1);
% coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
% dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
% T = dz1 + dz*round(coupon2mat/dz);
% t_sched2 = (dz1:dz:T)';
% t_sched2(end) = t_pillars(1);
% %t_scheduled
% %numPeriod = size(t_scheduled,1);
% Period_length2 = diff([0;t_sched2]);
% %P2 = Proba_survi_CIR(t_sched2, b, t_pillars, CIR_Param)
% %Spread_3y = (1-Recovery)*(1-P2(end))/sum(Period_length2.*P2)
% 
% S1 = CDS_spreads_AIG(1);
% %P_T1 = P2(end)
% P_T1_Inf = (1-Recovery - S1*sum(Period_length2(1:(end-1))))/ (S1*Period_length2(end)+(1-Recovery))
% P_T1_Sup = (1-Recovery)/ (S1*sum(Period_length2)+(1-Recovery))
% 
% A1 = K(1,1)
% A2 = S1*sum(Period_length2(1:(end-1)))
% 
% B1 = J(1) 
% %S_i_pi = CDS_spreads_AIG(1)+ (1-Recovery)*Interest;;
% %B1 = (Period_length(p(1+1))*S_i_pi + (1-Recovery))*exp(-Interest*T_i)
% %B1 = (Period_length2(end)*S_i_pi + (1-Recovery))*exp(-Interest*T_i)
% B2 = S1*Period_length2(end)+(1-Recovery)
% 

% 
% %TEST T2
% dz = 0.25;  %length in years bewteen two premium payment dates
% matDate = Pillars(2,1);
% coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
% dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
% T = dz1 + dz*round(coupon2mat/dz);
% t_sched2 = (dz1:dz:T)';
% t_sched2(end) = t_pillars(2);
% %t_scheduled
% %numPeriod = size(t_scheduled,1);
% Period_length2 = diff([0;t_sched2]);
% %P2 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param);
% %Spread_5y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2);
% 
% S2 = CDS_spreads_AIG(2);
% p1=13;
% p2 = length(t_sched2);
% Sum_up_T1 = sum(Period_length2(1:p1));
% Sum_T1_T2 = sum(Period_length2((p1+1):p2));
% Sum_T1_T2_minus_one = sum(Period_length2((p1+1):(p2-1)));
% 
% 
% %P_T2 = P2(end)
% P_T2_Inf = (1-Recovery - S2*(Sum_up_T1+P_T1_Sup*Sum_T1_T2_minus_one))/ (S2*Period_length2(end)+(1-Recovery))
% P_T2_Sup = (1-Recovery - S2*Sum_up_T1*P_T1_Inf)/ (S2*Sum_T1_T2+(1-Recovery))
% 
% %P_T2 = P2(end)
% P_T2_Inf_2 = (1-Recovery - S2*(Sum_up_T1+P_T1_Inf*Sum_T1_T2_minus_one))/ (S2*Period_length2(end)+(1-Recovery))
% P_T2_Sup_2 = (1-Recovery - S2*Sum_up_T1*P_T1_Sup)/ (S2*Sum_T1_T2+(1-Recovery))            
% 
% 
% %plot(t_scheduled, ones(size_t_sched)*P_T2_Inf,'-b')
% %plot(t_scheduled, ones(size_t_sched)*P_T2_Sup,'-r')  
% line([5 5], [P_T2_Inf P_T2_Sup],'LineWidth',2, 'Color', 'k')
% line([5.05 5.05], [P_T2_Inf_2 P_T2_Sup_2],'LineWidth',2, 'Color', 'r')

for i=2:Nb_pillars
    %for i=2:2        
    
    %Min : large bounds approach
    Bounds_Large(i,1) = (1-Recovery - K(i,1) - K(i,2:i)*Bounds_Large(1:(i-1),2))/J(i);
    
    
    
    
    %      A1 = (1-Recovery - K(i,1) - K(i,2:i)*Bounds_Large(1:(i-1),2))
    %      A2 = (1-Recovery - S2*(Sum_up_T1+P_T1_Sup*Sum_T1_T2_minus_one))
    
    %     B1 = J(i);
    %     B2 = (S2*Period_length2(end)+(1-Recovery));
    
    %Max : large bounds approach
    Bounds_Large(i,2) = (1-Recovery - K(i,1:(i-1))*Bounds_Large(1:(i-1),1))/(J(i)+K(i,i));
    
    %Min : large bounds approach
    Bounds_Sharp(i,1) = (1-Recovery - K(i,1) - K(i,2:i)*Bounds_Sharp(1:(i-1),1))/J(i);
    %Max : large bounds approach
    Bounds_Sharp(i,2) = (1-Recovery - K(i,1:(i-1))*Bounds_Sharp(1:(i-1),2))/(J(i)+K(i,i));    
    
    %t_scheduled
    %numPeriod = size(t_scheduled,1);
    %Period_length = diff([0;t_sched]);
    %P4 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param);
    %Spread_10y = (1-Recovery)*(1-P4(end))/sum(Period_length.*P4);


    %     P_T2 = P2(end)
    %     P_T2_Inf = (1-Recovery - S2*(Sum_up_T1+P_T1_Sup*Sum_T1_T2_minus_one))/ (S2*Period_length(end)+(1-Recovery))
    %     P_T2_Sup = (1-Recovery - S2*Sum_up_T1*P_T1_Inf)/ (S2*Sum_T1_T2+(1-Recovery))
    % 
    %     %P_T2 = P2(end)
    %     P_T2_Inf_2 = (1-Recovery - S2*(Sum_up_T1+P_T1_Inf*Sum_T1_T2_minus_one))/ (S2*Period_length(end)+(1-Recovery))
    %     P_T2_Sup_2 = (1-Recovery - S2*Sum_up_T1*P_T1_Sup)/ (S2*Sum_T1_T2+(1-Recovery))            


    %plot(t_scheduled, ones(size_t_sched)*P_T2_Inf,'-b')
    %plot(t_scheduled, ones(size_t_sched)*P_T2_Sup,'-r')  
    
    t_sched_temp = t_sched{i};
    T_i = t_sched_temp(p(i+1));
    line([T_i T_i], [Bounds_Large(i,1) Bounds_Large(i,2)],'LineWidth',2, 'Color', 'k')
    line([T_i+eps T_i+eps], [Bounds_Sharp(i,1) Bounds_Sharp(i,2)],'LineWidth',2, 'Color', 'r')
    
end


%
% hax=axes; 
% T1 = 3; %your point goes here 
% line([T1 T1],get(hax,'YLim'))
% T2 = 5;
% line([T2 T2],get(hax,'YLim'))








