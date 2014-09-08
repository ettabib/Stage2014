%% Plot Survival distr function for a Gamma-Levy OU intensity model
% mean-reversion parameter is a piecewise-constant function of time


bp = 0.0001;

%Lambda is fixed at 1/Taille_Saut
Taille_Saut = 10*bp;
lambda = 1/Taille_Saut;
%Sigma is fiwed at 1
sigma = 1;

% definition grille a (vitesse retour moyenne) et grille c
nb_a = 10;
nb_c = 5;
nb_s = 10; %nb sigma param

min_a = 0.01;
max_a = 2;
pas_a = (max_a-min_a)/nb_a;

min_s = 0.1;
max_s = 6;
pas_s = (max_s-min_s)/nb_s;

min_c = 0.1;
max_c = 10;
pas_c = (max_c-min_c)/nb_c;

a = min_a:pas_a:max_a
s = min_s:pas_s:max_s
c = min_c:pas_c:max_c

a = [1];
%s = [1];
c = [1];

Size_a = size(a,2);
Size_s = size(s,2);
Size_c = size(c,2);



%% Assumption on recovery rate and interest rate
Recovery = 0.40;
%Interest = 0.03;
Interest = 0;

Intensity_model = 2; %Gamma-OU

%%  Bootstrap procedure
load('CDX9_17_Dec_2007.mat'); %Load Index_data containing all relevant information regarding this date
% and in particular CDS spreads of all consituents at different pillars

toDate = Index_data.Today;                                   % today
coupDate = Index_data.coupDate;                             % first coupon payment date
Pillars = Index_data.CDS_Pillars;

Index_AIG = 5;
CDS_spreads_AIG = bp*Index_data.CDS_Spreads(Index_AIG,:)

Nb_pillars = size(Pillars, 1);

%Convert cell vector of dates pillars into vectors of real number t_pillars
%corresponding to lengths between today and pillars dates 
t_pillars = zeros(Nb_pillars, 1);
for i=1:Nb_pillars
    t_pillars(i) = (datenum(Pillars(i),'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
end


CDS_market_spreads = CDS_spreads_AIG';

Res = cell(Size_a, Size_s, Size_c);


for i = 1:Size_a
    for j = 1:Size_s
        for k = 1:Size_c
            %a_c = [a(i); c(j)];
            OU_Levy_Param = [a(i); s(j); c(k); lambda];
            %Bootstraped_Int = Bootstrap_intensities(toDate, coupDate, Pillars, CDS_market_spreads, OU_Levy_Param, Recovery, Interest);

            Bootstraped_Int = Bootstrap_intensities_ver_2(toDate, coupDate, Pillars, CDS_market_spreads, Recovery, Interest, Intensity_model, OU_Levy_Param)

            b = Bootstraped_Int;
            Res{i,j,k} = b;

            %Eta(i, :) = Bootstraped_Int';

            %disp(Index_data.Names(Index_AIG,1));
            %disp(num2str(Eta(i, :)));
        end
    end
end



%%Plot Survival Proba
t_scheduled=(0:0.1:10)';
hold on;
for i = 1:Size_a
    for j = 1:Size_s
        for k = 1:Size_c
            b = Res{i,j};
            OU_Levy_Param = [b(1); a(i); s(j); c(k); lambda];

            P = Proba_survi_OU_Gamma(t_scheduled, b, t_pillars, OU_Levy_Param);
            plot(t_scheduled, P);
            %Price = CDS_price(toDate, coupDate, Pillars,  b, OU_Levy_Param, Recovery, Interest);
            %Price.Fair_spread;

              %CDS maturity 3y
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
%             P2 = Proba_survi_OU_Gamma(t_sched, b, t_pillars, OU_Levy_Param);
%             Spread_3y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2)
%             S1 = CDS_spreads_AIG(1)
%             
%             %CDS maturity 5y
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
%             Spread_5y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2)
%             S2 = CDS_spreads_AIG(2)
            
            
        end        
    end
end


%% Test de Coherence - CASE Interest = 0
%         %CDS maturity 3y
dz = 0.25;  %length in years bewteen two premium payment dates
matDate = Pillars(1,1);
coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
T = dz1 + dz*round(coupon2mat/dz);
t_sched = (dz1:dz:T)';
t_sched(end) = t_pillars(1);
%t_scheduled
%numPeriod = size(t_scheduled,1);
Period_length = diff([0;t_sched]);
% P2 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param)
% Spread_3y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2)

S1 = CDS_spreads_AIG(1);
P_T1 = P2(end)
P_T1_Inf = (1-Recovery - S1*sum(Period_length(1:(end-1))))/ (S1*Period_length(end)+(1-Recovery))
P_T1_Sup = (1-Recovery)/ (S1*sum(Period_length)+(1-Recovery))

%plot(t_scheduled, ones(size_t_sched)*P_T1_Inf,'-b')
%plot(t_scheduled, ones(size_t_sched)*P_T1_Sup,'-r')
line([3 3], [P_T1_Inf P_T1_Sup],'LineWidth',2, 'Color', 'k')

%         Sum_Coeff = Spread_3y*sum(Period_length)+ (1-Recovery);
%         Coeff_T_3_y = Spread_3y*Period_length(end)+(1-Recovery);
%         Rel_Coeff_T_3_y = Coeff_T_3_y/Sum_Coeff

 %CDS maturity 5y
dz = 0.25;  %length in years bewteen two premium payment dates
matDate = Pillars(2,1);
coupon2mat = (datenum(matDate,'dd-mmm-yyyy')-datenum(coupDate,'dd-mmm-yyyy'))/365;
dz1 = (datenum(coupDate,'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
T = dz1 + dz*round(coupon2mat/dz);
t_sched = (dz1:dz:T)';
t_sched(end) = t_pillars(2);
%t_scheduled
%numPeriod = size(t_scheduled,1);
Period_length = diff([0;t_sched]);
% P2 = Proba_survi_CIR(t_sched, b, t_pillars, CIR_Param);
% Spread_5y = (1-Recovery)*(1-P2(end))/sum(Period_length.*P2);


S2 = CDS_spreads_AIG(2);
p1=13;
p2 = length(t_sched);
Sum_up_T1 = sum(Period_length(1:p1));
Sum_T1_T2 = sum(Period_length((p1+1):p2));
Sum_T1_T2_minus_one = sum(Period_length((p1+1):(p2-1)));


P_T2 = P2(end)
P_T2_Inf = (1-Recovery - S2*(Sum_up_T1+P_T1_Sup*Sum_T1_T2_minus_one))/ (S2*Period_length(end)+(1-Recovery))
P_T2_Sup = (1-Recovery - S2*Sum_up_T1*P_T1_Inf)/ (S2*Sum_T1_T2+(1-Recovery))


%plot(t_scheduled, ones(size_t_sched)*P_T2_Inf,'-b')
%plot(t_scheduled, ones(size_t_sched)*P_T2_Sup,'-r')  
line([5 5], [P_T2_Inf P_T2_Sup],'LineWidth',2, 'Color', 'k')










    

   