%% Bootstrap procedure to find piece-wise constant parameters associated with a generalized OU-Levy model or a CIR model
% We consider an intensity model with an extended OU-Levy intensity
% process:)
% model dX_t = a(b(t)-X_t)dt + \sigma X_t^{\alpha}dL_{ct} with X0=x
% We assume that b(t) is a piecewise constant function of time with a grid
% point corresponding to CDS pillars. The aim is to bootstrap the
% piecewise-constant function b(t) from a CDS curve.
% Parameters a and c are fixed a priori
% \alpha is either 0 or 1/2
% When \alpha = 0 : L is a Lévy-subordinator (OU process driven by a Levy
% subordinator)
% When \alpha = 1/2 : L is a brownian motion (CIR process)


%% Assumption on recovery rate and interest rate
Recovery = 0.40;
Interest = 0.03;
bp = 0.0001;

%% Choice of the underlying stochastic intensity model

% Intensity_model = 1; %CIR process
% Intensity_model = 2; %Gamma-Levy OU process
% Intensity_model = 3; %IG-Levy OU process
% Intensity_model = 4; %CGMY-Levy OU process

% CIR process
% Intensity_model = 1; 
% disp(' ');
% disp('--------- CIR version of the model ---------');
% disp(' ');
% a = 3 
% sigma = 0.05
% Param = [a; sigma]; %Parametre du CIR sans le point initial

% Gamma-Levy OU process
Intensity_model = 2;
disp(' ');
disp('--------- Gamma-OU version of the model ---------');
disp(' ');
a = 3;
sigma = 500*bp;
%sigma = 1;
c = 100 %Nb moyen de saut par an
%c = 20 %Nb moyen de saut par an
%Taille_Saut = 3*bp;
%Taille_Saut = 1000*bp;
Taille_Saut = 100*bp;
lambda = 1/Taille_Saut;
Param = [a; sigma; c; lambda]; %Parametre du Gamma-OU sans le point initial


% IG-Levy OU process
% Intensity_model = 3;
% disp(' ');
% disp('--------- IG-OU version of the model ---------');
% disp(' ');
% a = 3
% sigma = 500*bp
% %sigma = 1;
% c = 20 %Nb moyen de saut par an
% %Taille_Saut = 3*bp;
% Taille_Saut = 100*bp;
% lambda = 1/Taille_Saut
% Param = [a; sigma; c; lambda]; %Parametre du Gamma-OU sans le point initial

%CGMY-Levy OU process
% Intensity_model = 4;
% disp(' ');
% disp('--------- CGMY-OU version of the model ---------');
% disp(' ');
% a = 3
% sigma = 500*bp
% %sigma = 1;
% c = 20 %Nb moyen de saut par an
% M = 10; %M positive parameter
% Y = -10; %Y<1
% Param = [a; sigma; c; M; Y]; %Parametre du Gamma-OU sans le point initial


%% Choice of the bootstrap method
%Boostrap = 1; %The starting point of the intensity process is set at the constant intensity that best fit the shortest pillar spread
Boostrap = 2; %The starting point of the CIR proces is fixed at the first mean-reversion level and bootstapped on the shortest pillar spread

%% Load General CDX Index Data
load('CDX9_17_Dec_2007_OU_Levy.mat'); %Load Index_data containing all relevant information regarding this date
% and in particular CDS spreads of all consituents at different pillars

toDate = Index_data.Today;                                   % today
coupDate = Index_data.coupDate;                             % first coupon payment date
Pillars = Index_data.CDS_Pillars;
All_CDS_spreads = bp*Index_data.CDS_Spreads;

Nb_pillars = size(Pillars, 1);

%Convert cell vector of dates pillars into vectors of real number t_pillars
%corresponding to lengths between today and pillars dates 
t_pillars = zeros(Nb_pillars, 1);
for i=1:Nb_pillars
    t_pillars(i) = (datenum(Pillars(i),'dd-mmm-yyyy')-datenum(toDate,'dd-mmm-yyyy'))/365;
end



%% Bootstrap intensities for all CDS in the index
%Nb_CDS = size(All_CDS_spreads,1);
Index_CDS = 1:5;
Nb_CDS = size(Index_CDS,2);
%Eta : Vector of bootsraped intensities
Eta = zeros(Nb_CDS, Nb_pillars);



if (Boostrap == 1)
    
  

    %%  The starting point of the underlying process is set at the constant that best fit the shortest pillar spread
    % It minimizes the least square error between model spreads and market
    % spreads given that the underlying intensity is assumed to be constant 
    %Load Index_data, i.e., the object containing relevant information on
    %CDX9_17_Dec_2007_First_Pillar
    disp(' ');    
    disp(' ------------- Bootstrap version 1 --------------');
    disp(' ');
    disp('Starting point fixed at the implied constant intensity');
    disp(' ');   

    load('CDX9_17_Dec_2007_First_Pillar.mat'); %Load Starting point of individual CIR intensity processes
    Constant_Intensities = Index_data.Eta;

    for i=1:Nb_CDS  
        CDS_market_spreads = All_CDS_spreads(Index_CDS(i),:)'; %On recupere la courve de credit associe au nom i
        Param_with_starting_point = [Constant_Intensities(Index_CDS(i)); Param];   %On rajoute le premier parametre (point initiale) par celui obtenu apres une calibration à intensite constante
        Bootstraped_Int = Bootstrap_intensities_ver_1(toDate, coupDate, Pillars, CDS_market_spreads, Recovery, Interest, Intensity_model, Param_with_starting_point );
        Eta(Index_CDS(i), :) = Bootstraped_Int';
        disp(Index_data.Names(Index_CDS(i),1));
        disp(num2str(Eta(Index_CDS(i), :)));
        disp(' ');
        %disp(Constant_Intensities(i));
    end
    
    res.Initial_point = Constant_Intensities(:,1);
    
else 

    %%  The starting point of the CIR process is set at the first mean-reversion coeff and bootstrapped from the first pillar spread
    disp(' ');   
    disp(' ------------- Bootstrap version 2 --------------');
    disp(' ');
    disp('Starting point equal to first mean-reversion coeff and bootstrapped from the first pillar spread');
    disp(' ');   
    
    for i=1:Nb_CDS  
        CDS_market_spreads = All_CDS_spreads(Index_CDS(i),:)';
        Bootstraped_Int = Bootstrap_intensities_ver_2(toDate, coupDate, Pillars, CDS_market_spreads, Recovery, Interest, Intensity_model, Param);
        Eta(Index_CDS(i), :) = Bootstraped_Int';
        disp(Index_data.Names(Index_CDS(i),1));
        disp(num2str(Eta(Index_CDS(i), :)));
        disp(' ');
    end
    
    res.Initial_point = Eta(:,1); %starting point equal to the first mean-reversion param
    
end
    
res.Eta = Eta;
res.CDS_pillars = t_pillars';


%% PLOT Implied AIG Intensity paths
AIG = 5;
b1 =  res.Eta(5,1);
b2 =  res.Eta(5,2);
b3 =  res.Eta(5,3);

b_scheme = res.Eta(AIG,:);

if (Boostrap == 1)
    Param_with_starting_point = [Constant_Intensities(AIG); Param];
else
    Param_with_starting_point = [b1; Param];
end

% %TEST GAMMA COMPENSATED
Compensator = Param_with_starting_point(3)*Param_with_starting_point(4)/(Param_with_starting_point(2)*Param_with_starting_point(5))
b_scheme = b_scheme - Compensator

NbPath = 40;
Nb_time_step = 1000;
Plot_Path(NbPath, Nb_time_step, b_scheme, t_pillars, Param_with_starting_point, Intensity_model)

