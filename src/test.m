bp = 0.0001;
% Recovery
Recovery = 0.40;
% Interest
Interest = 0.03;

%% Loading data AIG spreads 4 maturities
% CDS_spreads_AIG = bp*Index_data.CDS_Spreads(Index_AIG,:);    %spreads
% T  = Index_data.CDS_Pillars_in_years';
% n  = size(T,2);
% T  = [ 0 ; T' ];
% S  = CDS_spreads_AIG;
% S  = [ 0 ; S ]' ;

%% Import the data
[~, ~, raw] = xlsread('/Users/mohammad/Documents/MyProject/Stage2014/data/France_1.xlsx','Sheet1');
raw = raw(2:end,:);
%Dates = cell2mat(raw(:,1)) + datenum('01-Jan-1904');
%raw(:,1) = {Dates}; 

%raw(:,datecol) = raw(:,datecol) + datenum('01-Jan-1904');

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Loading data AIG spreads 10 maturities
T  =  1:(size(data,2)-1);
n  = size(T,2);
T  = [ 0 ; T' ];
S  = data(1,(2:11)) * bp;

%% Test
DATE=596;
[Qmin,Qmax] = Proba_Survie(1:10,data(DATE,(2:11))' * bp,Recovery,Interest);
Y = [Qmin(2:(n+1)) , Qmax(1:n)];
figure
stairs(Y)

%% Date where There is AOA
DateAOA = 0;
AOA = false;
for i=2:size(data,1),
   [Qmin,Qmax] = Proba_Survie(1:10,data(i,(2:11))' * bp,Recovery,Interest);
   if((Qmin(2:(n+1)) <= Qmax(1:n))),        
       if(not(AOA)),           
            DateAOA = i; 
       end
       AOA = true;      
   else
       AOA = false ;
   end
end
DateAOA = DateAOA + 1;

%% Scatter for the cotation 18/08/2005
[Qmin,Qmax] = Proba_Survie(1:10,data(4,(2:11))' * bp,Recovery,Interest);
Y = [Qmin(2:(n+1)) , Qmax(1:n)];
figure
stairs(Y)
title('Scatter for the cotation 18/08/2005');
legend('Qmin','Qmax');

%% Scatter for the cotation 27/11/2007
[Qmin,Qmax] = Proba_Survie(1:10,data(597,(2:11))' * bp,Recovery,Interest);
Y = [Qmin(2:(n+1)) , Qmax(1:n)];
figure
stairs(Y)
title('Scatter for the cotation 27/11/2007');
legend('Qmin','Qmax');

%% Curves evolution
TabCurveQmax = zeros(size(data,1)-1,11);
TabCurveQmin = zeros(size(data,1)-1,11);
TabDiff = zeros(size(data,1)-1,11);
NormQmax     = zeros(size(data,1));
% figure
for i=2:(size(data,1)),
   [Qmin,Qmax] = Proba_Survie(1:10,data(i,(2:11))' * bp,Recovery,Interest);
   TabCurveQmax(i-1,:) = Qmax;
   TabCurveQmin(i-1,:) = Qmin;
   TabDiff(i-1,:)      = Qmax - Qmin; 
%    Y = [Qmin(2:(n+1)) , Qmax(1:n)];
%    stairs(Y); 
    %NormQmax(i-1) = abs(sum(abs(Qmax(1:10)-Qmax(2:11))) - mean(Qmax)); 
%     NormQmax(i-1) =  var(Qmax,0,1);
%      NormQmax(i-1) =  mean(Qmax,1);
    NormQmax(i-1) = mean(Qmax - Qmin);
%     NormQmax(i-1) = mean(data(i,(2:11)));
end
% Define the norm : N(C) = sum d - mean (d)

% MeanQmax = mean(TabCurveQmax,1);  
% MeanQmin = mean(TabCurveQmin,1);
% Y = [MeanQmin(2:(n+1)) , MeanQmax(1:n)];

figure
% plot(TabCurveQmax);
plot(NormQmax(1:1336));
%% Ploting all spreads
plot(var(data(:,2:11),0,2))