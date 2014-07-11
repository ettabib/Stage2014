bp = 0.0001;
% Recovery
Recovery = 0.40;
% Interest
Interest = 0.03;

%% Loading data AIG spreads 4 maturities
CDS_spreads_AIG = bp*Index_data.CDS_Spreads(Index_AIG,:);    %spreads
T  = Index_data.CDS_Pillars_in_years';
n  = size(T,2);
T  = [ 0 ; T' ];
S  = CDS_spreads_AIG;
S  = [ 0 ; S ]' ;

%% Import the data
[~, ~, raw] = xlsread('/Users/mohammad/Documents/MyProject/Stage2014/data/France_1.xlsx','Sheet1');

raw = raw(2:end,:);

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Loading data AIG spreads 10 maturities
T  =  1:(size(data,2)-1);
n  = size(T,2);
T  = [ 0 ; T' ];
S  = data(1,(2:11)) * bp;

%% Test
%[Qmin,Qmax] = Proba_Survie(Index_data.CDS_Pillars_in_years',[ 58 ; 54 ; 52 ; 49 ] * bp,Recovery,Interest);
[Qmin,Qmax] = Proba_Survie(1:10,data(1145,(2:11))' * bp,Recovery,Interest);
Y = [Qmin(2:(n+1)) , Qmax(1:n)];
figure
stairs(Y)


%% Date where The marke fit AOA
for i=2:size(data,1),
   [Qmin,Qmax] = Proba_Survie(1:10,data(i,(2:11))' * bp,Recovery,Interest);
   if((Qmin(2:(n+1)) <= Qmax(1:n)))
       fprintf('first date of AOA is %1.0f \n',i);       
       break 
   end
end

