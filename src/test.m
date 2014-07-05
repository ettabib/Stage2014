bp = 0.0001;
% Recovery
Recovery = 0.40;
% Interest
Interest = 0.03;



%% 

%load('data/CDX9_17_Dec_2007.mat')

CDS_spreads_AIG = bp*Index_data.CDS_Spreads(Index_AIG,:);    %spreads

T  = Index_data.CDS_Pillars_in_years';
%T  =  1:10;
n  = size(T,2);
T  = [ 0 ; T' ];
S = [ 49 ; 52 ; 58 ; 54 ] * bp ;
% S  = CDS_spreads_AIG;
%S  = data(44,(2:11)) * bp;
S  = [ 0 ; S ]' ;


%% Test
%[Qmin,Qmax] = Proba_Survie(Index_data.CDS_Pillars_in_years',[ 58 ; 54 ; 52 ; 49 ] * bp,Recovery,Interest);
[Qmin,Qmax] = Proba_Survie(1:10,data(1308,(2:11))' * bp,Recovery,Interest);

