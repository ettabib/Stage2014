bp = 0.0001;
% Recovery
Recovery = 0.40;
% Interest
Interest = 0.03;

%% 

load('data/CDX9_17_Dec_2007.mat')

toDate = Index_data.Today;                                   % today
coupDate = Index_data.coupDate;                              % first coupon payment date
Pillars = Index_data.CDS_Pillars;

Index_AIG = 5;
CDS_spreads_AIG = bp*Index_data.CDS_Spreads(Index_AIG,:);    %spreads


t0 = 0;
m  = 100;
T  = Index_data.CDS_Pillars_in_years';
n  = size(T,2);
d  = 1/4;
t  = 0:T(end); 
p  = zeros(n,1);
for i=1:n,
    p(i) = find(t == T(i));
end
S  = CDS_spreads_AIG;

%% Loading datas
%DATA = xlsread('/Users/mohammad/Documents/MyProject/Stage2014/data/France.xlsx');


%Preallocation
M=zeros(n,1);N=zeros(n,1);Qmin=zeros(n,1);Qmax=zeros(n,1);
%% Computing recursively 

for i=1:n,
    if(i==1),
        M(1) = 1 - Pd(t0,T(1),Interest);
        N(1) = d*sum(Pd(t0,t(1:(p(1)-1)),Interest));
    else
        M(i) = Pd(t0,T(i - 1),Interest) - Pd(t0,T(i),Interest);
        N(i) = d * sum(Pd(t0,t(p(i - 1) : (p(i) - 1)),Interest));    
    end
    
    V    = (1 - Recovery) * M(1 : i) + S(i) * N(1 : i);
    
    Qmin(T(i)) = 1 - Recovery - sum(V .* Qmax([1 ; (1:(i - 1))']));
    if(i==1),
        Qmin(T(i)) = Qmin(T(i)) /  (1 - Recovery + S(i) * d);
    else
        Qmin(T(i)) = Qmin(T(i)) / (Pd(t0,T(i - 1),Interest) * (1 - Recovery + S(i) * d));
    end
    
    V    = (1 - Recovery) * M(1:(i - 1)) + S(i) * N(1:(i - 1));
    
    Qmax(T(i)) = 1 - Recovery - sum(V .* Qmin(T(1:(i - 1))));
    if(i==1),
        Qmax(T(i)) = Qmax(T(i)) / ((1 - Recovery) + S(i) .* (N(i) + d .* Pd(t0,T(i),Interest)));
    else
        Qmax(T(i)) = Qmax(T(i)) / (Pd(t0,T(i - 1),Interest) .* (1 - Recovery) + S(i) .* (N(i) + d .* Pd(t0,T(i),Interest)));
    end
end
