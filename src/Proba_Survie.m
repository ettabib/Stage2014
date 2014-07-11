function [ Qmin , Qmax ] = Proba_Survie( Maturity , Spread , Recovery , Interest)

T0 = 0;
tp0 = T0;
n  = size(Maturity,2);
T  = [ 0 ; Maturity' ];
S  = [ 0 ; Spread ]' ;
d  = 1/4;
pas =  d;
t  = tp0:pas:T(end);
p  = zeros(n+1,1);

for i=0:n,
    I = i + 1;
    p(I) = find(t == T(I));
end

M=zeros(n,1);N=zeros(n,1);Qmin=zeros(n,1);Qmax=zeros(n,1);
Qmax(1)=1;Qmin(1)=1;

for i=1:n,
    I = i + 1;
    
    M(I) = Pd(T0,T(I - 1),Interest) - Pd(T0,T(I),Interest);
    N(I) = d * sum(Pd(T0,t(p(I - 1) : (p(I) - 1)),Interest));
    
    K    = 2 : I ;
    V    = (1 - Recovery) .* M(K) + S(I) .* N(K);
    
    Km1    = 1:(I - 1);
    Qmin(I) = 1 - Recovery - sum(V .* Qmax(Km1));
    Qmin(I) = Qmin(I) / (Pd(T0,T(I),Interest) * (1 - Recovery + S(I) * d));
    
    K = 2:(I - 1);
    V    = (1 - Recovery) * M(K) + S(I) * N(K);
    
    Qmax(I) = 1 - Recovery - sum(V .* Qmin(K));
    Qmax(I) = Qmax(I) / (Pd(T0,T(I - 1),Interest) * (1 - Recovery) + S(I) * (N(I) + d * Pd(T0,T(I),Interest)));
    
   % assert(Qmax(I-1)>=Qmin(I));   
    
end

end


