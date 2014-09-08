function Y=zeta_Gamma_exact(s,y,a,sigma,lambda)
% s=1;
% y=0;
% a=2;
% sigma = 0.01;
% lambda = 0.01;
% Y=zeta_Gamma_exact(s,y,a,sigma,lambda)

% s=6;
% y=0;
% a=2;
% sigma = 0.05;
% lambda = 1;
% Y=zeta_Gamma_exact(s,y,a,sigma,lambda)

% calcul de façon exacte la fonction zeta dans le modele Gamma OU

alpha0= 1+ (sigma*y)/lambda;

beta0=(sigma*(1-a*y))/(sigma*(1-a*y) + a*(lambda + sigma*y));


gamma0=(sigma*(1-a*y))/(a*(lambda + sigma*y));

coef=1/a;

int1=dilog(1-beta0);

int2=dilog(1-beta0*exp(-a*s)) + log( 1 + gamma0*(1-exp(-a*s)))*log(beta0*exp(-a*s));

integrale=int1 -int2;


Y=log(alpha0)*s  + coef*integrale;




end