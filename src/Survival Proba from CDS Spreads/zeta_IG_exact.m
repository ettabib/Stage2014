function Y=zeta_IG_exact(s,y,a,sigma,lambda)
% s=1;
% y=0;
% a=2;
% sigma = 0.01;
% lambda = 0.01;
% Y=zeta_IG_exact(s,y,a,sigma,lambda)



%calcule zeta de façon exacte dans le modele IG OU

coef= 1/a;

lametoile=lambda^2 + 2*sigma*y;

alpha=2*sigma*(1/a -y);

coef1=sqrt(lametoile + alpha);

coef2=sqrt(lametoile + alpha*(1-exp(-a*s)));

int1=sqrt(lametoile) -coef1*atanh(sqrt(lametoile)/coef1);

int2=coef2 -coef1*atanh(coef2/coef1);

integrale=int1 - int2;


Y=2*coef*integrale - lambda*s;


end