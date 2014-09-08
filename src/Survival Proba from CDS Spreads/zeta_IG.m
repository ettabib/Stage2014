function Y =zeta_IG(s,y,a,sigma,lambda)
% s=1;
% y=0;
% a=2;
% sigma = 0.01;
% lambda = 0.01;
% Y=zeta_IG(s,y,a,sigma,lambda)

% s = 2;
% y = 0;
% a = 0.5;
% lambda = 2;
% Y =zzeta(s,y,a,lambda)

%sigma=1;
%e=1/a;
%c=4;
%alpha0=1 + sigma*y/lambda;
%gam=(sigma*(1- a*y))/(sigma*(1-a*y) + a*(lambda + sigma*y));

gamma1 = sigma/a;



% param1 = (1-a*y)/(a*(y + lambda));
% param2 = param1*exp(-a*s);
% param3 = 1+gamma0*(1-exp(-a*s));


% dilog(param1)
% dilog(param2)
% 
% 
% Terme1 = log(alpha0)*s
% Terme2 = (1/a)*(dilog(param1))
% Terme3 = -(1/a)*(dilog(param2))
% Terme4 = -(1/a)*log(param3)*log(param2)
% Terme41 = log(param3)
% Terme42 = log(param2)


%Y=log(alpha0)*s + (1/a)*(dilog(param1) - dilog(param2) - log(param3)*log(param2));

%%Integration numerique
Nb_pas_tps = 1000;
pas = s/Nb_pas_tps;
x=(0:pas:s);
f=gamma1*(1-exp(-a*x));
integrand =(2*f +lambda^(2)).^(0.5) -lambda ;
integrand_gauche = integrand(1:(end-1));
integrand_droite = integrand(2:end);

Approx_int = 0.5*pas*sum(integrand_gauche+integrand_droite);


Y= Approx_int;
end


