function result = phi_CIR_2(s, param)
%function result = phi(s, y, param)
% T1 = 11;
% t  = 27;
% res1 = phi(t,0)
% res2 = phi(T1, phi(t-T1, 0))

a = param(1); %eta = a
sigma = param(2); %nu = sigma
h = sqrt(a^2+2*sigma^2);
hpa = h + a;
hma = h - a;

result = 2*(1-exp(-h*s))./(hpa + hma*exp(-h*s));

% B=(eta+sqrt(eta*eta+2*nu*nu))/2;
% C=(eta-sqrt(eta^2+2*nu*nu))/(-2);
% D=-1;
% A=(D*(nu*nu+eta*B)-C*(2*B-eta))/(B*D-C);
% result = (1+D*exp(-A*s))./(B+C*exp(-A*s));

% B=(eta+sqrt(eta*eta+2*nu*nu))/2;
% C=(1-B*y)*(eta+nu*nu*y-sqrt(eta^2+2*nu*nu))/(2*eta*y+nu*nu*y*y-2);
% D=(B+C)*y-1;
% A=(D*(nu*nu+eta*B)-C*(2*B-eta))/(B*D-C);
% result = (1+D*exp(-A*s))/(B+C*exp(-A*s));