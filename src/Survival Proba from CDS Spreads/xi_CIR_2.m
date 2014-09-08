function result = xi_CIR_2(s,param)
%function result = xi(s, y,param)
% T1 = 11;
% t  = 27;
% res1 = xi(t,0)
% res2 = xi(T1,phi(t-T1, 0)) + xi(t-T1,0)


a = param(1); %eta = a
sigma = param(2); %nu = sigma
h = sqrt(a^2+2*sigma^2);
hpa = h + a;
hma = h - a;

Q2 = (hpa+hma*exp(-h*s))/(2*h);

result = 2*a*((s/hpa) + (1/sigma^2)*log(Q2));

% B=(eta+sqrt(eta*eta+2*nu*nu))/2;
% C=(eta-sqrt(eta^2+2*nu*nu))/(-2);
% D=-1;
% A=(D*(nu*nu+eta*B)-C*(2*B-eta))/(B*D-C);
% Q1=(C-B*D)/(A*C);
% Q2=(B+C*exp(-A*s))/(B+C);
% result = eta*(Q1*log(Q2)+s)./B;

% B=(eta+sqrt(eta*eta+2*nu*nu))/2;
% C=(1-B*y)*(eta+nu*nu*y-sqrt(eta^2+2*nu*nu))/(2*eta*y+nu*nu*y*y-2);
% D=(B+C).*y-1;
% A=(D*(nu*nu+eta*B)-C*(2*B-eta))/(B*D-C);
% Q1=(C-B*D)/(A*C);
% Q2=(B+C*exp(-A*s))/(B+C);
% result = eta*(Q1*log(Q2)+s)/B;