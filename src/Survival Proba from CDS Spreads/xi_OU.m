function Y=xi_OU(s, a)

%Y = s - (1/a)*(1-exp(-a*s)) + y*(1- exp(-a*s));
Y = s - (1/a)*(1-exp(-a*s));

end