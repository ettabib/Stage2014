%% Testing distrToX function

lambda = 1/10;
t = 0:0.001:1000;
F = 1 - exp(- lambda * t );
% a = 1;
% t0 = -10;
%F = 1 / pi * atan( (t - t0) / a) + 0.5;
D = 0;
for i=1:1000
   D(i) = distrToX(F,t);
end

histfit(D);