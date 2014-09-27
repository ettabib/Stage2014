%% simulation of random variable from distribution function
% inverse algorithm 
% input :
%           F  : the distribution fuction
%           t  : time discretisation, same size as F
% output :
%            a random simulation for the distribution F

function [ X ] = distrToX(F,t)
i = 1;
L = size(F,2);
u = rand();

while (F(i) < u  && i < L )
    i = i+1;
end

X = t(i);

end