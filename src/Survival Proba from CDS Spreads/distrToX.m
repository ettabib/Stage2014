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
u = rand() * min(max(F),1) + max(min(F),0);
X = t(find(F==max(F(F < u)))+1);
end