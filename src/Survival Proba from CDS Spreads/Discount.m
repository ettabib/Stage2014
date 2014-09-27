function [ y ] = Discount( t , T , R)
% computing the discount factor
%input :
%   t : date of discount    
%   R : the rate
%   T : the maturity
y = exp(- R * (T - t));

end

