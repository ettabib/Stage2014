function [ res ] = Pd( t0,T,rs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    res = exp(-rs * (T - t0));

end

