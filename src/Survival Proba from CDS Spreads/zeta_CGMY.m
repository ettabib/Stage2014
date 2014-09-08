function Z = zeta_CGMY(s,y,a,sigma,M,Y)



%%Integration numerique




Approx_int =quad(@(x) (M+ sigma*phii(s-x,y,a)).^Y ,0,s);


Z= (M^Y)*gamma(-Y)*s - gamma(-Y)*Approx_int;
end


