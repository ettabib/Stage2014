function Z = zeta_CGMY(s,y,a,sigma,M,Y)



%%Integration numerique
Nb_pas_tps = 1000;
pas = s/Nb_pas_tps;
x=(0:pas:s);
integrand = (M+ sigma*phii(s-x,y,a)).^Y ;
integrand_gauche = integrand(1:(end-1));
integrand_droite = integrand(2:end);

Approx_int = 0.5*pas*sum(integrand_gauche+integrand_droite);


Z= (M^Y)*gamma(-Y)*s - gamma(-Y)*Approx_int;
end


