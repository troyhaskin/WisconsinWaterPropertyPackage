function x = QualityFromDensity(rho,rhol,rhog)
    v  = 1./rho;
    vl = 1./rhol;
    vg = 1./rhog;
    
    x = (v - vl)./(vg-vl);
end