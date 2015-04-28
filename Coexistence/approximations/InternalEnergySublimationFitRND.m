function iND = InternalEnergySublimationFitRND(delta)
    
    %   Polynomial coefficients
    p = [...
        +1.0432097155419256E-04;
		+3.2280354066173116E-03;
		+2.0525126191073412E-02;
		+1.4691227083450359E-01;
		+7.7448830797229196E+00];
    mu = [...
        -1.5343944820087060E+01
		+3.5952346109503961E+00];
    
    %   Evaluate
    iND = HornersMethod(log(delta),p,mu);
    
end