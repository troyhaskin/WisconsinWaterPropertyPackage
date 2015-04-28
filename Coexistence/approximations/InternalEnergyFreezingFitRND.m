function iND = InternalEnergyFreezingFitRND(delta)
    
    %   Polynomial coefficients
    p = [...
        +9.6406567149800798E-04
		-1.9014982636480854E-02
		-1.2643465497582723E-01
		-1.1658552299207586E-01];
    mu = [...
        +3.2517510533271032E+00
        +1.3080586467684155E-01];
    
    %   Evaluate
    iND = HornersMethod(delta,p,mu);

end