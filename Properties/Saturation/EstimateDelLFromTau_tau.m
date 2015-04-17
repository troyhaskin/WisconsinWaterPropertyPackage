function delLdtau = EstimateDelLFromTau_tau(tau)
    
    b    = [1.99274064,1.09965342,-5.10839303E-1,-1.75493479,-45.5170352,-6.74694450E5];
    
    VarTheta = 1 - 1./tau;
    
    % Derivatives wrt to VarTheta
    Part = 0;
    Part = Part + b(1) * (  1/3) * VarTheta.^( -2/3);
    Part = Part + b(2) * (  2/3) * VarTheta.^( -1/3);
    Part = Part + b(3) * (  5/3) * VarTheta.^(  2/3);
    Part = Part + b(4) * ( 16/3) * VarTheta.^( 13/3);
    Part = Part + b(5) * ( 43/3) * VarTheta.^( 40/3);
    Part = Part + b(6) * (110/3) * VarTheta.^(107/3);
    
    % Multiply by derivative of VarTheta wrt tau
    delLdtau = Part ./ (tau.^2);
    
end