function delL = EstimateDelLFromTau(tau)
    
    b    = [1.99274064,1.09965342,-5.10839303E-1,-1.75493479,-45.5170352,-6.74694450E5];
    
    VarTheta = 1 - 1./tau;
    
    Part = 1;
    Part = Part + b(1) * VarTheta.^(  1/3);
    Part = Part + b(2) * VarTheta.^(  2/3);
    Part = Part + b(3) * VarTheta.^(  5/3);
    Part = Part + b(4) * VarTheta.^( 16/3);
    Part = Part + b(5) * VarTheta.^( 43/3);
    delL = Part + b(6) * VarTheta.^(110/3);
    
end