function delL = EstimateDelLFromTau(tau)
    
    VarTheta = 1 - 1./tau;
    
    Part = 1;
    Part = Part +  1.99274064    * VarTheta.^(  1/3);
    Part = Part +  1.09965342    * VarTheta.^(  2/3);
    Part = Part -  5.10839303E-1 * VarTheta.^(  5/3);
    Part = Part -  1.75493479    * VarTheta.^( 16/3);
    Part = Part - 45.5170352     * VarTheta.^( 43/3);
    delL = Part -  6.74694450E5  * VarTheta.^(110/3);
    
end