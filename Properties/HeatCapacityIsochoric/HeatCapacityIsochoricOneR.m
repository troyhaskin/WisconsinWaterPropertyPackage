function cv = HeatCapacityIsochoricOneR(delta,tau)
    
    cvStar = DimensioningHeatCacvStar();
    cvND   = HeatCapacityIsochoricOneRND(delta,tau)   ;
    
    cv = cvND * cvStar;
    
end