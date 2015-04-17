function cv = HeatCapacityIsochoricOneR(delta,tau)
    
    cvStar = DimensioningHeatCapacity();
    cvND   = HeatCapacityIsochoricOneRND(delta,tau)   ;
    
    cv = cvND * cvStar;
    
end