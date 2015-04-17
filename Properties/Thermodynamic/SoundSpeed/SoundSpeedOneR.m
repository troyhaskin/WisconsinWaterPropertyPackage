function w = SoundSpeedOneR(delta,tau)
    
    wStar = DimensioningSoundSpeed();
    wND   = SoundSpeedOneRND(delta,tau);
    
    w = wND * wStar;
    
end