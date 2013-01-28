function Pstar = DimensioningPressure()
    
    [R,rhoc,Tc] = Nondimensionalizers();
    
    Pstar = R * rhoc * Tc;
    
end