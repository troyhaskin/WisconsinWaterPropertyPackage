function delG = DoNoIterationValueReducedDensityGas()
    %   Any density above this value will use the curve fits and not 
    %   iteration for the saturation properties due to an inability to 
    %   get the system to robustly converge above the value;
    delG = 0.54962871618379; 
end