function tau = DoNoIterationTau()
    %   Any temperature above this value will use the curve fits and not 
    %   iteration for the saturation properties due to an inability to 
    %   get the system to robustly converge above the value;
    tau = 0.9;%1.0110875;
end
