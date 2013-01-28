function MaskedVec = SmartMask(Vector,Mask)
    
    if isscalar(Vector)
        MaskedVec = Vector;
        
    elseif length(Vector) == length(Mask)
        MaskedVec = Vector(Mask);
        
    else
        MaskedVec = Vector;
    end
    
end