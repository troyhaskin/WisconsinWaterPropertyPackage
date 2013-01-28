function AugmentedVec = AugmentVector(Vector,Augmentor)
   
    if isscalar(Vector)
        AugmentedVec = Vector;
        
    elseif length(Vector) == length(Augmentor)
        AugmentedVec = Vector(Augmentor);
        
    elseif length(Vector) < length(Augmentor)
        Temp         = Vector;
        AugmentedVec = zeros(size(Augmentor));
        
        AugmentedVec(Augmentor) = Temp;
        
    else
        error('Input vector bigger than the augmentation matrix.');
    end
    
end