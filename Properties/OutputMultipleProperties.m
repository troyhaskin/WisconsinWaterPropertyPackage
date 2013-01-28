function CellOutput = OutputMultipleProperties(Psi,SinglePropertySize)

    Nprops     = size(Psi,2);
    CellOutput = cell(1,2);
    
    for k = 1:Nprops
        CellOutput{k} = RestoreShape(Psi(:,k),SinglePropertySize);
    end
    
end