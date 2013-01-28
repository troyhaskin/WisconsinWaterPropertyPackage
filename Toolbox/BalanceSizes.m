function [Numeric1, Numeric2] = BalanceSizes(Numeric1,Numeric2,~)
    
%     if (nargin < 3)
%         FillValue = 0;
%     end
    
    Numeric1ScalarNumeric2Array  =     isscalar(Numeric1)  && not(isscalar(Numeric2));
    Numeric1ArrayNumeric2Scalar  = not(isscalar(Numeric1)) &&     isscalar(Numeric2);
%     Numeric1VectorNumeric2Vector =     isvector(Numeric1)  &&     isvector(Numeric2);
    
   
    if Numeric1ScalarNumeric2Array 
        % Expand scalar Numeric1 to the size of array Numeric2
        Numeric1 = ones(size(Numeric2)) * Numeric1;
    
    elseif Numeric1ArrayNumeric2Scalar 
        % Expand scalar Numeric2 to the size of array Numeric1
        Numeric2 = ones(size(Numeric1)) * Numeric2;
        
%     elseif Numeric1VectorNumeric2Vector  
%         % Expand the shorter vector to the size the other with padding equal to FillValue
%         SizeNumeric1 = size(Numeric1);
%         SizeNumeric2 = size(Numeric2);
%         
%         if any(SizeNumeric1 > SizeNumeric2)
%             SizeDiff                = SizeNumeric1 - SizeNumeric2;
%             SizeDiff(SizeDiff == 0) = 1;
%             Padding                 = repmat(FillValue,SizeDiff);
%             Numeric2                = [];
%             TempVector = zeros(size(Numeric1)) + FillValue;
%             TempVector(1:LengthNumeric2) = Numeric2;
%             Numeric2 = TempVector;
%             
%             
%         else
%         end
        
        
    elseif (not(isscalar(Numeric1)) && not(isscalar(Numeric2)))
        if any(size(Numeric1) ~= size(Numeric2))
            warning('Thermodynamics:BalanceSizes:NotScalarUnequalSize',...
                  'Cannot balance array with different dimensions.');
        end
    end
    
end