function [Prop1,Prop2,varargout] = ReshapeProperties(Prop1,Prop2)
    
    
    varargout{1} = size(Prop1);
    varargout{2} = size(Prop2);
    
    Prop1 = Prop1(:);
    Prop2 = Prop2(:);
    
    if (length(Prop1) ~= length(Prop2)) && ...
            not(isscalar(Prop1))         && ...
            not(isscalar(Prop2))
        error('Prop2hermodynamics:UnequalStates',...
            ['Prop2he number of elements for the two specified '...
            'properties must be equal.']);
    end
    
end
