function varargout = PropertyReshape(varargin)
    
    switch(length(varargin))
        case(1)
            % Two independent properties passed:
            %      - return original size and columnfied versions.
            
            % Pull independent properties
            Prop = varargin{:};
            
            % Outputs
            varargout{1} = Prop(:)      ;
            varargout{2} = size(Prop)   ;
        
        case(2)
            % Two independent properties passed:
            %      - return original size and columnfied versions.
            
            % Pull independent properties
            [Prop1,Prop2] = varargin{:};
            
            % Due a state count (element) check.
            if (length(Prop1) ~= length(Prop2)) && ...
                    not(isscalar(Prop1))         && ...
                    not(isscalar(Prop2))
                error('Thermodynamics:UnequalStates'                    ,...
                    ['The number of elements for the two specified '    ,...
                    'properties must be equal or one must be a '       ,...
                    'scalar value.']);
            end
            
            % Outputs
            varargout{1} = Prop1(:)     ;
            varargout{2} = Prop2(:)     ;
            varargout{3} = size(Prop1)  ;
            varargout{4} = size(Prop2)  ;
            
        case(3)
            % One dependent property and two sizes passed:
            %      - return reshaped dependent property to the 'largest' size passed
            
            % Pull the inputs
            [Psi,SizeProp1,SizeProp2] = varargin{:};
            
            % Reshape the property
            switch(sum(SizeProp1) >= sum(SizeProp2))
                case true
                    varargout{1} = reshape(Psi,SizeProp1);
                case false
                    varargout{1} = reshape(Psi,SizeProp2);
            end
            
        otherwise
            error('Thermodynamics:WrongNumberOfInputs'          ,...
                  'PropertyReshape takes 2 or 3 inputs only.')  ;
    end
    
    
end
