function varargout = VectorChunk(Vector,ChunkSize)
    
    if (nargin ~= 2) || isempty(ChunkSize)
        varargout{1} = Vector;
        return;
    else
        varargout = cell(1,nargout);
    end
    
    [Nrow,Ncol] = size(Vector);
    
    if Nrow > 1
        Nvector = Nrow;
        Empty   = zeros(0,1);
    else
        Nvector = Ncol;
        Empty   = zeros(1,0);
    end
    
    % Scalar case
    if isscalar(ChunkSize)
        Bottom = 1;
        for k = 1:nargout
            
            if (Bottom <= Nvector) && (ChunkSize ~= 0)
                Top          = min(Bottom + ChunkSize - 1,Nvector);
                varargout{k} = Vector(Bottom:Top);
                Bottom       = Top + 1;
            else
                varargout{k} = Empty;
            end
            
        end
        
        
    else
        
        Nchunks = length(ChunkSize);
        Bottom  = 1;
        
        % Vector form
        for k = 1:Nchunks
            
            if (Bottom <= Nvector) && (ChunkSize(k) ~= 0)
                Top          = min(Bottom + ChunkSize(k) - 1,Nvector);
                varargout{k} = Vector(Bottom:Top);
                Bottom       = Top + 1;
            else
                varargout{k} = Empty;
            end
            
        end
        
        % Scalar continuation for remaining outputs
        if Nchunks < nargout
            for k = Nchunks+1:nargout
                
                if (Bottom <= Nvector) && (ChunkSize(Nchunks) ~= 0)
                    Top          = min(Bottom + ChunkSize(Nchunks) - 1,Nvector);
                    varargout{k} = Vector(Bottom:Top);
                    Bottom       = Top + 1;
                else
                    varargout{k} = Empty;
                end

            end
        end
        
    end
    
end
