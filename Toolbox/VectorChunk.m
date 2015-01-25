function varargout = VectorChunk(Vector,ChunkSize)
    
    if (nargin == 1) || isempty(ChunkSize)
        varargout{1} = Vector;
        return;
    elseif not(isvector(Vector))
        varargout{1} = Vector;
        return;
    end
    
    
    
    %   Determine the type of vector input (row/column) such that an appropriately shaped
    %   empty vector can be created which will not (hopefully) interfere with concatenation.
    [Nrow,Ncol] = size(Vector);
    
    if (Ncol == 1)
        Nvector = Nrow;
        Empty   = zeros(0,1);
        
    else
        Nvector = Ncol;
        Empty   = zeros(1,0);
        
    end
    
    
    
    %   Set-up
    ChunkSizeLength = length(ChunkSize) ;   %   Number of chunks (subvectors of Vector)
    
    if (ChunkSizeLength > 1)
        
        %   Cell array for the chunks
        if (nargout == ChunkSizeLength)
            % Primary use case
            Chunks    = cell(1,nargout) ; % Storage
            Nchunks   = nargout         ; % Loop limit
            NoPacking = true()          ; % Remainder to pack in a cell array
            
            
        elseif (nargout >  ChunkSizeLength)
            %   More outputs than given chunks.
            
            %   Check if scalar continuation is needed
            if (Nvector <= sum(ChunkSize))
                % No.
                Chunks    = cell(1,ChunkSizeLength) ; % Storage
                Nchunks   = ChunkSizeLength         ; % Loop limit
                NoPacking = true()                  ; % Remainder to pack in a cell array
                
                
            else % (Nvector >  sum(ChunkSize))
                % Yes.
                Nremaining = Nvector - sum(ChunkSize(1:end-1));
                Nrequired  = (ChunkSizeLength - 1) + ceil(Nremaining/ChunkSize(end));
                
                if (Nrequired <= nargout)
                    %	With scalar continuation, there are enough outputs for the subvectors
                    Chunks    = cell(1,Nrequired)   ; % Storage
                    Nchunks   = Nrequired           ; % Loop limit
                    NoPacking = true()              ; % Remainder to pack in a cell array

                else % (Nrequired >  nargout)
                    %	With scalar continuation, there are not enough outputs for the subvectors
                    Chunks    = cell(1,Nrequired)   ; % Storage
                    Nchunks   = Nrequired           ; % Loop limit
                    NoPacking = false()             ; % Remainder to pack in a cell array
                end
                
            end


        else % (nargout <  ChunkSizeLength)
            Chunks    = cell(1,ChunkSizeLength) ; % Storage
            Nchunks   = ChunkSizeLength         ; % Loop limit
            NoPacking = false()                 ; % Remainder to pack in a cell array
        end

        ChunkSize = @(k) ChunkSize(k)   ;   % Not exactly needed here; more for scalar case

    else
        Chunks    = cell(1,nargout) ; % Storage
        Nchunks   = nargout         ; % Loop limit
        NoPacking = true()          ; % Remainder to pack in a cell array
        ChunkSize = @(k) ChunkSize  ; % For loop indexing
    end



    Bottom = 1;
    for k = 1:Nchunks

        if (Bottom <= Nvector) && (ChunkSize(k) ~= 0)
            Top       = min(Bottom + ChunkSize(k) - 1,Nvector); % Top must not exceed vector's length
            Chunks{k} = Vector(Bottom:Top)                    ; % Assign to varargout
            Bottom    = Top + 1                               ; % Update lower limit
        else
            Chunks{k} = Empty ;  %
        end

    end
    
    
    varargout = cell(1,nargout);
    if NoPacking
        varargout(1:Nchunks) = Chunks;
    else
        varargout(1:nargout-1) = Chunks(1:nargout-1);
        varargout(end)         = {Chunks(nargout:end)};
    end
 
end
