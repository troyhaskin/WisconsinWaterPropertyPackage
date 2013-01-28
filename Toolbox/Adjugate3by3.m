function [A,B,C,D,E,F,G,H,K] = Adjugate3by3(varargin)
    
    
    if (length(varargin) == 1)      % ------------- Single 3by3 matrix passed in
        CellMatrix = num2cell(varargin{1});
        [a,b,c,d,e,f,g,h,k] = CellMatrix{:};
        
    elseif length(varargin) == 9    % ------------- Individual element-wise vectors
        [a,b,c,d,e,f,g,h,k] = varargin{:};
        
    else                            % ------------- *** ERROR ***
        error('Thermodynamics:Toolbox:WrongInputCount',...
              'Invert3by3 takes 1 matrix or 9 vector element-wise arguments only.');
    end
          
    A = (e.*k- f.*h);
    B = (f.*g- d.*k);
    C = (d.*h- e.*g);
    
    D = (c.*h- b.*k);
    E = (a.*k- c.*g);
    F = (g.*b- a.*h);
    
    G = (b.*f- c.*e);
    H = (c.*d- a.*f);
    K = (a.*e- b.*d);
    
    if (nargout == 1) || (nargout == 0)
        A = [A , B , C  ;...
             D , E , F  ;...
             G , H , K  ];
    end
    
end