function [A,B,C,D,E,F,G,H,K] = Invert3by3(varargin)
    
    if (length(varargin) == 1)
        CellMatrix = num2cell(varargin{1});
        [a,b,c,d,e,f,g,h,k] = CellMatrix{:};
        
    elseif length(varargin) == 9
        [a,b,c,d,e,f,g,h,k] = varargin{:};
    else
        error('Thermodynamics:Toolbox:WrongInputCount',...
              'Invert3by3 takes 1 matrix or 9 vector element-wise arguments only.');
    end
       
    
    iDet = 1./(a.*(e.*k- f.*h) + b .*(f.*g- d.*k) + c.*(d.*h- e.*g));
    
    A = (e.*k- f.*h).*iDet;
    B = (f.*g- d.*k).*iDet;
    C = (d.*h- e.*g).*iDet;
    
    D = (c.*h- b.*k).*iDet;
    E = (a.*k- c.*g).*iDet;
    F = (g.*b- a.*h).*iDet;
    
    G = (b.*f- c.*e).*iDet;
    H = (c.*d- a.*f).*iDet;
    K = (a.*e- b.*d).*iDet;
    
    if (nargout == 1) || (nargout == 0)
        A = [A , B , C  ;...
             D , E , F  ;...
             G , H , K  ];
    end
    
end