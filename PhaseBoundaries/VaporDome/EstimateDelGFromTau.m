function delG = EstimateDelGFromTau(tau)
%   Approximate saturated gas density given a saturation temperature.
%   This is a very primitive function so no input checking (such as ensuring the tau
%   is within its saturation limits) is performed.
% 
%   Source:
%       Wagner, Wolfgang, and A. Pruss. "International equations for the saturation 
%       properties of ordinary water substance. Revised according to the international 
%       temperature scale of 1990. Addendum to J. Phys. Chem. Ref. Data 16, 893 (1987)." 
%       Journal of Physical and Chemical Reference Data 22.3 (1993): 783-787.
%

    VarTheta = 1 - 1./tau;
    
    delG = exp(...
                -2.03150240E+0 * VarTheta.^( 1/3)  ...
                -2.68302940E+0 * VarTheta.^( 2/3)  ...
                -5.38626492E+0 * VarTheta.^( 4/3)  ...
                -1.72991605E+1 * VarTheta.^( 3  )  ...
                -4.47586581E+1 * VarTheta.^(37/6)  ...
                -6.39201063E+1 * VarTheta.^(71/6)  ...
            );
 
end