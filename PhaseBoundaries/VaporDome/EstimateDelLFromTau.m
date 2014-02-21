function delL = EstimateDelLFromTau(tau)
%   Approximate saturated liquid density given a saturation temperature.
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
    
    delL =      1.99274064E+0 * VarTheta.^(  1/3) ...
            +   1.09965342E+0 * VarTheta.^(  2/3) ...
            -   5.10839303E-1 * VarTheta.^(  5/3) ...
            -   1.75493479E+0 * VarTheta.^( 16/3) ...
            -   4.55170352E+1 * VarTheta.^( 43/3) ...
            -   6.74694450E+5 * VarTheta.^(110/3) ...
            +   1;
    
end