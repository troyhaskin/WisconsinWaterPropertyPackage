function [delL,delG] = NearTripleStabilityDensityLimitsR()
    %   Near the triple point, the saturation line is a multi-valued function of reduced
    %   density.  As such, this limit is used in algorithms where this may cause issues;
    %   the method of resolving those issues is left to the caller.

    delL = 3.1049457143789; % A shift of 5E-12 below the critical reduce density
    delG = 0.0            ; 
end