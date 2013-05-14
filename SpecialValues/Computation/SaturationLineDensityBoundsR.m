function [delL,delG] = SaturationLineDensityBoundsR()
    %   Reduced densities outside the internal defined by these values cannot lie
    %   on the saturation line at all.

    delL     = MaximumSaturationDensityR();
    [~,delG] = TriplePointDensitiesR();
end