function [rhoLo,rhoHi] = saturableDensities()
    [delLo,delHi] = saturableDeltas();
    rhoLo         = delLo * CriticalDensity();
    rhoHi         = delHi * CriticalDensity();
end