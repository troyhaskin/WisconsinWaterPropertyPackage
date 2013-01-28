function [delLt,delGt] = TriplePointReducedDensities()

    rhoc          = CriticalDensity()       ;
    [rholt,rhogt] = TriplePointDensities()  ;

    delLt = rholt / rhoc;
    delGt = rhogt / rhoc;

end