function [P,T] = SublimationLine(varargin)
    
    %   Call low-level function
    [Pnd,tau] = SublimationLineRND(varargin{:});
    
    %   Dimensionalize
    T = CriticalTemperature()./tau      ;
    P = Pnd * DimensioningPressure()    ;

end