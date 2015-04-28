function [P,T] = MeltingLine(varargin)
    
    [Pnd,tau] = MeltingLineRND(varargin{:});
    
    P = Pnd * DimensioningPressure();
    T = CriticalTemperature() ./ tau;
    
end