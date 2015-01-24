function [P,rhoL,rhoG] = SaturationStateGivenTsat(T,varargin)
    
    [P,rhoL,rhoG] = SaturationStateGivenTemperature(T,varargin{:});
    
end
