function Handle = PhaseBoudariesPT(Npoints)
    
    if (nargin < 1) || isempty(Npoints)
        Npoints = 500;
    end
    
    Tt = TriplePointTemperature();
    Tc = CriticalTemperature();
    
    Tih  = GeometricSpace(251.165,Tt,0.960, Npoints)';
    Tiii = linspace(251.165 , 256.164 , Npoints);
    Tv   = linspace(256.164 , 273.31  , Npoints);
    Tvi  = linspace(273.31  , 355     , Npoints);
    Tvii = linspace(355     , 715     , Npoints);
    Tsub = GeometricSpace(50      , Tt ,0.98     , Npoints);
    Tsat = linspace(Tt      , Tc      , Npoints);
    
	Pih  = PressureMeltIh   (Tih)   ;
    Piii = PressureMeltIII  (Tiii)  ;
    Pv   = PressureMeltV    (Tv)    ;
    Pvi  = PressureMeltVI   (Tvi)   ;
    Pvii = PressureMeltVII  (Tvii)  ;
    Psub = PressureSublimate(Tsub)  ;
    
    Tmelt = [Tih,Tiii,Tv,Tvi,Tvii];
    Pmelt = [Pih,Piii,Pv,Pvi,Pvii];
    
    [Psat,~,~] = SaturationStateGivenTsat(Tsat);
    
    
    Handle = semilogy(Tmelt,Pmelt/1E6,'b',Tsub,Psub/1E6,'r',Tsat,Psat/1E6,'k','LineWidth',2);
    xlabel('Temperature [K]');
    ylabel('Pressure [MPa]');
    axis([200,700,1E-7,5E4]);
    
end