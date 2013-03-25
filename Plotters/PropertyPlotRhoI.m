function varargout = PropertyPlotRhoI(IsoTherms)
    
    [rhoLt,rhoGt] = TriplePointDensities();
    [iLt  ,iGt  ] = TriplePointInternalEnergies();
    x = linspace(0,1,200)';
    rhot = 1./(x./rhoLt + (1-x)./rhoGt);
    it   = x.*iLt   + (1-x).*iGt;
    
    rho  = GeometricSpace(0.95*rhoGt,999.96,1.025,500);
    Tsat = linspace(TriplePointTemperature,CriticalTemperature(),500)';

    [~,rhoL,rhoG] = SaturationStateGivenTsat(Tsat);
    iL = InternalEnergyOne(rhoL,Tsat);
    iG = InternalEnergyOne(rhoG,Tsat);
    
    Gray      = 0.4*[1,1,1];
    LineWidth = 1.4;
    figure();
    semilogx([rhoG;rhoL(end:-1:1)],[iG;iL(end:-1:1)],'Color',Gray,'LineWidth',LineWidth);
    hold('on');
    semilogx(rhot,it,'--','Color',Gray,'LineWidth',LineWidth);
    
    if (nargin < 1) || isempty(IsoTherms)
        Nisotherms = 12;
        Tiso = linspace(1.02*TriplePointTemperature,0.95*CriticalTemperature(),Nisotherms)';
    else
        Tiso = IsoTherms;
    end
    for k = 1:Nisotherms
        i = InternalEnergy(rho,Tiso(k));
        semilogx(rho,i,'-','Color',Gray,'LineWidth',LineWidth);
    end
    
    box('on');
    grid('on');
    axis([2E-3,2E3,-0.2E6,3E6]);
    xlabel('Density [kg/m^3');
    ylabel('Internal Energy [J/kg]');
    
    hold('off');
    
    if nargout >= 1
        varargout{1} = gca;
    end
end





