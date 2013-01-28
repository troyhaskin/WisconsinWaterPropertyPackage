function Pnd_tauSat = ClausiusClapeyronRRND(Psat,tauSat,delL,delG,il,ig)
    
    [R,rhoc,Tc] = Nondimensionalizers();
    
    % Latent internal energy - if not passed, it will be calculated.
    if (nargin < 5)
        il  = InternalEnergyOneR(delL,tauSat);
        ig  = InternalEnergyOneR(delG,tauSat);
    end
    ilg = ig - il                        ;
    
    % Reduced density change
    nulg = 1./delG - 1./delL;

    % Clausius–Clapeyron relation
    Pnd_tauSat = -(ilg./nulg + Psat/rhoc) ./ (R * Tc * tauSat);    
    
    %   Standard form of the Clausius–Clapeyron relation:
    %       P_Tsat = hlg / (Tsat * vlg)
    %
    %   This form was not used because it is not dimensionless, and the form 
    %   above uses internal energies since it is less computationally intensive.
    
end