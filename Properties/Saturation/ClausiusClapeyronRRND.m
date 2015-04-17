function Pnd_tauSat = ClausiusClapeyronRRND(Pnd,tauSat,delL,delG,iNDL,iNDG)

    % Latent internal energy - if not passed, it will be calculated.
    if (nargin < 5)
        iNDL  = InternalEnergyOneRND(delL,tauSat);
        iNDG  = InternalEnergyOneRND(delG,tauSat);
    end
    
    % Latent internal energy and volume difference
    iNDLG = iNDG - iNDL         ;
    vLG   = (1./delG - 1./delL) ;

    % Clausius–Clapeyron relation
    Pnd_tauSat = -(iNDLG./vLG + Pnd) ./ tauSat;

    %   Standard form of the Clausius–Clapeyron relation:
    %       P_Tsat = hlg / (Tsat * vlg)
    %
    %   This form was not used because it is not dimensionlessand uses internal energies
    %   since it is less computationally intensive.

end