function Pnd_tauSat = ClausiusClapeyronRRND(Pnd,tauSat,delL,delG,iNDL,iNDG)

    % Latent internal energy - if not passed, it will be calculated.
    if (nargin < 5)
        iNDL  = InternalEnergyOneRND(delL,tauSat);
        iNDG  = InternalEnergyOneRND(delG,tauSat);
    end
    iNDLG = iNDG - iNDL                        ;

    % Clausius–Clapeyron relation
    Pnd_tauSat = ( iNDLG ./ (1./delL - 1./delG) - Pnd) ./ tauSat;

    %   Standard form of the Clausius–Clapeyron relation:
    %       P_Tsat = hlg / (Tsat * vlg)
    %
    %   This form was not used because it is not dimensionless, and the form 
    %   above uses internal energies since it is less computationally intensive.

end