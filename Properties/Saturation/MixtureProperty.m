function phiMix = MixtureProperty(weight,phiL,phiG)
    phiMix = (1-weight).*phiL + weight.*phiG;
end