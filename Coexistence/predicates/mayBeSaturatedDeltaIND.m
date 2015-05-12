function mayBeSaturated = mayBeSaturatedDeltaIND(delta,iND)
    
    [delGmin,delLmax] = saturableDeltas()   ;
    delGMd    = +4.8776578879917196E-02     ;
    delL1     = +2.8459913232353582E+00     ;
    delL2     = +3.0782066729942112E+00     ;
    [delLt,~] = TriplePointDensitiesR()     ;

    
    %   First gas region
    delLo = delGmin ;
    delHi = delGMd  ;
    p     = [-3.2499910703192775E-02;-7.4671888646641177E-02;-2.1636053030460331E-02;...
             +9.0483066918309166E-02;+3.9796278206848351E-01;+8.3561702716351043E+00];
    mu   = [-6.6330916880020361E+00;+3.4987820465749464E+00];
    fit  = @(d) HornersMethod(log(d),p,mu);
    checkIfBelow(delLo,delHi,fit);
    
    %   Second gas region
    delLo = delGMd  ;
    delHi = 1.0     ;
    p     = [-7.8653690909833998E-02;+1.2262659029693300E-01;+7.3900753243268225E-02;...
		     -1.0481383180212968E-01;+2.6280313141655121E-02;-8.7408880015500223E-01;...
		     +7.8123098789296312E+00];
    mu    = [+4.7995899158744704E-01;+3.8842425427048866E-01];
    fit   = @(d) HornersMethod(log(d),p,mu);
    checkIfBelow(delLo,delHi,fit);
    
    %   First liquid region
    delLo = 1.0     ;
    delHi = delL1   ;
    p     = [-3.0624027886252880E-03;+9.1976024435002756E-03;-2.4292421224378975E-02;...
		     -1.5269278369899392E-01;-8.1984060638154010E-01;+5.8550232706532022E+00];
    mu    = [+1.5969827410431401E+00;+4.6472809781185104E-01];
    fit   = @(d) HornersMethod(d,p,mu) + 6.1045559161394536E-03;
    checkIfBelow(delLo,delHi,fit);
    
    %   Second liquid region
    delLo = delL1   ;
    delHi = delL2   ;
    p     = [-1.8133000982982282E-03;-6.0537095366927111E-03;-9.8053104622959682E-03;...
		     -3.7152584808103573E-02;-4.3127426628907917E-01;+1.4004968109677314E+00];
    mu    = [+2.9765993371528601E+00;+6.8095682282090120E-02];
    fit   = @(d) HornersMethod(d,p,mu) + 7.4979900206617600E-04;
    checkIfBelow(delLo,delHi,fit);

    %   Third liquid region
    delLo       = delL2   ;
    delHi       = delLmax ;
    [~,iLAtMax] = iNDsAtSaturableDeltas();
    c           = [+5.1770062815506157E-01;+3.2967387240093049E+00;+2.1813809229584362E+00];
    shift       = +5.9565801837926924E-04;
    fit   = @(d) c(2)*(delLmax-d).^c(1).*exp(c(3)*(delLmax-d)) + iLAtMax + shift;
    checkIfBelow(delLo,delHi,fit);
    
    %   Fourth liquid region
    delLo       = delLt     ;
    delHi       = delLmax   ;
    [~,iLAtMax] = iNDsAtSaturableDeltas();
    c           = [+4.9763662117192797E-01;-2.7486285142106492E+00;-2.2862334141468203E+01];
    shift       = -1.1372162114144102E-05;
    fit   = @(d) c(2)*(delLmax-d).^c(1).*exp(c(3)*(delLmax-d)) + iLAtMax + shift;
    checkIfAbove(delLo,delHi,fit);
    
    
    function [] = checkIfBelow(delLo,delHi,fit)
        inRegion = (delLo <= delta) & (delta <= delHi);
        if any(inRegion)
            mayBeSaturated(inRegion) = iND(inRegion) <= fit(delta(inRegion));
        end
    end
    function [] = checkIfAbove(delLo,delHi,fit)
        inRegion = (delLo <= delta) & (delta <= delHi);
        if any(inRegion)
            isAbove                  = fit(delta(inRegion)) <= iND(inRegion);
            mayBeSaturated(inRegion) = mayBeSaturated(inRegion) & isAbove   ;
        end
    end

end