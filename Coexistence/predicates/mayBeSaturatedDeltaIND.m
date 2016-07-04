function mayBeSaturated = mayBeSaturatedDeltaIND(delta,iND)
    
    %   Set-up constants
    [delGmin,delLmax] = saturableDeltas()   ;
    [delLt,~] = TriplePointDensitiesR()     ;

    %   Allocate
    mayBeSaturated = (delGmin <= delta) & (delta <= delLmax);


    %   First Region
    delLo = delGmin ;
    delHi = +8.3150922949780570E-01;
    p     = [...
        -6.7383117203846383E-04;-3.7822522555875965E-04;+1.7245159813175481E-02;...
        +3.3644952181128424E-02;-7.3764268933753735E-02;-1.6183465099368938E-01;...
        +2.4524989269529882E-01;+3.7075792356302606E-01;-7.9978201755559697E-01;...
        -3.5398835980831334E-01;+8.6809450709023608E+00];
    mu   = [+1.1425287830295507E+00;+1.2883257928093858E-01];
    in   = (delLo <= delta) & (delta <= delHi);
    iFit = HornersMethod((-log(delta(in))).^(1/10),p,mu);
    mayBeSaturated(in) = iND(in) <= iFit;
    
    %   Second Region
    delLo = delHi  ;
    delHi = +1.7200399120550229E+00;
    p     = [...
        +7.9749237168808234E-04;-6.8470131903780108E-04;-6.6350801762660797E-03;...
        +4.9258616082293611E-03;+1.2116996586194539E-02;-4.7949738871147354E-02;...
        -8.4080465340015939E-02;-4.3708634957394443E-01;+6.4162827595504375E+00];
    mu    = [+2.1175754244628892E-01;+2.5535936865002246E-01];
    in   = (delLo <= delta) & (delta <= delHi);
    iFit = HornersMethod(log(delta(in)),p,mu);
    mayBeSaturated(in) = iND(in) <= iFit;
    
    %   Third Region
    delLo = delHi   ;
    delHi = delLmax ;
    p     = [...
        -1.7411654429974298E-04;-1.5302447049417923E-03;-2.2278202903091180E-03;...
        -7.2701301421645593E-03;+2.1904606117070133E-01;+1.6709907582231713E+00;...
        +2.2219368629781862E+00];
    mu    = [+1.2870040788074313E+00;+7.6116441164427562E-01];
    in    = (delLo <= delta) & (delta <= delHi);
    iFit  = HornersMethod(sqrt(delLmax^2- delta(in).^2),p,mu);
    mayBeSaturated(in) = iND(in) <= iFit;
    
    %   Fourth Region
    delLo = delLt   ;
    delHi = delLmax ;
    p     = [+9.7084541797564503E-05;-1.6302461315719365E-02;+2.8030350009359235E-02];
    mu    = [+2.5131829335951629E-02;+1.4643342753719566E-02];
    in    = (delLo <= delta) & (delta <= delHi);
    iFit  = HornersMethod(sqrt(delLmax^2- delta(in).^2),p,mu);
    mayBeSaturated(in) = iND(in) <= iFit;

end