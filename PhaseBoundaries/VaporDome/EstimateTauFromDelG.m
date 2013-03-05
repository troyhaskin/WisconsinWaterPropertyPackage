function tau = EstimateTauFromDelG(delG)
    
    %   A collection of brute-force polynomial fits that returns an estimate of
    %   the inverse, reduced temperature.  The fits are made from solving the 
    %   full equilibrium system are are intended as guess values when density
    %   is the provided saturation value for solving the system.
    %
    %   The fits are in fives zones where the divisions were somewhat arbitrary
    %   and mostly made based on qualitatively accuracy requirements:
    %       1.) 273.16K - 350.00K
    %       2.) 350.00K - 500.00K
    %       3.) 500.00K - 625.00K
    %       4.) 625.00K - 640.00K
    %       5.) 640.00K - 647.096K
    %
    
    %
    delGTrip = 1.50763221266E-005;
    delG350K = 8.08349894630E-004;
    delG500K = 4.09903928885E-002;
    delG625K = 3.67361100488E-001;
    delG640K = 5.49628716184E-001;
    delGCrit = 1                 ;
    
    % Setup
    [delG,SizeDelG] = Columnify(delG)           ; % Creates a column vector with a restore size
    Zone1           = (delG >  delGTrip) & (delG <  delG350K);
    Zone2           = (delG >  delG350K) & (delG <  delG500K);
    Zone3           = (delG >  delG500K) & (delG <  delG625K);
    Zone4           = (delG >  delG625K) & (delG <  delG640K);
    Zone5           = (delG >  delG640K) & (delG <  delGCrit);
    
    tau        = delG * 0;
    tau(Zone1) = TauGuessZone1(delG(Zone1));
    tau(Zone2) = TauGuessZone2(delG(Zone2));
    tau(Zone3) = TauGuessZone3(delG(Zone3));
    tau(Zone4) = TauGuessZone4(delG(Zone4));
    tau(Zone5) = TauGuessZone5(delG(Zone5));
    
    tau = RestoreShape(tau,SizeDelG);
end



function tau = TauGuessZone1(delG)
    
    c = [   -2.9589357885E-002;
            +2.1153805660E-001;
            -5.1614870033E-001;
            +3.0702575646E-001;
            +6.5892530382E-001;
            -9.3089278132E-001;
            -1.5616072092E-001;
            +7.3124422721E-001;
            -1.4761970657E-001;
            -1.7572616574E-001;
            +4.2374288601E-003;
            +8.7895657301E-002;
            -1.3732169020E-001;
            +1.9931539137E+000];
    
    mu  = [+2.7611837294E-004,+2.7783040898E-004];
    tau = HornersMethod(delG,c,mu);
end

function tau = TauGuessZone2(delG)
    
    c = [   -2.9635671671E-002;
            +1.9947630593E-001;
            -4.3961595870E-001;
            +1.5061778274E-001;
            +7.0932833252E-001;
            -7.2515906404E-001;
            -3.4527518867E-001;
            +6.5568162217E-001;
            -2.4380445442E-002;
            -1.7540787161E-001;
            -2.2545601455E-002;
            +9.2422377408E-002;
            -1.4221311618E-001;
            +1.4421248612E+000];
    
    mu  = [+1.4474673363E-002,+1.4133645680E-002];
    tau = HornersMethod(delG,c,mu);
end

function tau = TauGuessZone3(delG)
    
    c = [   -8.7297617453E-005;
            +5.8065455012E-004;
            -1.2903780773E-003;
            +6.5567462711E-004;
            +1.1781358270E-003;
            -4.5048428211E-004;
            -2.9265209845E-003;
            +4.3610908150E-003;
            -4.7170364894E-003;
            +8.4553733391E-003;
            -1.7246666802E-002;
            +3.6638584326E-002;
            -8.0872318273E-002;
            +1.1168202303E+000];
    
    mu  = [+1.5910152412E-001,+1.1139982970E-001];
    tau = HornersMethod(delG,c,mu);
end

function tau = TauGuessZone4(delG)
    
    c = [   -4.0683139662E-011;
            +2.9937814049E-010;
            +4.4486162761E-010;
            -7.5537462272E-010;
            +7.1290081629E-009;
            +4.8765813026E-008;
            +5.7627941197E-008;
            +1.3794800859E-007;
            -2.5322921359E-006;
            +1.4032304054E-005;
            -1.4255648424E-004;
            +1.3759996416E-003;
            -8.6217895680E-003;
            +1.0217845675E+000];
    
    mu  = [+4.4866559221E-001,+6.3776869503E-002];
    tau = HornersMethod(delG,c,mu);
end

function tau = TauGuessZone5(delG)
    
    c = [   +2.7669679037E-010;
            -2.9697997784E-009;
            +1.1723760516E-008;
            -1.7308137458E-008;
            -9.9055403256E-009;
            +6.6427937407E-008;
            -1.5542094136E-007;
            +7.0842654039E-007;
            -5.1202509187E-006;
            +3.1538225930E-005;
            -3.0709676426E-004;
            +1.7924322721E-003;
            -4.4720435192E-003;
            +1.0039251895E+000];
    
    mu  = [+6.7709296663E-001,+1.1948666743E-001];
    tau = HornersMethod(delG,c,mu);
end
