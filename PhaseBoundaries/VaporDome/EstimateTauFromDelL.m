function tau = EstimateTauFromDelL(delL)
    
    %   A collection of brute-force polynomial fits that returns an estimate of
    %   the inverse, reduced temperature.  The fits are made from solving the 
    %   full equilibrium system and are intended as guess values when density
    %   is the provided saturation value for solving the system.
    %
    %   The fits are in four zones:
    %     1.)   Zone one consists of the low-temperature region of the  
    %           coexistence curve where the water expands with decreasing 
    %           temperature (approximately 273.16K to 281.3K).  This zone is 
    %           multi-valued and any fit for it (one is commented out) is 
    %           replaced by a simple average of the tau at maximum density to 
    %           the triple point tau.
    %     2.)   Zone 2 is from 281.3K to 450K, where the break is made
    %           for more accuracy and lower oscillations as the curve begins
    %           to turn over.
    %     3.)   Zone 3 runs from 450K to 644K, where the break is made such that
    %           high accuracy near the critical point can be made.
    %     4.)   Zone 4 runs from 644K to 1.00001*647.096K and offers a good approximation
    %           of tau almost near the critical point.
    %     5.)   Zone 5 run from 1.00001*647.096K to 647.096 and offers a good approximation
    %           of tau near the critical point.
    
    
    % Optional argument
    delLTrip = 3.10494571438;
%     delLMax  = 3.10535755003;
    delL450K = 2.76503493715;
    delL644K = 1.36693897706;
    delLNear = 1.00001      ;
    delLCrit = 1            ;
    
    % Setup
    [delL,SizeDelL] = Columnify(delL)           ; % Creates a column vector with a restore size
    Zone1           =                      (delL >= delLTrip);
    Zone2           = (delL <  delLTrip) & (delL >  delL450K);
    Zone3           = (delL <  delL450K) & (delL >  delL644K);
    Zone4           = (delL <  delL644K) & (delL >= delLNear);
    Zone5           = (delL <  delLNear) & (delL >= delLCrit);
    
    tau = delL * 0;
    tau(Zone1) = TauGuessZone1(delL(Zone1)); % Constat valuedue to the multi-valuedness of this region
    tau(Zone2) = TauGuessZone2(delL(Zone2));
    tau(Zone3) = TauGuessZone3(delL(Zone3));
    tau(Zone4) = TauGuessZone4(delL(Zone4));
    tau(Zone5) = TauGuessZone5(delL(Zone5));
    
    tau = RestoreShape(tau,SizeDelL);
end

function tau = TauGuessZone1(delL)
    
%     c = [   -4.0202980877E-003;
%             -1.6886050072E-002;
%             -1.4246617901E-002;
%             +2.4000490254E-002;
%             +3.4811519678E-002;
%             -8.6313961917E-003;
%             -2.3341337800E-002;
%             -1.9179798468E-003;
%             +2.0145784116E-003;
%             -9.8818274228E-003;
%             +2.3555504613E+000];
%     
%     mu  = [  +3.1052042484E+000,+1.4988549232E-004];
%     tau = HornersMethod(delL,c,mu);
    
    tau = delL * 0 + 2.3346;
end

function tau = TauGuessZone2(delL)
    
    c = [	+1.0647328624E-001;
            +4.0987236559E-001;
            +1.1131508711E-001;
            -1.2313131781E+000;
            -1.0516339232E+000;
            +1.3394245906E+000;
            +1.5238087266E+000;
            -6.1186168125E-001;
            -8.5056415875E-001;
            +1.2858258845E-001;
            +2.2309047729E-001;
            +6.8522032462E-002;
            +2.3129084425E-001;
            +1.7105303521E+000];
    
    mu  = [+2.9645711505E+000,+1.2368112754E-001];
    tau = HornersMethod(delL,c,mu);
end

function tau = TauGuessZone3(delL)
    
    c = [   +7.0753113463E-005;
            +2.5642740773E-004;
            -8.1434474675E-005;
            -9.9787445842E-004;
            -1.3259216684E-004;
            +1.7413858285E-003;
            +7.8162890255E-004;
            -1.8058518000E-004;
            +1.9855512269E-003;
            +6.5728590948E-003;
            +1.8125076406E-002;
            +6.1727750366E-002;
            +1.4865974225E-001;
            +1.1396478920E+000];
    
    mu  = [+2.2450124894E+000,+4.5327613677E-001];
    tau = HornersMethod(delL,c,mu);
end

function tau = TauGuessZone4(delL)
    
    c = [   +7.3707960642E-015;
            +4.4992448994E-014;
            +1.9031955602E-013;
            -5.8425488682E-013;
            +7.6284699964E-012;
            -1.0467918199E-010;
            +2.0489610065E-009;
            -2.8293355832E-008;
            +5.4103722097E-007;
            -4.1720684046E-006;
            +7.8073138922E-005;
            +8.8625434547E-004;
            +2.5831142141E-003;
            +1.0023628234E+000];
    
    mu  = [+1.2846636314E+000,+1.1053989397E-001];
    tau = HornersMethod(delL,c,mu);
end

function tau = TauGuessZone5(delL)
    c = [   -1.2427888617613720E-13;
            +1.6349800613368345E-11;
            -1.6826603201245960E-09;
            +1.7914607206586920E-07;
            +1.5044848361890351E-06;
            +4.0702911932209303E-06;
            +1.0000036400293815E+00 ];
    mu = [+1.0309140731308570E+00,+1.1619681188514367E-02];
    
    tau = HornersMethod(delL,c,mu);
end