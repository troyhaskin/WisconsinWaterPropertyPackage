function H = ViscosityDiluteGasLimitConstants()
    % The constants are assigned in the documented order
    % but are stored in reverse to be compatible with
    % Horner's Method in the call.
    H(4) = +1.677520E+0;
    H(3) = +2.204620E+0;
    H(2) = +6.366564E-1;
    H(1) = -2.416050E-1;
end