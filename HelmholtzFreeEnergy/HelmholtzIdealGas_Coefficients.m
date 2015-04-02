function [n_o,gamma_o] = HelmholtzIdealGas_Coefficients()

    %   These commented coefficients are the exact values from the IAPWS-95 
    %   reference.  The coefficients actually used have been adjusted to 
    %   yield exact 0s for the saturated liquid internal energy and entropy at
    %   the triple point.
    %
    %     n_o_Exact(1) = -8.3204464837497;
    %     n_o_Exact(2) =  6.6832105275932;
    %
    n_o(1) = -8.3204464837420350; % = n_o_Exact(1) + 7.6650e-12
    n_o(2) =  6.6832105275899885; % = n_o_Exact(2) - 3.2117e-12
    n_o(3) =  3.00632          ;
    n_o(4) =  0.012436         ;
    n_o(5) =  0.97315          ;
    n_o(6) =  1.27950          ;
    n_o(7) =  0.96956          ;
    n_o(8) =  0.24873          ;

    gamma_o(1) =  0.0       ;
    gamma_o(2) =  0.0       ;
    gamma_o(3) =  0.0       ;
    gamma_o(4) =  1.28728967;
    gamma_o(5) =  3.53734222;
    gamma_o(6) =  7.74073708;
    gamma_o(7) =  9.24437796;
    gamma_o(8) = 27.50751050;

end