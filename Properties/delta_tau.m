function del_t = delta_tau(delta,tau)
    % Isobaric derivative of delta w.r.t. tau
    
    P_delta = PressureOneR_delta(delta,tau);    
    P_tau   = PressureOneR_tau  (delta,tau);
    
    del_t = -P_tau ./ P_delta;
    
end