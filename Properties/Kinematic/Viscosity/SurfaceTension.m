function sigma = SurfaceTension(T)
    tau = 1 - T/CriticalTemperature()   ;
    B     = 235.8E-3                    ; %[N/m]
    b     = -0.625                      ; %[-]
    mu    = 1.256                       ; %[-]
    sigma = B * tau.^mu .* (1 + b*tau)  ;
end