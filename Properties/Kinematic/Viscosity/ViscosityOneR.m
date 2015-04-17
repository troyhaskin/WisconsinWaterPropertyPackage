function mu = ViscosityOneR(delta,tau)
    
    muRef = DimensioningViscosity()            ; %[Pa-s]
    muND  = ViscosityOneRND(delta,tau)   ; %[-]
    
    mu    = muND * muRef;%[Pa-s]
end