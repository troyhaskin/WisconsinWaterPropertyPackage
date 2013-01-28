function Delta = GetDelta(deltaMod,Theta,B,a)
    Delta = Theta.^2 + B * deltaMod.^a;
end