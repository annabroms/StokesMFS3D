function lambda0 = getLambda0(F,T,Kin)
    D = Kin'*Kin;
    ab = D\[F;T];
    lambda0 = Kin*ab;
end