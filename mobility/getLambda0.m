function [lambda0,A] = getLambda0(F,T,Kin)
    D = Kin'*Kin;
    A = D\[F;T]; % rbm corresponding to the force and torque.
    lambda0 = Kin*A;
end