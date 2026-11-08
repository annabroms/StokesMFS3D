function zeta = stagnant_cap_exact_drag(phi,muIn,muOut,radius)
%STAGNANT_CAP_EXACT_DRAG Sadhal-Johnson axial resistance, note equation (50).
if nargin == 0
    self_test();
    return
end
assert(all(phi>=0 & phi<=pi,'all'));
assert(muIn>=0 && muOut>0 && radius>0);
ratio = muIn/muOut;
angleTerm = 2*phi+sin(phi)-sin(2*phi)-sin(3*phi)/3;
zeta = 4*pi*muOut*radius*((2+3*ratio)/(2*(1+ratio)) + ...
    angleTerm/(4*pi*(1+ratio)));
end

function self_test()
%SELF_TEST Verify clean-drop, solid-sphere, and clean-bubble limits.
mu = 1.7; a = 2.3; ratio = 4;
assert(abs(stagnant_cap_exact_drag(0,ratio*mu,mu,a)- ...
    2*pi*mu*a*(2+3*ratio)/(1+ratio))<1e-12);
assert(abs(stagnant_cap_exact_drag(pi,ratio*mu,mu,a)-6*pi*mu*a)<1e-12);
assert(abs(stagnant_cap_exact_drag(0,0,mu,a)-4*pi*mu*a)<1e-12);
angle = pi/3;
A = 2*angle+sin(angle)-sin(2*angle)-sin(3*angle)/3;
assert(abs(stagnant_cap_exact_drag(angle,0,mu,a)- ...
    mu*a*(4*pi+A))<1e-12);
fprintf('stagnant_cap_exact_drag self-test passed.\n');
end
