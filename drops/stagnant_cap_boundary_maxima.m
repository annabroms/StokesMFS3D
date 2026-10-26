function maxima = stagnant_cap_boundary_maxima(report,validation,solution)
%STAGNANT_CAP_BOUNDARY_MAXIMA Four imposed-condition maxima by surface region.
% Order: mobile velocity continuity, mobile no penetration, mobile
% tangential traction jump, cap exterior relative tangential velocity.
% Velocity conditions are divided by U_*, traction by mu_out*U_*/a.
if nargin == 0
    self_test();
    return
end
assert(validation.count==numel(report.point.continuity));
assert(any(validation.mobile) && any(validation.cap));
assert(solution.velocityScale>0 && solution.tractionScale>0);
mobile = validation.mobile;
cap = validation.cap;
maxima = [ ...
    max(report.point.continuity(mobile))/solution.velocityScale, ...
    max(report.point.normal(mobile))/solution.velocityScale, ...
    max(report.point.mobileTraction(mobile))/solution.tractionScale, ...
    max(report.point.capOuter(cap))/solution.velocityScale];
assert(all(isfinite(maxima)) && all(maxima>=0));
end

function self_test()
%SELF_TEST Verify each region, maximum, and physical normalization scale.
validation.cap = [true;true;false;false];
validation.mobile = ~validation.cap;
validation.count = 4;
report.point.continuity = [100;100;2;6];
report.point.normal = [100;100;4;8];
report.point.mobileTraction = [100;100;12;20];
report.point.capOuter = [3;9;100;100];
solution.velocityScale = 2;
solution.tractionScale = 4;
actual = stagnant_cap_boundary_maxima(report,validation,solution);
assert(isequal(actual,[3,4,5,4.5]));
fprintf('stagnant_cap_boundary_maxima self-test passed.\n');
end
