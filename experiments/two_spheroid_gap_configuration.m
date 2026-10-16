function pair = two_spheroid_gap_configuration(caseName,semiAxes,delta)
%TWO_SPHEROID_GAP_CONFIGURATION Build an exact-gap oriented spheroid pair.
%
% The local major axis is the first coordinate direction.  The two contact
% points are placed at +/-delta/2 along a shared support normal, so the
% Euclidean surface-to-surface distance is exactly delta.  The available
% configurations are "tip", "side", and a non-symmetric "oblique" pair.

caseName = lower(string(caseName));
semiAxes = semiAxes(:);
if ~isscalar(caseName) || ~ismember(caseName,["tip","side","oblique"])
    error('two_spheroid_gap_configuration:BadCase', ...
        'caseName must be "tip", "side", or "oblique".');
end
if numel(semiAxes) ~= 3 || any(~isfinite(semiAxes)) || ...
        any(semiAxes <= 0) || semiAxes(1) <= semiAxes(2) || ...
        semiAxes(2) ~= semiAxes(3)
    error('two_spheroid_gap_configuration:BadSemiaxes', ...
        'semiAxes must define a prolate spheroid [a,b,b] with a>b>0.');
end
if ~isscalar(delta) || ~isreal(delta) || ~isfinite(delta) || delta <= 0
    error('two_spheroid_gap_configuration:BadGap', ...
        'delta must be a positive finite scalar.');
end

degrees = pi/180;
switch caseName
    case "tip"
        rotation1 = eye(3);
        rotation2 = eye(3);
        gapNormal = [1;0;0];
        label = "tip-to-tip";
    case "side"
        rotation1 = eye(3);
        rotation2 = eye(3);
        gapNormal = [0;1;0];
        label = "side-by-side";
    case "oblique"
        rotation1 = rotation_z(18*degrees)*rotation_y(14*degrees);
        rotation2 = rotation_z(101*degrees)*rotation_y(-29*degrees)* ...
            rotation_z(-37*degrees);
        gapNormal = [1;0.62;0.37];
        gapNormal = gapNormal/norm(gapNormal);
        label = "oblique, non-symmetric";
end

assert_rotation(rotation1);
assert_rotation(rotation2);

contactPoint1 = -0.5*delta*gapNormal;
contactPoint2 = 0.5*delta*gapNormal;
supportOffset1 = ellipsoid_support_offset( ...
    semiAxes,rotation1,gapNormal);
supportOffset2 = ellipsoid_support_offset( ...
    semiAxes,rotation2,-gapNormal);
centre1 = contactPoint1-supportOffset1;
centre2 = contactPoint2-supportOffset2;

pair.caseName = caseName;
pair.label = label;
pair.semiAxes = semiAxes;
pair.delta = delta;
pair.gap = norm(contactPoint2-contactPoint1);
pair.gapNormal = gapNormal;
pair.contactPoints = [contactPoint1,contactPoint2];
pair.centres = [centre1,centre2];
pair.rotations = {rotation1,rotation2};

for bodyIndex = 1:2
    localContact = pair.rotations{bodyIndex}.'*( ...
        pair.contactPoints(:,bodyIndex)-pair.centres(:,bodyIndex));
    surfaceLevel = sum((localContact./semiAxes).^2);
    assert(abs(surfaceLevel-1) < 1000*eps, ...
        'Constructed contact point is not on the spheroid surface.');
end
assert(abs(pair.gap-delta) < 100*eps*max(1,delta));
end


function offset = ellipsoid_support_offset(semiAxes,rotation,direction)
%ELLIPSOID_SUPPORT_OFFSET Locate the support point relative to its centre.
localDirection = rotation.'*direction;
denominator = norm(semiAxes.*localDirection);
offset = rotation*((semiAxes.^2).*localDirection)/denominator;
end


function rotation = rotation_y(angle)
%ROTATION_Y Return a right-handed rotation about the y axis.
rotation = [cos(angle),0,-sin(angle); 0,1,0; ...
    sin(angle),0,cos(angle)];
end


function rotation = rotation_z(angle)
%ROTATION_Z Return a right-handed rotation about the z axis.
rotation = [cos(angle),-sin(angle),0; ...
    sin(angle),cos(angle),0; 0,0,1];
end


function assert_rotation(rotation)
%ASSERT_ROTATION Check the orthogonality and handedness of one frame.
assert(norm(rotation.'*rotation-eye(3),'fro') < 100*eps);
assert(abs(det(rotation)-1) < 100*eps);
end
