function A = stokes_mfs_traction_matrix(basis,x,n)
%STOKES_MFS_TRACTION_MATRIX Mixed image-source traction at unit strength.
% Rows and each basis block use node-interleaved Cartesian components.
if nargin == 0
    self_test(); A = []; return
end
M = size(x,1);
C = 3*sum(arrayfun(@(b) size(b.points,1),basis));
A = zeros(3*M,C);
col0 = 0;
for b = 1:numel(basis)
    block = basis(b);
    y = block.points;
    N = size(y,1);
    r = cell(1,3);
    for j = 1:3
        r{j} = x(:,j)-y(:,j).';
    end
    r2 = r{1}.^2+r{2}.^2+r{3}.^2;
    if any(r2==0,'all')
        error('stokes_mfs_traction_matrix:CoincidentPoints', ...
            'A source coincides with a target.');
    end
    ir3 = r2.^(-3/2);
    ir5 = ir3./r2;
    ir7 = ir5./r2;
    rn = r{1}.*n(:,1)+r{2}.*n(:,2)+r{3}.*n(:,3);
    if block.family=="S_y"
        cols = col0+(1:3*N);
        A(:,cols) = block.scale*stokes_directional_stokeslet_matrix( ...
            y,x,block.normalDirections,block.derivativeDistances, ...
            "traction",n);
        col0 = col0+3*N;
        continue
    end
    if block.family=="Q_s"
        cols = col0+(1:3*N);
        A(:,cols) = block.scale*stokes_axial_quadrupole_matrix( ...
            y,x,block.normalDirections,"traction",n);
        col0 = col0+3*N;
        continue
    end
    if ismember(block.family,["T","T_n","T_n2","T_t"])
        q = block.normalDirections;
        if ~isequal(size(q),[N,3])
            error('stokes_mfs_traction_matrix:BadDirection', ...
                'Each stresslet needs one normal direction.');
        end
        a = r{1}.*q(:,1).'+r{2}.*q(:,2).'+r{3}.*q(:,3).';
        qn = n*q.';
    end
    for k = 1:3
        for j = 1:3
            switch block.family
                case "S"
                    value = -3/(4*pi)*r{k}.*rn.*r{j}.*ir5;
                case "R"
                    w = zeros(M,N);
                    wn = zeros(M,N);
                    for h = 1:3
                        v = rotlet_component(h,k,r);
                        wn = wn+n(:,h).*v;
                        if h == j, w = v; end
                    end
                    value = -3/(8*pi)*(w.*rn+r{j}.*wn).*ir5;
                case {"T","T_n","T_n2","T_t"}
                    value = -q(:,k).'.*n(:,j).*ir3/(2*pi) ...
                        -3/(4*pi)*(r{j}.*(qn.*r{k}+n(:,k).*a) + ...
                        (q(:,j).'.*r{k}+(j==k)*a).*rn).*ir5 ...
                        +15/(2*pi)*a.*r{k}.*r{j}.*rn.*ir7;
                case "D"
                    value = 3/(2*pi)*((j==k)*rn+r{j}.*n(:,k)+ ...
                        r{k}.*n(:,j)).*ir5- ...
                        15/(2*pi)*r{k}.*r{j}.*rn.*ir7;
                otherwise
                    error('stokes_mfs_traction_matrix:BadFamily', ...
                        'Unknown source family.');
            end
            A(j:3:3*M,col0+k:3:col0+3*N) = block.scale*value;
        end
    end
    col0 = col0+3*N;
end
end

function v = rotlet_component(j,k,r)
% Return one component of the rotlet cross product for a unit coefficient.
if j==k
    v = zeros(size(r{1}));
elseif j==1 && k==2 || j==2 && k==3 || j==3 && k==1
    v = r{6-j-k};
else
    v = -r{6-j-k};
end
end

function self_test()
% Compare every column of the selected families with direct physical traction.
base = struct('family',"S",'points',[0.1 0.2 -0.1], ...
    'normalDirections',zeros(0,3),'scale',1.3);
basis = repmat(base,8,1);
families = ["S","R","T","D","T_n","T_t","T_n2","Q_s"];
for k = 1:numel(families)
    basis(k).family = string(families(k));
    basis(k).points = [0.1*k 0.2 -0.1];
    if ismember(k,[3 5 6 7 8]), basis(k).normalDirections = [0 1 0]; end
end
x = [1.1 0.2 0.3;0.5 1.2 0.4];
n = x./vecnorm(x,2,2);
A = stokes_mfs_traction_matrix(basis,x,n);
grid.body.basis = basis;
for k = 1:size(A,2)
    z = zeros(size(A,2),1); z(k)=1;
    t = image_traction(grid,{z},x,n,1);
    assert(norm(A(:,k)-reshape(t.',[],1),inf)<2e-12);
end
% Integrate each kernel separately: enclosing a singularity changes moments.
[z,w] = gauss(64);
phi = (0:127)*(2*pi/128);
x = [reshape(sqrt(1-z.^2)*cos(phi),[],1), ...
    reshape(sqrt(1-z.^2)*sin(phi),[],1),repmat(z,128,1)];
w = repmat(w(:)*(2*pi/128),128,1);
maxMoment = 0;
for position = [0.2,1.4]
    for k = 1:numel(basis)
        basis(k).points = [position 0 0];
    end
    A = stokes_mfs_traction_matrix(basis,x,x);
    for k = 1:size(A,2)
        b = ceil(k/3); c = zeros(1,3); c(mod(k-1,3)+1) = basis(b).scale;
        t = reshape(A(:,k),3,[]).';
        force = sum(w.*t,1);
        torque = sum(w.*cross(x,t,2),1);
        expectedForce = zeros(1,3); expectedTorque = zeros(1,3);
        if position<1
            if basis(b).family=="S"
                expectedForce = -c;
                expectedTorque = -cross(basis(b).points,c);
            elseif basis(b).family=="R"
                expectedTorque = -c;
            end
        end
        maxMoment = max([maxMoment,norm(force-expectedForce), ...
            norm(torque-expectedTorque)]);
    end
end
if maxMoment>5e-10
    error('stokes_mfs_traction_matrix:SelfTest','Kernel force/torque balance failed.');
end
derivativeMoment = basis(1);
derivativeMoment.family = "S_y";
derivativeMoment.points = [0.2 0 0];
derivativeMoment.normalDirections = [0 1 0];
derivativeMoment.derivativeDistances = 0.1;
derivativeMoment.scale = 1;
t = reshape(stokes_mfs_traction_matrix(derivativeMoment,x,x)* ...
    [0;0;1],3,[]).';
assert(norm(sum(w.*t,1))<5e-10);
assert(norm(sum(w.*cross(x,t,2),1)-[-0.1 0 0])<5e-10);
derivative = basis(1);
derivative.family = "S_y";
derivative.points = [0.2 0 0];
derivative.normalDirections = [0 1 0];
derivative.derivativeDistances = 0.1;
targets = [1 0 0;0 1 0];
assert(norm(stokes_mfs_traction_matrix(derivative,targets,targets)- ...
    derivative.scale*stokes_directional_stokeslet_matrix(derivative.points,targets, ...
    derivative.normalDirections,derivative.derivativeDistances, ...
    "traction",targets),'fro')<1e-14);
fprintf('stokes_mfs_traction_matrix self-test passed; moment error %.3e.\n',maxMoment);
end
