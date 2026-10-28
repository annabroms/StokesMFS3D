function [D,maxResidual,maxIterations,maxForceError] = ...
    image_stress_divergence(q,opt,h)
%IMAGE_STRESS_DIVERGENCE Differentiate complete unit-force stresslet responses.
% opt.N_image and opt.nearThreshold must keep enrichment fixed throughout
% the stencil. Positions of image and cap nodes follow the perturbed geometry.
% D(:,:,p) = sum_(j,alpha) d(S_p M e_jalpha)/d q_jalpha.
% Body loads are [Fx;Fy;Fz;Tx;Ty;Tz]; nodal coefficients are interleaved xyz.
% This is the translational divergence; spheres have no orientation dependence.

P = size(q,1);
D = zeros(3,3,P);
maxResidual = 0;
maxIterations = 0;
maxForceError = 0;
for body = 1:P
    for axis = 1:3
        F = zeros(6*P,1);
        F(6*(body-1)+axis) = 1;
        stress = cell(2,1);
        for side = 1:2
            shiftedQ = q;
            shiftedQ(body,axis) = shiftedQ(body,axis)+(2*side-3)*h;
            [~,it,~,err,details] = solve_mobility_images(shiftedQ,F,opt);
            stress{side} = apply_KS_transpose(details.grid,details.lambda);
            maxResidual = max(maxResidual,err);
            maxIterations = max(maxIterations,it);
            maxForceError = max(maxForceError, ...
                norm(details.realisedForceTorque-F,inf));
        end
        D = D+(stress{2}-stress{1})/(2*h);
    end
end
end
