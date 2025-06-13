function [position_init] = set_position(npart,L,tol,dim,offset)
%SET_POSITION(npart,L,tol,dim,offset) randomizes spherical particles in a box of
%length L not closer to each other than tol. No smart algorithm is used. If dim < 2, 
% the spheres are organized in a "layer" with some offset from the plane.

if nargin <5
    dim = 3
end


position_init = [];
maxit = 20;

for k = 1:npart
    if dim == 3
        new_pos = [1+(L-2*1)*rand(1,1) 1+(L-2*1)*rand(1,1) ...
            1+(L-2*1)*rand(1,1)];
    else
        new_pos = [1+(L-2*1)*rand(1,1) 1+(L-2*1)*rand(1,1) ...
            -offset+2*offset*rand(1,1)];
    end

    itr = 0;

    while isclose(new_pos,position_init,tol,L) %&& itr < maxit 

        if dim == 3
            new_pos = [1+(L-2*1)*rand(1,1) 1+(L-2*1)*rand(1,1) ...
                1+(L-2*1)*rand(1,1)];
        else
            new_pos = [1+(L-2*1)*rand(1,1) 1+(L-2*1)*rand(1,1) ...
                -offset+2*offset*rand(1,1)];
        end
        itr = itr+1;

    end

    position_init = [position_init; new_pos];

end



function res = isclose(new_pos,position_init,tol,L)

    for a = 1: size(position_init,1)

        if norm(position_init(a,:)-new_pos)<2+tol

            res = 1;

            return;

        end

    end

    res = 0;

end





end