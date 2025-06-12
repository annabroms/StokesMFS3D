function res = rotate_vector(x,R)

x_reshape = reshape(x, 3, []);
rot = R*x_reshape;

res = rot(:);

end