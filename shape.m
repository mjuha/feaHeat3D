function [N,dN,jac,inv_jac] = shape(gp,xe)

% local coordinate
r = gp(1);
s = gp(2);
t = gp(3);

% shape functions
N = [ 1-r-s-t, r, s, t ];
N_r = [ -1, 1, 0, 0 ];
N_s = [-1, 0, 1, 0];
N_t = [-1, 0, 0, 1];

dN = zeros(4,3);

x_r = N_r * xe(:,1);
x_s = N_s * xe(:,1);
x_t = N_t * xe(:,1);
%
y_r = N_r * xe(:,2);
y_s = N_s * xe(:,2);
y_t = N_t * xe(:,2);
%
z_r = N_r * xe(:,3);
z_s = N_s * xe(:,3);
z_t = N_t * xe(:,3);

jacobian = [x_r, x_s, x_t;  y_r, y_s, y_t; z_r, z_s, z_t];
jac = det(jacobian);

% check jacobian
if jac < 1.0e-14
    error('Negative jacobian, element too distorted!');
end

inv_jac = inv(jacobian);
% Note: inv_jac = [r_x, r_y; s_x, s_y]

for i=1:4
    dN(i,:) = [N_r(i), N_s(i), N_t(i)] * inv_jac; %#ok<MINV>
end

end
