function [fe] = computeSideLoad(el,xe)

global convectionLoad

fe = zeros(4,1);

% Gauss - Legendre rule
gp = [1/3, 1/3];
w = 1.0;

% compute residual: side loads
face = convectionLoad(el,2);
h = convectionLoad(el,6); % coefficient
Ta = convectionLoad(el,7); % ambient temperature

% loop over gauss points
if face == 1 % local nodes 1-2-3
    r = gp(1);
    s = gp(2);
    % t = 0
    Nshape = [ 1-r-s, r, s, 0 ];
    N_r = [ -1, 1, 0, 0 ];
    N_s = [ -1, 0, 1, 0 ];
    N_t = [0, 0, 0, 0];
elseif face == 2 % local nodes 1-2-4
    r = gp(1);
    t = gp(2);
    % s = 0
    Nshape = [ 1-r-t, r, 0, t ];
    N_r = [ -1 ,1, 0, 0 ];
    N_s = [0, 0, 0, 0];
    N_t = [ -1 ,1, 0, 1 ];
elseif edge == 3 % local nodes 1-3-4
    % r = 0
    s = gp(1);
    t = gp(2);
    Nshape = [ 1-s-t, 0, s, t ];
    N_r = [0, 0, 0, 0];
    N_s = [ -1, 0, 1, 0 ];
    N_t = [ -1, 0, 0, 1 ];
elseif edge == 4 % local nodes 2-3-4
    % t = 1 - r - s
    r = gp(1);
    s = gp(2);
    Nshape = [ 0, r, s, 1-r-s ];
    N_r = [ 0, 1, 0, -1 ];
    N_s = [ 0, 0, 1, -1 ];
    N_t = [ 0, 0, 0, 0 ];
else
    error('Wrong face, check input!');
end

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

fe = fe + Nshape' * (h*Ta) * w * jac;

end