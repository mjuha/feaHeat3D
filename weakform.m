function [fe,me] = weakform(el,xe,de,matNum)

global convectionLoad MAT TS

dt = TS{1};
alpha = TS{3};

% 1 point formula - degree of precision 1
gp =  [ 1/4, 1/4, 1/4;];
w = 1/6;

% get material properties
prop = cell2mat(MAT(matNum));
rho = prop(1); % density
cp = prop(2); % heat capacity
k = prop(3); % conductivity

% stress-strain displacement matrix
B = zeros(3,4);
% loop over gauss points
[N,dN,jac] = shape(gp,xe);
for j=1:4 % loop over local nodes
    B(:,j) = dN(j,:);
end
ke = B' * k * B * w * jac;
%
me = N' * (rho*cp) * N * w * jac;

fe = zeros(4,1);
if size(convectionLoad,1) > 0
    index = find(convectionLoad(:,1)==el,4); % up to four faces
    flag = size(index,1);
    % compute side load
    if flag > 0
        faces = size(index,1);
        for i=1:faces
            [ke1,fe1] = computeSideLoad(index(i),xe);
            fe = fe1;
            ke = ke + ke1;
        end
    end
end

fe = fe - ke * de;
%
me = me + alpha*dt*ke;

end
