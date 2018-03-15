function feaHeat3D(filename,varargin)
% feaHeat3D.m
% This is the main file that implement a tetrahedra
% element for solving transient heat transfer problems
%
% Author: Dr. Mario J. Juha
% Date: 11/09/2017
% Mechanical Engineering
% Universidad de La Sabana
% Chia -  Colombia
%
% Clear variables from workspace
clearvars

global nel neq nzmax coordinates U elements nn LM irow icol ID TS isTimeDBC

% Specify file name
%filename = '\\Client\C$\Users\marioju\Documents\Work\valvula\example.inp';

% check if user have specified rregression tests
if nargin > 1
    error('Only accept filename!')
%     % extract the regression values
%     normDisp = varargin{1};
%     normRot = varargin{2};
%     regressionTol = varargin{3};
end

% read data
fprintf('************************\n')
fprintf('Reading input data\n')
fprintf('************************\n\n')
outfile = readData(filename);

%WriteVTKFile(outfile,0)

% at time t = 0
fprintf('\nComputing time 0.0\n\n')
% ===========================
% assembling stiffness matrix
% ===========================
K = zeros(1,nzmax);
F = zeros(neq,1);
% set counter to zero
count = 0;
for i=1:nel
    xe = coordinates(elements(i,2:5),:);
    de = U(1,elements(i,2:5));
    matNum = elements(i,1); % element material number
    [fe,me] = weakform(i,xe,de',matNum);
    for k=1:4
        i_index = LM(k,i);
        if (i_index > 0)
            F(i_index) = F(i_index) + fe(k);
            for m=1:4
                j_index = LM(m,i);
                if (j_index > 0)
                    count = count + 1;
                    K(count) = me(k,m);
                end
            end
        end
    end
end
fprintf('************************\n')
fprintf('Solving system of equations\n')
fprintf('************************\n\n')
M = sparse(irow,icol,K,neq,neq);
F = M\F;
% assign solution
for r=1:nn
    i_index = ID(r);
    if (i_index > 0)
        U(2,r) = F(i_index);
    end
end
WriteVTKFile(outfile,0)



% fprintf('************************\n')
% fprintf('Assembling stiffness matrix and force vector\n')
% fprintf('************************\n\n')
fprintf('************************\n')
fprintf('Advancing in time\n')
fprintf('************************\n\n')
%
t = 0;
dt = TS{1}; % delta time
tf = TS{2}; % final time
alpha = TS{3}; % method
nts = TS{4}; % ouput every nts
counter = 0;
count1 = 0;
while ( t <= tf )
    counter = counter + 1;
    t = t + dt;
    fprintf('Computing time %f\n\n',t)
    % ===========================
    % assembling stiffness matrix
    % ===========================
    K = zeros(1,nzmax);
    F = zeros(neq,1);
    % set counter to zero
    count = 0;
    if isTimeDBC
        DBC_InTime(t)
    end
    % predictor value
    dp = U(1,:) + (1-alpha)*dt*U(2,:);
    for i=1:nel
        xe = coordinates(elements(i,2:5),:);
        de = dp(elements(i,2:5));
        matNum = elements(i,1); % element material number
        [fe,me] = weakform(i,xe,de',matNum);
        for k=1:4
            i_index = LM(k,i);
            if (i_index > 0)
                F(i_index) = F(i_index) + fe(k);
                for m=1:4
                    j_index = LM(m,i);
                    if (j_index > 0)
                        count = count + 1;
                        K(count) = me(k,m);
                    end
                end
            end
        end
    end
    fprintf('************************\n')
    fprintf('Solving system of equations\n')
    fprintf('************************\n\n')
    M = sparse(irow,icol,K,neq,neq);
    F = M\F;
    % assign solution
    for r=1:nn
        i_index = ID(r);
        if (i_index > 0)
            U(2,r) = F(i_index);
            % corrector
            U(1,r) = dp(r) + alpha*dt*U(2,r);
        end
    end
    if mod(counter,nts) == 0
        count1 = count1 + 1; 
        WriteVTKFile(outfile,count1)
    end
end
% fprintf('************************\n')
% fprintf('Computing stresses and printing final mesh and results\n')
% fprintf('************************\n\n')
% 
% computeStressStrain
% 
% WriteVTKFile(outfile,1)
end
