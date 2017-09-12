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

global nel neq nzmax coordinates U elements nn LM irow icol ID pointForce

% Specify file name
filename = '\\Client\C$\Users\marioju\Documents\Work\feaTransientHeat3D\example.inp';

% read data
fprintf('************************\n')
fprintf('Reading input data\n')
fprintf('************************\n\n')
outfile = readData(filename);

WriteVTKFile(outfile,0)

% fprintf('************************\n')
% fprintf('Assembling stiffness matrix and force vector\n')
% fprintf('************************\n\n')
% % ===========================
% % assembling stiffness matrix
% % ===========================
% K = zeros(1,nzmax);
% F = zeros(neq,1);
% % set counter to zero
% count = 0;
% for i=1:nel
%     xe = coordinates(elements(i,2:4),:);
%     de = U(:,elements(i,2:4));
%     [fe,ke] = weakform(i,xe,de);
%     for k=1:6
%         i_index = LM(k,i);
%         if (i_index > 0)
%             F(i_index) = F(i_index) + fe(k);
%             for m=1:6
%                 j_index = LM(m,i);
%                 if (j_index > 0)
%                     count = count + 1;
%                     K(count) = ke(k,m);
%                 end
%             end
%         end
%     end
% end
% % assign point force
% for i=1:size(pointForce,1)
%   i_index = ID(pointForce(i,2),pointForce(i,1));
%   F(i_index) = F(i_index) + pointForce(i,3);
% end
% fprintf('************************\n')
% fprintf('Solving system of equations\n')
% fprintf('************************\n\n')
% M = sparse(irow,icol,K,neq,neq);
% F = M\F;
% % assign solution
% for r=1:nn
%     for s=1:2
%         i_index = ID(s,r);
%         if (i_index > 0)
%             U(s,r) = F(i_index);
%         end
%     end
% end
% fprintf('************************\n')
% fprintf('Computing stresses and printing final mesh and results\n')
% fprintf('************************\n\n')
% 
% computeStressStrain
% 
% WriteVTKFile(outfile,1)
