function WriteVTKFile( outfiledest,istep  )
% ======================================================================
% This file is part of feaHeat3D.

%    feaTri is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as
%    published by the Free Software Foundation, either version 3 of the
%    License, or (at your option) any later version.

%    feaTri is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.

%    You should have received a copy of the GNU Lesser General Public
%    License along with feaTri.
%    If not, see <http://www.gnu.org/licenses/>.
%
% Brief explanation:
%
%
% Author: Mario J. Juha, Ph.D
% Date:
%
% =====================================================================

global coordinates elements nn nel U

if istep < 10
    % file name
    fname = [outfiledest 'output000' num2str(istep) '.vtk'];
elseif istep < 100
    % file name
    fname = [outfiledest 'output00' num2str(istep) '.vtk'];
elseif istep < 1000
    % file name
    fname = [outfiledest 'output0' num2str(istep) '.vtk'];
else
    % file name
    fname = [outfiledest 'output' num2str(istep) '.vtk'];
end

% open file
fid = fopen(fname, 'w');

fprintf(fid, '# vtk DataFile Version 3.8\n');
fprintf(fid, 'Mesh\n');
fprintf(fid,'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid, '%s %d %s\n','POINTS ', nn, 'float');
for i=1:nn
    fprintf(fid, '%g  %g  %g\n', coordinates(i,1), coordinates(i,2), coordinates(i,3));
end
fprintf(fid, '%s %d %d\n','CELLS ', nel, 5*nel);
for i=1:nel
    fprintf(fid, '%d %d %d %d %d\n',4, elements(i,2:5)-1);
end
fprintf(fid, '%s %d\n','CELL_TYPES ', nel);
for i=1:nel
    fprintf(fid, '%d\n', 10);
end
fprintf(fid, '%s %d\n', 'POINT_DATA ', nn);
fprintf(fid, 'SCALARS T float 1\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
for i=1:nn
    fprintf(fid, '%g\n',U(1,i));
end
% fprintf(fid, '%s %d\n', 'CELL_DATA ', nel);
% fprintf(fid, 'TENSORS Stress float\n');
% %fprintf(fid, 'LOOKUP_TABLE default\n');
% switch STATE
%     case 'planeStress'
%         for i=1:nel
%             fprintf(fid, '%g %g %g\n',stress(1,i), stress(3,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',stress(3,i), stress(2,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',0.0, 0.0, 0.0 );
%         end
%     case 'planeStrain'
%         poisson = MAT(1,2);
%         for i=1:nel
%             fprintf(fid, '%g %g %g\n',stress(1,i), stress(3,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',stress(3,i), stress(2,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',0.0, 0.0, (1+poisson)*(stress(1,i) ...
%                 + stress(2,i)));
%         end
%     otherwise
%         error('Unknown state!\n')
% end
% 
% fprintf(fid, 'TENSORS Strain float\n');
% %fprintf(fid, 'LOOKUP_TABLE default\n');
% switch STATE
%     case 'planeStress'
%         poisson = MAT(1,2);
%         young = MAT(1,1);
%         for i=1:nel
%             fprintf(fid, '%g %g %g\n',strain(1,i), strain(3,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',strain(3,i), strain(2,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',0.0, 0.0, -(poisson/young)*...
%                 (stress(1,i)+stress(2,i)));
%         end
%     case 'planeStrain'
%         for i=1:nel
%             fprintf(fid, '%g %g %g\n',strain(1,i), strain(3,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',strain(3,i), strain(2,i), 0.0 );
%             fprintf(fid, '%g %g %g\n',0.0, 0.0, 0.0);
%         end
%     otherwise
%         error('Unknown state!\n')
% end

% close file
fclose(fid);


end

