function ComputeSparsity
% ======================================================================
% This file is part of FlowSolveTri.

%    FlowSolveTri is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.

%    FloSolveTri is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.

%    You should have received a copy of the GNU Lesser General Public 
%    License along with Foobar.  
%    If not, see <http://www.gnu.org/licenses/>.
%
% Brief explanation:
%  FlowSolveTri is based on a stabilized finite element implementation
%  of the Navier-Stokes equation that uses linear triangular elements
%  with equal interpolation for velocities and pressures.
%
%
% Author: Mario J. Juha, Ph.D
% Date: 04/28/2014
% Tampa - Florida
% =====================================================================
global LM nel irow icol nzmax

nzmax = 0;
for elem=1:nel
    for k=1:4
        i_index = LM(k,elem);
        if (i_index > 0)
            for m=1:4
                j_index = LM(m,elem);
                if (j_index > 0)
                    nzmax = nzmax + 1;
                end
            end
        end
    end
end

irow = zeros(1,nzmax);
icol = zeros(1,nzmax);

count = 0;
for elem=1:nel
    for k=1:4
        i_index = LM(k,elem);
        if (i_index > 0)
            for m=1:4
                j_index = LM(m,elem);
                if (j_index > 0)
                    count = count + 1;
                    irow(count) = i_index;
                    icol(count) = j_index;
                end
            end
        end
    end
end

end

