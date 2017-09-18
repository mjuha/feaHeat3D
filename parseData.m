function parseData(nodeSet, sideSet)

% global variables
global elements nn nel triNode ID LM fluxLoad convectionLoad
global U neq
global DBCSet HFCSet NBCSet

% post-process boundary conditions
% ===============
% ID array
% ===============
ID = ones(1,nn);
%================

[m,~] = size(triNode);
numericCell = nodeSet(:,1);
numericVector = cell2mat(numericCell);
count = 0;
% check nodeSet
for i=1:m
    phytag = triNode(i,1);
    % search in nodeSet
    [row,~] = find(numericVector==phytag);
    if  isempty(row)
        continue
    end
    % find in DBCSet first
    name = nodeSet{row,2};
    % find in DBCSet first
    [row,~] = find(strcmp(DBCSet,name),3);
    found = false;
    for j=1:length(row)
        dof = DBCSet{row(j),2};
        if strcmp(dof,'T')
            value = DBCSet{row(j),3};
            U(1,triNode(i,2:4)) = value;       
        elseif strcmp(dof,'TFunction')
            value = DBCSet{row(j),3};
            U(1,triNode(i,2:4)) = value(2);       
        else
            error('Unknown dof!')
        end
        ID(triNode(i,2:4)) = 0; % T
        found = true;
    end
    if found
        count = count + 1;
    end
end

% for flux BCs
[row,~] = size(HFCSet);
if row > 0
    fluxLoad = zeros(m-count,6);
    numericCell = sideSet(:,1);
    numericVector = cell2mat(numericCell);
    count = 1;
    % check sideSet
    for i=1:m
        phytag = triNode(i,1);
        % search in sideSet
        [row,~] = find(numericVector==phytag);
        if  isempty(row)
            continue
        end
        % find in NBCSet first
        name = sideSet{row,2};
        % find in NBCSet first
        [row,~] = find(strcmp(HFCSet,name),3);
        for j=1:length(row)
            value = HFCSet{row(j),2};
            fluxLoad(count,3:5) = triNode(i,2:4);
            fluxLoad(count,6) = value;
            count = count + 1;
        end
    end
else
    fluxLoad = [];
end

[row,~] = size(NBCSet);
if row > 0
    % for convection BCs
    convectionLoad = zeros(m-count,7);
    numericCell = sideSet(:,1);
    numericVector = cell2mat(numericCell);
    % check sideSet
    count = 1;
    for i=1:m
        phytag = triNode(i,1);
        % search in sideSet
        [row,~] = find(numericVector==phytag);
        if  isempty(row)
            continue
        end
        % find in NBCSet first
        name = sideSet{row,2};
        % find in NBCSet first
        [row,~] = find(strcmp(NBCSet,name),3);
        for j=1:length(row)
            coeff = NBCSet{row(j),2};
            value = NBCSet{row(j),3};
            convectionLoad(count,3:5) = triNode(i,2:4);
            convectionLoad(count,6) = coeff;
            convectionLoad(count,7) = value;
            count = count + 1;
        end
    end
else
    convectionLoad = [];
end


% find element where convection is applied
[m,~] = size(convectionLoad);
if m > 0
    fprintf('\nFinding elements where convection is applied...\n');
end
for i=1:m
    A = convectionLoad(i,3:5);
    for j=1:nel
        B = elements(j,2:5);
        C = ismember(B,A);
        if sum(C) == 3
            convectionLoad(i,1) = j;
            break
        end
    end
end
% find face where convection is applied
if m > 0
    fprintf('\nFinding faces where convection is applied...\n');
end
for i=1:m
    A = convectionLoad(i,3:5);
    j = convectionLoad(i,1);
    B = elements(j,2:5);
    C = ismember(B,A);
    tmp = find(C);
    % find face
    if ( tmp(1)== 1 && tmp(2) == 2 && tmp(3) == 3  ) || ...
            ( tmp(1)== 3 && tmp(2) == 1 && tmp(3) == 2 ) || ...
            ( tmp(1)== 2 && tmp(2) == 3 && tmp(3) == 1  ) || ...
            ( tmp(1)== 1 && tmp(2) == 3 && tmp(3) == 2  ) || ...
            ( tmp(1)== 3 && tmp(2) == 2 && tmp(3) == 1  ) || ...
            ( tmp(1)== 2 && tmp(2) == 1 && tmp(3) == 3  )
        convectionLoad(i,2) = 1; % face 1
    elseif ( tmp(1)== 1 && tmp(2) == 2 && tmp(3) == 4  ) || ...
            ( tmp(1)== 1 && tmp(2) == 4 && tmp(3) == 2 ) || ...
            ( tmp(1)== 2 && tmp(2) == 1 && tmp(3) == 4  ) || ...
            ( tmp(1)== 2 && tmp(2) == 4 && tmp(3) == 1  ) || ...
            ( tmp(1)== 4 && tmp(2) == 1 && tmp(3) == 2  ) || ...
            ( tmp(1)== 4 && tmp(2) == 2 && tmp(3) == 1  )
        convectionLoad(i,2) = 2; % face 2
    elseif ( tmp(1)== 1 && tmp(2) == 3 && tmp(3) == 4  ) || ...
            ( tmp(1)== 1 && tmp(2) == 4 && tmp(3) == 3 ) || ...
            ( tmp(1)== 3 && tmp(2) == 1 && tmp(3) == 4  ) || ...
            ( tmp(1)== 3 && tmp(2) == 4 && tmp(3) == 1  ) || ...
            ( tmp(1)== 4 && tmp(2) == 1 && tmp(3) == 3  ) || ...
            ( tmp(1)== 4 && tmp(2) == 3 && tmp(3) == 1  )
        convectionLoad(i,2) = 3; % face 3
    elseif ( tmp(1)== 2 && tmp(2) == 3 && tmp(3) == 4  ) || ...
            ( tmp(1)== 2 && tmp(2) == 4 && tmp(3) == 3 ) || ...
            ( tmp(1)== 3 && tmp(2) == 2 && tmp(3) == 4  ) || ...
            ( tmp(1)== 3 && tmp(2) == 4 && tmp(3) == 2  ) || ...
            ( tmp(1)== 4 && tmp(2) == 2 && tmp(3) == 3  ) || ...
            ( tmp(1)== 4 && tmp(2) == 3 && tmp(3) == 2  )
        convectionLoad(i,2) = 4; % face 4
    else
        error('Edge not found!\n')
    end
end
%
[row,~] = size(NBCSet);
if row > 0
    index = find(~convectionLoad(:,1),1);
    if ~isempty(index)
        convectionLoad(index:end,:) = [];
    end
end

% find element where flux is applied
[m,~] = size(fluxLoad);
if m > 0
    fprintf('\nFinding elements where heat flux is applied...\n');
end
for i=1:m
    A = fluxLoad(i,3:5);
    for j=1:nel
        B = elements(j,2:5);
        C = ismember(B,A);
        if sum(C) == 3
            fluxLoad(i,1) = j;
            break
        end
    end
end
% find face where flux is applied
if m > 0
    fprintf('\nFinding faces where heat flux is applied...\n');
end
for i=1:m
    A = fluxLoad(i,3:5);
    j = fluxLoad(i,1);
    B = elements(j,2:5);
    C = ismember(B,A);
    tmp = find(C);
    % find face
    if ( tmp(1)== 1 && tmp(2) == 2 && tmp(3) == 3  ) || ...
            ( tmp(1)== 3 && tmp(2) == 1 && tmp(3) == 2 ) || ...
            ( tmp(1)== 2 && tmp(2) == 3 && tmp(3) == 1  ) || ...
            ( tmp(1)== 1 && tmp(2) == 3 && tmp(3) == 2  ) || ...
            ( tmp(1)== 3 && tmp(2) == 2 && tmp(3) == 1  ) || ...
            ( tmp(1)== 2 && tmp(2) == 1 && tmp(3) == 3  )
        fluxLoad(i,2) = 1; % face 1
    elseif ( tmp(1)== 1 && tmp(2) == 2 && tmp(3) == 4  ) || ...
            ( tmp(1)== 1 && tmp(2) == 4 && tmp(3) == 2 ) || ...
            ( tmp(1)== 2 && tmp(2) == 1 && tmp(3) == 4  ) || ...
            ( tmp(1)== 2 && tmp(2) == 4 && tmp(3) == 1  ) || ...
            ( tmp(1)== 4 && tmp(2) == 1 && tmp(3) == 2  ) || ...
            ( tmp(1)== 4 && tmp(2) == 2 && tmp(3) == 1  )
        fluxLoad(i,2) = 2; % face 2
    elseif ( tmp(1)== 1 && tmp(2) == 3 && tmp(3) == 4  ) || ...
            ( tmp(1)== 1 && tmp(2) == 4 && tmp(3) == 3 ) || ...
            ( tmp(1)== 3 && tmp(2) == 1 && tmp(3) == 4  ) || ...
            ( tmp(1)== 3 && tmp(2) == 4 && tmp(3) == 1  ) || ...
            ( tmp(1)== 4 && tmp(2) == 1 && tmp(3) == 3  ) || ...
            ( tmp(1)== 4 && tmp(2) == 3 && tmp(3) == 1  )
        fluxLoad(i,2) = 3; % face 3
    elseif ( tmp(1)== 2 && tmp(2) == 3 && tmp(3) == 4  ) || ...
            ( tmp(1)== 2 && tmp(2) == 4 && tmp(3) == 3 ) || ...
            ( tmp(1)== 3 && tmp(2) == 2 && tmp(3) == 4  ) || ...
            ( tmp(1)== 3 && tmp(2) == 4 && tmp(3) == 2  ) || ...
            ( tmp(1)== 4 && tmp(2) == 2 && tmp(3) == 3  ) || ...
            ( tmp(1)== 4 && tmp(2) == 3 && tmp(3) == 2  )
        fluxLoad(i,2) = 4; % face 4
    else
        error('Edge not found!\n')
    end
end

%
[row,~] = size(HFCSet);
if row > 0
    index = find(~fluxLoad(:,1),1);
    if ~isempty(index)
        fluxLoad(index:end,:) = [];
    end
end

% Fill ID array
count = 0;
for j=1:nn
    if ID(j) ~= 0
        count = count + 1;
        ID(j) = count;
    end
end

% number of equation
neq = max(ID);

% =================
% Generate LM array
% =================
LM = zeros(4,nel);
for k=1:nel
    for j=1:4
        LM(j,k) = ID(elements(k,j+1));
    end
end

ComputeSparsity

end

