function DBC_InTime(t)

global nodeSet triNode DBCSet U TS

[m,~] = size(triNode);
numericCell = nodeSet(:,1);
numericVector = cell2mat(numericCell);
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
    for j=1:length(row)
        dof = DBCSet{row(j),2};   
        if strcmp(dof,'TFunction')
            value = DBCSet{row(j),3};
            if value(1) == 1
                val = TLinear(t,value(2),value(3),TS{2});
            end
            U(1,triNode(i,2:4)) = val;       
        end
    end
end


end