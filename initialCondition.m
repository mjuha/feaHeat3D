function initialCondition(IC)

global U nn

% solution vector
U = zeros(2,nn);

name = IC{1};

if strcmp(name,'Constant')
    value = IC{2};
end

for i = 1:nn
    U(1,i) = value;
end


end