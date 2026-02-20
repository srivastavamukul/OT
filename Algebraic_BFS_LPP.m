clc 
clear all

C = [2 3 4 7];
A = [ 2 3 -1 4;1 -2 6 -7]
B = [8 ;-3]

n = size(A,2);
m = size(A,1);

ncm = nchoosek(n,m);

if n<m 
    error( ...
        )
end

pair = nchoosek(1:n,m);

sol= [];
for i = 1 : ncm
    y = zeros(n,1);
    X = A(:,pair(i,:)) \ B;
    if all(X>0 & X~=inf & X~=-inf)
        y(pair(i,:)) = X;
        sol =[sol y];
    end
end

z = C*sol;

[zmax,zindex] = max(z);

bfs = sol(:,zindex);
zmax;


optimal_Value = [bfs ; zmax]
optimalbfs = array2table(optimal_Value.')
optimalbfs.Properties.VariableNames = {'x1','x2','x3','x4','z'}
