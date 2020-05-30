function [K,E] = stiffness_2D(materialObj,meshObj,elementObj,state)
young = materialObj.young;
enu = materialObj.enu;
m = meshObj.m;
n = meshObj.n;
x = meshObj.x;
node = meshObj.node;


%% Judge the state of stress and strain and compute the E matrix
if strcmp(state,'Plane Stress')
    E = (young/(1-enu^2))*[1,enu,0; ...
                           enu,1,0; ...
                           0,0,(1-enu)/2];
elseif strcmp(state,'Plane Strain') 
    E = [1-enu, enu, 0; ...
        enu, 1-enu,0; ...
        0, 0, (1-2*enu)/2] * (young/((1+enu)*(1-2*enu)));
end
%% loop for each elements
K = sparse(n*2,n*2);
for i = 1:m
    B = @(xi,eta) elementObj.StrainMatrix(meshObj,i,xi,eta);
    J = @(xi,eta) elementObj.Jacobian(meshObj,i,xi,eta);
    fun =@(xi,eta) B(xi,eta)' * E * B(xi,eta) * det(J(xi,eta));
    % here xi and eta must be a scalar !!!
    % Gauss Quadrature for triangles
    [intPoint,weight] = GQintegration(elementObj);
    testDem = B(0.5,0.5);
    K_elem = zeros(size(testDem,2));
    for j = 1:length(weight)
        K_elem = fun(intPoint(j,1),intPoint(j,2)) .* weight(j) + K_elem;
    end
    rowcol = zeros(2*length(node(:,i)),1);
    rowcol(1:2:end-1) = 2*node(:,i)-1;
    rowcol(2:2:end) = 2*node(:,i);
    K(rowcol,rowcol) = K(rowcol,rowcol) + K_elem ;
    
end