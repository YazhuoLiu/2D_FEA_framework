function [stress,strain] = compute_stress_strain(elementObj,meshObj,u,E)
m = meshObj.m;
n = meshObj.n;
node = meshObj.node;

strain_x = sparse(n,m);
strain_y = sparse(n,m);
strain_xy = sparse(n,m);
stress_x = sparse(n,m);
stress_y = sparse(n,m);
stress_xy = sparse(n,m);
for i = 1:m
    B = @(xi,eta) elementObj.StrainMatrix(meshObj,i,xi,eta);
    xi = elementObj.Node_xi;
    eta = elementObj.Node_eta;
    d = zeros(2*length(node(:,i)),1);
    d(1:2:end-1) = u(2*node(:,i)-1);
    d(2:2:end) = u(2*node(:,i));
    for j = 1:length(xi)
        jNodeStrain = B(xi(j),eta(j))*d;
        strain_x(node(j,i),i) = jNodeStrain(1);
        strain_y(node(j,i),i) = jNodeStrain(2);
        strain_xy(node(j,i),i) = jNodeStrain(3);
        jNodeStress = E * jNodeStrain;
        stress_x(node(j,i),i) = jNodeStress(1);
        stress_y(node(j,i),i) = jNodeStress(2);
        stress_xy(node(j,i),i) = jNodeStress(3);
    end
end
stress = zeros(n,3);
strain = zeros(n,3);
stress(:,1) = SparseAverage(stress_x);
stress(:,2) = SparseAverage(stress_y);
stress(:,3) = SparseAverage(stress_xy);
strain(:,1) = SparseAverage(strain_x);
strain(:,2) = SparseAverage(strain_y);
strain(:,3) = SparseAverage(strain_xy);

end


function meanValue = SparseAverage(A)

A_full = full(A);
meanValue = zeros(size(A_full,1),1);
for i = 1:size(A_full,1)
    if A(i,:) == zeros(1,size(A,2))
        meanValue(i)=0;
    else
        L = A(i,:)~=0;
        A_row = A(i,L);
        meanValue(i) = mean(A_row);
    end
end

end
