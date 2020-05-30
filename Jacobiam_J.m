function J = Jacobiam_J(elementType,nodecoordinates,xi,eta)
    if strcmp(elementType,'T3')
        J = [1 0 -1; 0 1 -1] * nodecoordinates';
    elseif strcmp(elementType,'T6')
        Nxi  = [4*xi-1 ,0       , 4*eta+4*xi-3, 4*eta , -4*eta       ,4-8*xi-4*eta];
        Neta = [0      ,4*eta-1 , 4*eta+4*xi-3, 4*xi  , 4-4*xi-8*eta ,-4*xi       ];
        J = [Nxi;Neta] * nodecoordinates';
    elseif strcmp(elementType,'Q4')
        Nxi  = 0.25 * [eta-1 , 1-eta , 1+eta , -1-eta ];
        Neta = 0.25 * [xi-1  , -1-xi , 1+xi  , 1-xi   ];
        J = [Nxi;Neta] * nodecoordinates';
    elseif strcmp(elementType,'Q8')
        Nxi = zeros(1,8);
        Neta = zeros(1,8);
        Nxi(1) = -0.25*(eta+2*xi)*(eta-1);
        Nxi(2) = 0.25*(eta-2*xi)*(eta-1);
        Nxi(3) = 0.25*(eta+2*xi)*(eta+1);
        Nxi(4) = -0.25*(eta-2*xi)*(eta+1);
        Nxi(5) = xi*(eta-1);
        Nxi(6) = 0.5-0.5*eta^2;
        Nxi(7) = -xi*(eta+1);
        Nxi(8) = 0.5*eta^2-0.5;
        Neta(1) = -0.25*(2*eta+xi)*(xi-1);
        Neta(2) = 0.25*(2*eta-xi)*(xi+1);
        Neta(3) = 0.25*(2*eta+xi)*(xi+1);
        Neta(4) = -0.25*(2*eta-xi)*(xi-1);
        Neta(5) = 0.5*xi^2-0.5;
        Neta(6) = -eta*(xi+1);
        Neta(7) = 0.5-0.5*xi^2;
        Neta(8) = eta*(xi-1);
        J = [Nxi;Neta] * nodecoordinates';
    end
end