function B = compute_B_matrix(elementType,elementnumber,nodecoordinates,xi,eta)

num = [1 0 0 0 ; 0 0 0 1 ; 0 1 1 0];

if strcmp(elementType,'T3')
    A = 0.5*det([1,1,1;nodecoordinates]);
    x = nodecoordinates(1,:);
    y = nodecoordinates(2,:);
    B = (1/(2*A))*[y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0; ...
        0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1); ...
        x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];
    
elseif strcmp(elementType,'T6')

    Nxi  = [4*xi-1 ,0       , 4*eta+4*xi-3, 4*eta , -4*eta       ,4-8*xi-4*eta];
    Neta = [0      ,4*eta-1 , 4*eta+4*xi-3, 4*xi  , 4-4*xi-8*eta ,-4*xi       ];
    J = [Nxi;Neta] * nodecoordinates';
    if det(J) == 0
        fprintf('Jacobian matrix is sigular at element %i \n',elementnumber)
    end
    Gamma = inv(J);
    G = [Gamma zeros(2);zeros(2) Gamma];
    Ndev = zeros(4,12);
    Ndev(1,1:2:end-1) = Nxi;
    Ndev(2,1:2:end-1) = Neta;
    Ndev(3,2:2:end) = Nxi;
    Ndev(4,2:2:end) = Neta;
    B = num * G * Ndev;
    
elseif strcmp(elementType,'Q4')

    Nxi  = 0.25 * [eta-1 , 1-eta , 1+eta , -1-eta ];
    Neta = 0.25 * [xi-1  , -1-xi , 1+xi  , 1-xi   ];
    J = [Nxi;Neta] * nodecoordinates';
    if det(J) == 0
        fprintf('Jacobian matrix is sigular at element %i \n',elementnumber)
    end
    Gamma = inv(J);
    G = [Gamma zeros(2);zeros(2) Gamma];
    Ndev = zeros(4,8);
    Ndev(1,1:2:end-1) = Nxi;
    Ndev(2,1:2:end-1) = Neta;
    Ndev(3,2:2:end) = Nxi;
    Ndev(4,2:2:end) = Neta;
    B = num * G * Ndev;
    
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
    if det(J) == 0
        fprintf('Jacobian matrix is sigular at element %i \n',elementnumber)
    end
    Gamma = inv(J);
    G = [Gamma zeros(2);zeros(2) Gamma];
    Ndev = zeros(4,16);
    Ndev(1,1:2:end-1) = Nxi;
    Ndev(2,1:2:end-1) = Neta;
    Ndev(3,2:2:end) = Nxi;
    Ndev(4,2:2:end) = Neta;
    B = num * G * Ndev;
    
end

end
