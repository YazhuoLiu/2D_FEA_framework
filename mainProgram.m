%% Plane FEM code
% A finite element method (FEM) code in Fortran for solving general 2-D 
% elastostatic problems using the T3\T6\Q4\Q8 element
% Author: Liu Yazhuo
% Date: 2020.4.30
clear all; clc
q = 1e7;
meshsz = 2;
state = 'Plane Strain';

spmd
    switch labindex
        case 1
            steel = materials(210e9,0.3);
        case 2
            elementObj = linear_strain_triangular(12);
        case 3
            fid=fopen('Report.txt','wt'); 
            fprintf(fid,'state: ');
            fprintf(fid,[state,'\n\n']);
            fclose(fid);
    end
end
steel = steel{1}; elementObj = elementObj{2};

meshObj = mesh(elementObj,meshsz);

spmd
    switch labindex
        case 1
            left_x = find(meshObj.x(1,:) == -1)'; % node numbers on left boundary
            f = sparse(2*meshObj.n,1);
            rols=2*left_x -1;
        case 2
            [K,E] = stiffness_2D(steel,meshObj,elementObj,state);
            K_bc = K;
        case 3
            boundaryElement = mesh.BoundaryElement(meshObj);
        case 4
            writeReport('Report.txt','node coordinates',meshObj.x');
            writeReport('Report.txt','element connectivity',meshObj.node');
    end
end
f = f{1}; rols = rols{1}; K = K{2}; K_bc = K_bc{2}; E = E{2}; 
boundaryElement = boundaryElement{3};

% fixed boundary condition
K_bc(rols,:) = zeros(size(K_bc(rols,:)));
K_bc(rols,rols) = eye(length(rols));
f(rols,1) = 0;

% right boundary x-direction distributed load
f = nodalForce(f,meshObj,boundaryElement,q);

% solve K_bc*u=f get deformation\stress\strain\nodel_force
u = K_bc\f;
% deformation
u_x_node = full(u(1:2:end-1));
u_y_node = full(u(2:2:end));

spmd
    switch labindex
        case 1
            Nodal_force = K*u;
            Nodal_force_x_node = full(Nodal_force(1:2:end-1));
            Nodal_force_y_node = full(Nodal_force(2:2:end));
        case 2
            [stress,strain] = compute_stress_strain(elementObj,meshObj,u,E);
        case 3
            interp_Drawing(meshObj.x,u_x_node,elementObj.gridMethod,'Directional-deformation-x-direction')
        case 4
            writeReport('Report.txt','node displacement in x direction',u_x_node);
            writeReport('Report.txt','node displacement in y direction',u_y_node);
    end
end
stress = stress{2}; strain = strain{2};

spmd
    switch labindex
        case 1
            interp_Drawing(meshObj.x,stress(:,1),elementObj.gridMethod,'stress-x-nominal')
        case 2
            interp_Drawing(meshObj.x,stress(:,2),elementObj.gridMethod,'stress-y-nominal')
        case 3
            interp_Drawing(meshObj.x,stress(:,3),elementObj.gridMethod,'stress-xy-shear')
        case 4
            writeReport('Report.txt','strain',strain);
            writeReport('Report.txt','stress',stress);
    end
end
