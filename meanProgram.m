%% Plane FEM code
% A finite element method (FEM) code in Fortran for solving general 2-D 
% elastostatic problems using the T3\T6\Q4\Q8 element
% Author: Liu Yazhuo
% Date: 2020.4.30
clear all; clc
young = 210e9;
enu = 0.3;
elementType = 'T6';
meshsz = 1;
q = 1e7;
state = 'Plane Strain';
INTorder = 12;
%% material and mesh parameter
% Variables:
% m       = total number of elements
% n       = total number of nodes
% x       = coordinates of nodes
% node    = element connectivity table
% bc      = boundary conditions
% u       = nodal displacement vector
% stress  = nodal stress vector
% strain  = nodal strain vector
% young   = Young's modulus
% enu     = Poisson's ratio
% meshsz  = maximum mesh size
%% get mesh information
[x,node,gridMethod,elementOrder] = generate2DMesh(elementType,meshsz);
m = size(node,2);
n = size(x,2);

%% Define boundary conditions
left_x = find(x(1,:) == -1)'; % node numbers on left boundary
left_y = left_x;
load_u_right = find(x(1,:) == 1)'; % node numbers on right boundary
f = sparse(2*length(x),1);

%% Compute the stiffness matrix K
[K,E] = stiffness_2D(young,enu,m,n,state,x,node,elementType,INTorder);
K_bc = K;
%% apply deformation boundary conditions
% fixed boundary condition
rols = zeros(length(left_x),1);
rols=2*left_x -1;


% rols = zeros(length(left_x)*2,1);
% rols(1:2:end-1,1)=2*left_x-1;
% rols(2:2:end,1)=2*left_x;

K_bc(rols,:) = zeros(size(K_bc(rols,:)));
K_bc(rols,rols) = eye(length(rols));
f(rols,1) = 0;

% deformation boundary condition
% x-direction
location3 = 2 * load_u_right - 1;
% y-direction
% location3 = 2 * load_u_right;

% x-direction distributed load
% boundaryElement = findBoundaryElement(x,node,elementOrder);
% f = nodalForce(f,x,node,boundaryElement,elementType,q);


K_bc(location3,:) = 0;
K_bc(location3,location3) = eye(length(location3));
f(location3,1)=0.001;

%% solve K_bc*u=f get deformation\stress\strain\nodel_force
u = K_bc\f;
% deformation
u_x_node = full(u(1:2:end-1));
u_y_node = full(u(2:2:end));
% nodel force
Nodel_force = K*u;
Nodel_force_x_node = full(Nodel_force(1:2:end-1));
Nodel_force_y_node = full(Nodel_force(2:2:end));
% stress\strain
[stress,strain] = compute_stress_strain(m,n,x,node,elementType,u,E);

%% Postprocesses: data visualize
interp_Drawing(x,u_x_node,gridMethod,'Directional deformation: x-direction')
interp_Drawing(x,stress(:,1),gridMethod,'stress: x-direction')




