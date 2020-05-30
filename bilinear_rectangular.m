classdef bilinear_rectangular
    properties
        elementType
        elementOrder
        gridMethod
        INTorder
        Node_xi = [-1;1;1;-1];
        Node_eta = [-1;-1;1;1];
    end
   
    methods
        
        function elementObj = bilinear_rectangular(INTorder)
            elementObj.elementType = 'Q4';
            elementObj.elementOrder = 1;
            elementObj.gridMethod = 'linear';
            elementObj.INTorder = INTorder;
        end
        
        function [s,wt] = GQintegration(elementObj)
            Order = elementObj.INTorder;
            if Order == 1
                s(1,1) = .0 ; wt(1) = 4;
            elseif Order == 4
                s(1,1)= -1/sqrt(3);s(1,2)= -1/sqrt(3);
                s(2,1)= 1/sqrt(3); s(2,2)= -1/sqrt(3);
                s(3,1)= 1/sqrt(3); s(3,2)= 1/sqrt(3);
                s(4,1)= -1/sqrt(3); s(4,2)= 1/sqrt(3);
                wt = ones(Order,1);
            elseif Order == 9
                s(1,1)= -sqrt(0.6);s(1,2)= -sqrt(0.6); wt(1) = 25/81;
                s(2,1)= -sqrt(0.6); s(2,2)= 0; wt(2) = 40/81;
                s(3,1)= -sqrt(0.6); s(3,2)= sqrt(0.6); wt(3) = 25/81;
                s(4,1)= 0; s(4,2)= -sqrt(0.6); wt(4) = 40/81;
                s(5,1)= 0; s(5,2)= 0; wt(5) = 64/81;
                s(6,1)= 0; s(6,2)=sqrt(0.6); wt(6) = 40/81;
                s(7,1)= sqrt(0.6); s(7,2) = -sqrt(0.6); wt(7) = 25/81;
                s(8,1)= sqrt(0.6); s(8,2)= 0; wt(8) = 40/81;
                s(9,1)= sqrt(0.6); s(9,2)= sqrt(0.6); wt(9) = 25/81;
            else
                fprintf('wrong number of integrating points for a quadrilateral.')
            end
        end
        
        function J = Jacobian(elementObj,meshObj,NoEle,xi,eta)
            nodecoordinates = meshObj.x(:,meshObj.node(:,NoEle));
            Nxi  = 0.25 * [eta-1 , 1-eta , 1+eta , -1-eta ];
            Neta = 0.25 * [xi-1  , -1-xi , 1+xi  , 1-xi   ];
            J = [Nxi;Neta] * nodecoordinates';
        end
        
        function B = StrainMatrix(elementObj,meshObj,NoEle,xi,eta)
            num = [1 0 0 0 ; 0 0 0 1 ; 0 1 1 0];
            nodecoordinates = meshObj.x(:,meshObj.node(:,NoEle));
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
        end
        
    end
end