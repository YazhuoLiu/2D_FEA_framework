classdef linear_strain_triangular
    properties
        elementType
        elementOrder
        gridMethod
        INTorder
        Node_xi = [1;0;0;0.5;0;0.5];
        Node_eta = [0;1;0;0.5;0.5;0];
    end
   
    methods
        
        function elementObj = linear_strain_triangular(INTorder)
            elementObj.elementType = 'T6';
            elementObj.elementOrder = 2;
            elementObj.gridMethod = 'natural';
            elementObj.INTorder = INTorder;
        end
        
        function [s,wt] = GQintegration(elementObj)
            Order = elementObj.INTorder;
            if Order == 6
                s(1,1) = 0.816847572980459 ; s(1,2)= 0.091576213509771;
                s(2,1) = s(1,2); s(2,2) = s(1,1); s(3,1) = s(1,2); s(3,2) = s(1,2);
                s(4,1) = 0.108103018168070 ;  s(4,2)=.445948490915965;
                s(5,1) = s(4,2); s(5,2) = s(4,1); s(6,1) = s(4,2); s(6,2) = s(4,2);
                wt(1)=.109951743655322; wt(2)=wt(1); wt(3)=wt(1);
                wt(4)=.223381589678011; wt(5)=wt(4); wt(6)=wt(4);
                wt = 0.5 * wt;
            elseif Order == 1
                s(1,1)=1./3.  ; s(1,2)=1./3.  ;  wt(1)= .5;
            elseif Order == 3
                s(1,1)=.5 ;  s(1,2)=.5 ;  s(2,1)=.5;
                s(2,2)=0.;  s(3,1)=0.  ;  s(3,2)=.5;
                wt(1)=1./3.  ;  wt(2)=wt(1) ; wt(3)=wt(1); wt = .5*wt;
            elseif Order == 7
                s(1,1)=1./3. ; s(1,2)=1./3.;s(2,1)=.797426985353087 ;s(2,2)=.101286507323456;
                s(3,1)=s(2,2) ;  s(3,2)=s(2,1) ; s(4,1)=s(2,2) ;  s(4,2)=s(2,2);
                s(5,1)=.470142064105115 ;   s(5,2)=.059715871789770;
                s(6,1)=s(5,2) ; s(6,2)=s(5,1);  s(7,1)=s(5,1);  s(7,2)=s(5,1);
                wt(1)=.225 ; wt(2)=.125939180544827 ;  wt(3)=wt(2);  wt(4)=wt(2);
                wt(5)=.132394152788506;  wt(6)=wt(5)      ;  wt(7)=wt(5)     ;wt = .5*wt;
            elseif Order == 12
                s(1,1)=.873821971016996 ; s(1,2)=.063089014491502;
                s(2,1)=s(1,2) ;  s(2,2)=s(1,1);  s(3,1)=s(1,2) ;  s(3,2)=s(1,2);
                s(4,1)=.501426509658179 ;  s(4,2)=.249286745170910;
                s(5,1)=s(4,2); s(5,2)=s(4,1)   ;  s(6,1)=s(4,2) ;  s(6,2)=s(4,2);
                s(7,1)=.636502499121399 ;      s(7,2)=.310352451033785;
                s(8,1)=s(7,1) ;  s(8,2)=.053145049844816 ;  s(9,1)=s(7,2) ; s(9,2)=s(7,1);
                s(10,1)=s(7,2) ; s(10,2)=s(8,2) ; s(11,1)=s(8,2);   s(11,2)=s(7,1);
                s(12,1)=s(8,2) ;  s(12,2)=s(7,2);
                wt(1)=.050844906370207 ; wt(2)=wt(1); wt(3)=wt(1);
                wt(4)=.116786275726379 ; wt(5)=wt(4); wt(6)=wt(4);
                wt(7)=.082851075618374 ; wt(8:12)=wt(7)           ; wt = .5*wt;
            elseif Order == 16
                s(1,1)=1./3. ;  s(1,2)=1./3.  ;  s(2,1)=.658861384496478;
                s(2,2)=.170569307751761 ; s(3,1)=s(2,2)   ;  s(3,2)=s(2,1);
                s(4,1)=s(2,2)  ; s(4,2)=s(2,2);
                s(5,1)=.898905543365938 ; s(5,2)=.050547228317031;
                s(6,1)=s(5,2);  s(6,2)=s(5,1) ; s(7,1)=s(5,2)  ;  s(7,2)=s(5,2);
                s(8,1)=.081414823414554; s(8,2)=.459292588292723;
                s(9,1)=s(8,2)  ;  s(9,2)=s(8,1);  s(10,1)=s(8,2) ;  s(10,2)=s(8,2);
                s(11,1)=.008394777409958; s(11,2)=.263112829634638;
                s(12,1)=s(11,1)    ;  s(12,2)=.728492392955404;
                s(13,1)=s(11,2) ;   s(13,2)=s(11,1)  ;  s(14,1)=s(11,2); s(14,2)=s(12,2);
                s(15,1)=s(12,2) ;  s(15,2)=s(11,1) ;  s(16,1)=s(12,2) ;  s(16,2)=s(11,2);
                wt(1)=.144315607677787 ; wt(2)=.103217370534718 ; wt(3)=wt(2); wt(4)=wt(2);
                wt(5)=.032458497623198 ; wt(6)=wt(5)   ;  wt(7)=wt(5);
                wt(8)=.095091634267284 ; wt(9)=wt(8)   ;  wt(10)=wt(8);
                wt(11)=.027230314174435 ; wt(12:16) = wt(11)  ;     wt = .5*wt;
            else
                fprintf('wrong number of integrating points for a triangle')
            end 
        end
        
        function J = Jacobian(elementObj,meshObj,NoEle,xi,eta)
            nodecoordinates = meshObj.x(:,meshObj.node(:,NoEle));
            Nxi  = [4*xi-1 ,0       , 4*eta+4*xi-3, 4*eta , -4*eta       ,4-8*xi-4*eta];
            Neta = [0      ,4*eta-1 , 4*eta+4*xi-3, 4*xi  , 4-4*xi-8*eta ,-4*xi       ];
            J = [Nxi;Neta] * nodecoordinates';
        end
        
        function B = StrainMatrix(elementObj,meshObj,NoEle,xi,eta)
            num = [1 0 0 0 ; 0 0 0 1 ; 0 1 1 0];
            nodecoordinates = meshObj.x(:,meshObj.node(:,NoEle));
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
        end
        
    end
end