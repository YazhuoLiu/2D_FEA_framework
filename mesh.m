classdef mesh
    properties
        x
        node
        elementType
        elementOrder
        meshsz
        m
        n
    end
    
    methods 
        function meshObj = mesh(elementObj,meshsz)
            meshObj.meshsz = meshsz;
            meshObj.elementType = elementObj.elementType;
            meshObj.elementOrder = elementObj.elementOrder;
            elementType = elementObj.elementType;
            if strcmp(elementType(1),'T')
                model = createpde(2);
                % 2*2 plate
                R1 = [3 4 -1 -1 1 1 -1 1 1 -1]';
                gm = [R1];
                ns = char('R1');
                ns = ns';
                sf = 'R1';
                % 2*2 squire with a hole
%               R1 = [3 4 -1 -1 1 1 -1 1 1 -1]';
%               C1 = [1 0 0 0.25]';
%               C1 = [C1;zeros(length(R1) - length(C1),1)];
%               gm = [R1,C1];
%               ns = char('R1','C1');
%               ns = ns';
%               sf = 'R1-C1';
                g = decsg(gm,sf,ns);
                geometryFromEdges(model,g);
                % selecte element order
                if strcmp(elementType(2),'3')
                    elementORDER = 'linear';
                elseif strcmp(elementType(2),'6')
                    elementORDER = 'quadratic';
                else
                    fprintf('Invalid element type, program terminated')
                end
                meshInfo = generateMesh(model,'Hmax',meshsz,'GeometricOrder',elementORDER);
                meshObj.x = meshInfo.Nodes;
                meshObj.node = meshInfo.Elements;
                pdeplot(model);

            elseif strcmp(elementType(1),'Q')
    
                if strcmp(elementType(2),'4')
                    Q4_NODE = load('Q4_NODE.txt');
                    meshObj.x = Q4_NODE(:,[2,3])';
                    Q4_ELEMENT = load('Q4_ELEMENT.txt');
                    meshObj.node = Q4_ELEMENT(:,2:5)';
                    
                elseif strcmp(elementType(2),'8')
                    Q8_NODE = load('Q8_NODE.txt');
                    meshObj.x = Q8_NODE(:,[2,3])';
                    Q8_ELEMENT = load('Q8_ELEMENT.txt');
                    meshObj.node = Q8_ELEMENT(:,2:9)';
                    
                else
                    fprintf('Invalid element type, program terminated')
                end
            end
            meshObj.m = size(meshObj.node,2);
            meshObj.n = size(meshObj.x,2);
        end
    end
    
    methods(Static)
        
        function boundaryElement = BoundaryElement(meshObj)
            node = meshObj.node;
            x = meshObj.x;
            elementOrder = meshObj.elementOrder;

            % node at right boundary
            BoundaryNode = find(x(1,:) == 1);
            coordinate = x(:,BoundaryNode);
            nodeInfo = [BoundaryNode;coordinate]';
            nodeInfo = sortrows(nodeInfo,3);
            nodeInfo = nodeInfo(:,1);

            if elementOrder == 1
                boundaryElement = zeros(length(nodeInfo)-1,3);
                for i = 1:length(nodeInfo)-1
                    [row,col] = find(node==nodeInfo(i));
                % position = [row,col];
                    for j = 1:length(col)
                        k = find(node(:,col(j))==nodeInfo(i+1));
                        if isempty(k)
                            continue
                        else
                            elementNumber = col(j);
                            IsoNodeNumber = [row(j),k];
                            boundaryElement(i,1)= elementNumber;
                            boundaryElement(i,[2,3]) = IsoNodeNumber;
                        end
                    end 
                end
            elseif elementOrder == 2
                boundaryElement = zeros((length(nodeInfo)-1)/2,4);
                for i = 1:2:length(nodeInfo)-1
                    [row,col] = find(node==nodeInfo(i));
                % position = [row,col];
                    for j = 1:length(col)
                        k = find(node(:,col(j))==nodeInfo(i+1));
                        if isempty(k)
                            continue
                        else
                            elementNumber = col(j);
                            temp = find(node(:,col(j))==nodeInfo(i+2));
                            IsoNodeNumber = [row(j),k,temp];
                            boundaryElement((i+1)/2,1)= elementNumber;
                            boundaryElement((i+1)/2,[2,3,4]) = IsoNodeNumber;
                        end
                    end 
                end
            end
        end
        
        
        
    end
end