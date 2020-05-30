function boundaryElement = findBoundaryElement(x,node,elementOrder)
m = size(node,2);
n = size(x,2);

% node at right boundary
rightBoundaryNode = find(x(1,:) == 1);
coordinate = x(:,rightBoundaryNode);
nodeInfo = [rightBoundaryNode;coordinate]';
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