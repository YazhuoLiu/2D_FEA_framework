function [x,node,gridMethod,elementOrder] = generate2DMesh(elementType,meshsz)

if strcmp(elementType(1),'T')
    model = createpde(2);
    
    % 2*2 plate
    R1 = [3 4 -1 -1 1 1 -1 1 1 -1]';
    gm = [R1];
    ns = char('R1');
    ns = ns';
    sf = 'R1';
    
    
    % 2*2 squire with a hole
%     R1 = [3 4 -1 -1 1 1 -1 1 1 -1]';
%     C1 = [1 0 0 0.25]';
%     C1 = [C1;zeros(length(R1) - length(C1),1)];
%     gm = [R1,C1];
%     ns = char('R1','C1');
%     ns = ns';
%     sf = 'R1-C1';
    
    g = decsg(gm,sf,ns);
    geometryFromEdges(model,g);
    
    
    % selecte element order
    if strcmp(elementType(2),'3')
        elementORDER = 'linear';
        gridMethod = 'linear';
        elementOrder = 1;
    elseif strcmp(elementType(2),'6')
        elementORDER = 'quadratic';
        gridMethod = 'natural';
        elementOrder = 2;
    else
        fprintf('Invalid element type, program terminated')
        quit cancel;
    end
    meshInfo = generateMesh(model,'Hmax',meshsz,'GeometricOrder',elementORDER);
    x = meshInfo.Nodes;
    node = meshInfo.Elements;
    pdeplot(model);

elseif strcmp(elementType(1),'Q')
    
    if strcmp(elementType(2),'4')
        Q4_NODE = load('Q4_NODE.txt');
        x = Q4_NODE(:,[2,3])';
        Q4_ELEMENT = load('Q4_ELEMENT.txt');
        node = Q4_ELEMENT(:,2:5)';
        gridMethod = 'linear';
        elementOrder = 1;
    elseif strcmp(elementType(2),'8')
        Q8_NODE = load('Q8_NODE.txt');
        x = Q8_NODE(:,[2,3])';
        Q8_ELEMENT = load('Q8_ELEMENT.txt');
        node = Q8_ELEMENT(:,2:9)';
        gridMethod = 'natural';
        elementOrder = 2;
    else
        fprintf('Invalid element type, program terminated')
    end

else
    fprintf('Invalid element type, program terminated')
end
end