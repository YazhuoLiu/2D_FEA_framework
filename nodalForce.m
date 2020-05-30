function f = nodalForce(f,meshObj,boundaryElement,q)
x = meshObj.x;
node = meshObj.node;
elementType = meshObj.elementType;

for i = 1:size(boundaryElement,1)
    ele = boundaryElement(i,1);
    isoeleNode = boundaryElement(i,2:end);
    ele_Node = node(isoeleNode,ele);
    xy = x(:,ele_Node);
    J = sqrt(sum((xy(:,1)- xy(:,end)).^2));
    
    rols = ele_Node*2-1;
    if strcmp(elementType,'T3')
        f(rols) = f(rols) + J*[q/2,q/2]';
    elseif strcmp(elementType,'Q4')
        f(rols) = f(rols) + J*[q/2,q/2]';
    elseif strcmp(elementType,'T6')
        f(rols) = f(rols) + J*[q/6,2*q/3,q/6]';
    elseif strcmp(elementType,'Q8')
        f(rols) = f(rols) + J*[q/6,2*q/3,q/6]';
    end
end

end
