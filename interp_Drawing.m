function [] = interp_Drawing(x,nodelValue,gridMethod,titles)
    x_coord = x(1,:)';
    y_coord = x(2,:)';
    [xq,yq] = meshgrid(-1:0.05:1,-1:0.05:1);
    
%     for i = 1:size(xq,1)
%        for j = 1:size(xq,2)
%            if xq(i,j)^2+yq(i,j)^2 < 0.25^2
% %                interp_Value(i,j) = NaN;
%                xq(i,j) = NaN;
%                yq(i,j) = NaN;
%            end
%        end
%     end
    
    interp_Value = griddata(x_coord,y_coord,nodelValue,xq,yq,gridMethod);
    figure
    set(gcf,'unit','centimeters','position',[3,3,35,12])
    
    
    
    
    subplot(1,2,1)
    % surfc(xq,yq,interp_Value,'EdgeColor','none')
    surfc(xq,yq,interp_Value)
    % s.Surface.EdgeColor = 'none'
    colorbar
    title([titles,' (space diagram)'])
    xlabel('x');ylabel('y');zlabel(titles);
    
    subplot(1,2,2)
    contourf(xq,yq,interp_Value)
    colorbar
    title([titles,' (ichnography)'])
    xlabel('x');ylabel('y')
    axis equal
end