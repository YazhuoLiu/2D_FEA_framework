classdef materials
    properties
        young
        enu
    end
    
    methods
        function materialObj = materials(young,enu)
            materialObj.young = young;
            materialObj.enu = enu;
        end
    end
end