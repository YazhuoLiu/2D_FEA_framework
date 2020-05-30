function [] = writeReport(fileName,tittle,data)
    fid=fopen(fileName,'a'); 
    fprintf(fid,tittle);
    fprintf(fid,'\n');
    for i = 1:size(data,1)
        for j = 1:size(data,2)
            fprintf(fid,'%f\t',data(i,j));
        end
        fprintf(fid,'\n\n');
    end
    fclose(fid);
end
