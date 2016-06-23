function fileName = saveName(fileName, extension)
%matSaveName creates a filename to save data under that will not overwrite
%other files with the same name.
%   filename - path string without an extension
%   extension - file extension, including the period

%finding files with the same name
temp = ls;
L = length(fileName);
l_ext = length(extension);
Z = cell(0);
j = 1;
count = 0;
if(size(temp,1) > 1) %pc
    V = size(temp,2);
    Q = size(temp,1);
    while(j <= Q)
        if(strcmp(temp(j,1:L),fileName))
            count = count+1;
            Z{count} = temp(j,L+1:V);
        end
        j = j + 1;
    end
else %mac
    while(j < (length(temp) - L))
        if(strcmp(temp(j:j+L-1),fileName))
            if(strcmp(temp(j+L:j+L+l_ext-1), extension))
                if(j == 1)
                    count = count + 1;
                    Z{count} = temp(j+L:end);
                elseif((j ~= 1) && isspace(temp(j-1)))
                    count = count + 1;
                    Z{count} = temp(j+L:end);
                end
            end
        end
        j = j + 1;
    end
end
%determing what the numerical suffix of the file needs to be
temp = zeros(1,length(Z));
for j = 1:length(Z)
    k = 1;
    done = 0;
    while((k < (length(Z{j}) - 3)) && ~done)
        if(strcmp(Z{j}(k:k+l_ext), extension))
            done = 1;
            if(~isempty(Z{j}(2:k-1)))
                temp(j) = str2double(Z{j}(2:k-1));
            end
        end
        k = k + 1;
    end
end
%changing filename and saving
temp = max(temp)+1;
if(temp == 1)
    temp = temp+1;
end
if(isempty(temp))
    temp = fileName;
else
    temp = [fileName,num2str(temp)];
end

fileName = [temp, extension];

end
