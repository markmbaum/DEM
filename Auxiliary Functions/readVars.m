function [r, x, y, theta] = readVars(file_name)

%open file and get radii
ifile = fopen(file_name, 'r');
l = fgetl(ifile);
l = fgetl(ifile);
l = strsplit(l, ',');
r = zeros(length(l), 1);
for i = 1:length(r)
    r(i) = str2num(l{i});
end
N = length(r);
%count the number of lines to get the number of iterations
%read junk lines
for i = 1:3
    l = fgetl(ifile);
end
%count lines of data
count = 0;
l = fgetl(ifile);
while(l ~= -1)
    count = count + 1;
    l = fgetl(ifile);
end
count = count/N;
%close and reopen file
fclose(ifile);
ifile = fopen(file_name, 'r');
%read junk lines
for i = 1:5
    l = fgetl(ifile);
end
%read data
x = zeros(count, N);
y = zeros(count, N);
theta = zeros(count, N);
for i = 1:count
    for j = 1:N
        v = fscanf(ifile, '%f,%f,%f\n', 3);
        x(i,j) = v(1);
        y(i,j) = v(2);
        theta(i,j) = v(3);
    end
end

end
