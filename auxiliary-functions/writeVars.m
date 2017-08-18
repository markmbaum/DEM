function writeVars(ofile, x, y, theta)

for i = 1:length(x)
    fprintf(ofile, '%f,%f,%f\n', x(i), y(i), theta(i));
end

end
