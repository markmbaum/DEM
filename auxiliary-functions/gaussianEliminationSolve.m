function [x] = gaussianEliminationSolve(A,b)
%perform Gaussian Elimination with NO PIVOTING and backward substitution
%to solve a matrix equation

n = length(b);

%gaussian elimination
for i = 1:n-1
   for j = i+1:n
       m = A(j,i)/A(i,i);
       for k = i:n
           A(j,k) = A(j,k) - m*A(i,k);
       end
       b(j) = b(j) - m*b(i);
   end
end

%back-substitution
x = zeros(n,1);
for i = n:-1:1
   s = 0;
   for j = i+1:n
       s = s + A(i,j)*x(j);
   end
   x(i) = (b(i) - s)/A(i,i);
end

end