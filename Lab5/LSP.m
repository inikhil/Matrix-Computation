function [x] = LSP(t)

format long e;

m = length(t);
n = 20;

A = zeros(m,n);
count = 0;

for i=1:n
    A(:,i) = t.^count;
    count = count+1;
end


b = sin(pi*t/5) + (t/5);

disp(cond(A,2));

x = A \ b';

x = flipud(x);

err = polyval(x,t) - b;
disp(norm(err,2));

end
