function [L,U,P,Q]=gecp(A)
% pr*l*u*pc'=a or pr'*a*pc=l*u
[n,n]=size(A);
p=[1:n];
q=[1:n];
for i=1:n-1
%  find largest absolute entry in (i:n,i:n) submatrix
   am=max(max(abs(A(i:n,i:n))));
   [I,J]=find(abs(A(i:n,i:n)) == am);
   imax=I(1)+i-1;
   jmax=J(1)+i-1;
%  swap rows
   if (imax ~= i),
      p([i imax])= p([imax i]);
      A([i imax], :) = A([imax i], :);
   end
%  swap columns
   if (jmax ~= i),
      q([i jmax])=q([jmax i]);
      A(:,[i jmax] ) = A(:,[jmax i]);
   end
%  eliminate
   A(i+1:n,i) = A(i+1:n,i)/A(i,i);
   A(i+1:n,i+1:n) = A(i+1:n,i+1:n) - A(i+1:n,i)*A(i,i+1:n);
end
I=eye(n);
L=I+tril(A,-1);
U=triu(A);
P=I(p,:);
Q=I(:,q);
%err=[norm(pr*l*u*pc'-a,1)/norm(a,1); cond(a); cond(l); cond(u)];
end