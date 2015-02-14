function [pr,l,u,pc,err]=gexp(a)
% pr*l*u*pc'=a or pr'*a*pc=l*u
n=max(size(a));
pr=eye(n);pc=eye(n);
aa=a;
for i=1:n-1
%  find largest absolute entry in (i:n,i:n) submatrix
   am=max(max(abs(aa(i:n,i:n))));
   [I,J]=find(abs(aa(i:n,i:n)) == am);
   imax=I(1)+i-1;
   jmax=J(1)+i-1;
%  swap rows
   if (imax ~= i),
      temp = pr(:,i);
      pr(:,i) = pr(:,imax);
      pr(:,imax) = temp;
      temp = aa(i,:);
      aa(i,:) = aa(imax,:);
      aa(imax,:) = temp;
   end
%  swap columns
   if (jmax ~= i),
      temp = pc(:,i);
      pc(:,i) = pc(:,jmax);
      pc(:,jmax) = temp;
      temp = aa(:,i);
      aa(:,i) = aa(:,jmax);
      aa(:,jmax) = temp;
   end
%  eliminate
   aa(i+1:n,i) = aa(i+1:n,i)/aa(i,i);
   aa(i+1:n,i+1:n) = aa(i+1:n,i+1:n) - aa(i+1:n,i)*aa(i,i+1:n);
end
l=eye(n);l=l+tril(aa,-1);
u=triu(aa,0);
err=[norm(pr*l*u*pc'-a,1)/norm(a,1); cond(a); cond(l); cond(u)];
end












