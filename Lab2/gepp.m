function [L, U, P] = gepp(A)
    [n, n] = size(A);
    p = [1:n];
    for k=1:n-1
        [r, m] = max(abs(A(k:n, k)));
        m = m + k -1;
        if(A(m , k)==0)
            sprintf('A is singular');
        else if(m~=k)
            A([k m], :) = A([m k], :);
            p([k m]) = p([m k]);
            end
        end
        A(k+1:n, k) = A(k+1:n, k)/A(k, k);
        A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - A(k+1:n, k)*A(k, k+1:n);
    end
    L = tril(A, -1) + eye(n);
    U = triu(A);
    I = eye(n);
    P = I(p, :);
end
























































function [pr,l,u,pc,err]=gecp(a)
% pr*l*u*pc'=a or pr'*a*pc=l*u
n=max(size(a));
pr=eye(n);pc=eye(n);
aa=a;
for i=1:n-1,
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

















function r = cholsky(A)
    %chol(A)
    [n,n]=size(A);
    for i=1:n
       % for k=1:i-1
        %    A(i,i)=A(i,i)-A(k,i)*A(k,i);
        %end
         A(i,i)=A(i,i)-sum(A(1:(i-1),i).*A(1:(i-1),i));
        if(A(i,i)<=0)
            sprintf('A is not positive definite');
        else
            A(i,i)=sqrt(A(i,i));
            for j=i+1:n 
                %for k=1:i-1
                 %   A(i,j)=A(i,j)-A(k,i)*A(k,j);
                %end
                A(i,j)=A(i,j)-sum(A(1:(i-1),i).*A(1:(i-1),j));
                A(i,j)=A(i,j)/A(i,i);
            end
        end
    end
    r=triu(A);
end


function [L,U,P,Q]=gecp(A)
% pr*l*u*pc'=a or pr'*a*pc=l*u
[n,n]=size(A);
p=[1:n];
q=[1:n];
for i=1:n-1,
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
      A(:,[i jmax]) = A(:,[jmax i]);
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