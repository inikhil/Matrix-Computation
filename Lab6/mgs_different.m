function [Q, R] = mgs_different(Q)
    [n, m] = size(Q);
    R = zeros(m, m);
    for k=1:m
        R(k,k)= norm(Q(:,k),2);
        if(R(k,k) == 0)
            error('Condensed QR not possible!')
        end
        Q(:,k)=Q(:,k)/R(k,k);
        for j=k+1:m
            R(k,j)=Q(:,j)'*Q(:,k);
        end
        for j=k+1:m
            Q(:,j)=Q(:,j)-R(k,j)*Q(:,k);
        end
    end
end
    