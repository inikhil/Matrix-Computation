function [W] = Wilkinson(n)
        W = eye(n); %only diagonal entries 1%
        W(find(W==0)) = -1;
        W = tril(W);
        W(:,n) = 1;
end

