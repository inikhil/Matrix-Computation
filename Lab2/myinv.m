function X = myinv(A)
    [n, n] = size(A);
    [L, U, P] = gepp(A);
    I = eye(n);
    X = zeros(n);
    for i=1:n
        X(:, i) = q4_helper(A, I(:, i), L, U, P);
    end
end

function x = q4_helper(A, b, L, U, P)
    [n, n] = size(A);
    b = P*b;
    y = forward_col_lower(L, b);
    x = backward_col_upper(U, y); 
end