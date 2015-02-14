function[]=q0()
    load clown.mat;
    [U, S, V] = svd(X);
    colormap('gray');
    for k=5:5:195
        image(U(:, 1:k)*S(1:k, 1:k)*V(:,1:k)');
        S(k+1,k+1)/S(1,1)
    end
end