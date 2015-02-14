C = [];
N = 2:2:16;
for n=N
    H = hilb(n);
    C = [C; cond(H)];
end
figure();
semilogy(N, C);
C = [];
for n=N
    H = hilb(n);
    C = [C; cond(H, 1)]; %cond(X,p) returns the matrix condition number in p-norm:%
end
figure();
semilogy(N, C);
C = [];
for n=N
    H = hilb(n);
    C = [C; cond(H, inf)];
end
figure();
semilogy(N, C)