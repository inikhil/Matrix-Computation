function x= gecpsolve(A,b)
[L,U,p,q]=gecp(A);
b=p*b;
z=forward_col_lower(L, b);
y=backward_col_upper(U, z);
x=q*y;
end