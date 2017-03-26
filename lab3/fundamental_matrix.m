function F_denormalised = fundamental_matrix(x1, x2)

% normalise
[x1_normalised, T1] = normalise2dpts(x1);
[x2_normalised, T2] = normalise2dpts(x2);

%total points
n = size(x1_normalised, 2); % n is 8

W = zeros(n,9);
for i=1:n
    u = x1_normalised(1,i);
    u_p = x2_normalised(1,i);
    v = x1_normalised(2,i);
    v_p = x2_normalised(2,i);
    
    W(i, :) = [u*u_p, v*u_p, u_p, u*v_p, v*v_p, v_p, u, v, 1];
end

[U, D, V] = svd(W);

f_rank3 = V(:,9);
% f_rank3 = reshape(f_rank3, [3,3]);  % TODO: review why reshape gives the transposed matrix ¿?
f_rank3 = [f_rank3(1) f_rank3(2) f_rank3(3);
           f_rank3(4) f_rank3(5) f_rank3(6);
           f_rank3(7) f_rank3(8) f_rank3(9)];

[U, D, V] = svd(f_rank3);

D_hat = D;
D_hat(3,3) =  0;

F = U * D_hat * V';

F_denormalised = T2' * F * T1; 