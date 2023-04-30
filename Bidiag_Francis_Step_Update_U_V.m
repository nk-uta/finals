function [ U, Bi, V ] = Bidiag_Francis_Step_Update_U_V(U, Bi, V)

m = size(Bi, 1);
tau00 = Bi(1, 1) ^ 2;
tau10 = Bi(1, 2) * Bi(1, 1);
tau_m1_m1 = Bi(m-1, m) ^ 2 + Bi(m, m) ^ 2;

column = [tau00 - tau_m1_m1
		tau10];
G = Givens_rotation( column);
Bi(:, 1:2) = Bi(:, 1:2) * G;

% Ideally this update is supposed to be VT = GT * VT. But realize that VT = GT * VT is same as V = V * G
% This is like VT is being hit with Givens rotation from the left to converge to U * Sigma * V'
V(:, 1:2) = V(:, 1:2) * G;

for i = 1:size(Bi, 1) - 1
		row = [Bi(i, i),
					Bi(i + 1, i)];
		G_hat = Givens_rotation( row);
		Bi(i:i+1,:) = G_hat' * Bi(i:i+1,:);	
		U(:, i:i+1) = U(:, i:i+1) * G_hat;
		
		if i < m - 1
		column = [Bi(i, i + 1)
					Bi(i, i + 2)];
		G_tilde = Givens_rotation( column);
		Bi(:, i+1:i+2) = Bi(:, i+1:i+2) * G_tilde;
		V(:,i+1:i+2) = V(:, i+1:i+2) * G_tilde;
		end
end	

Bi;
U;
V;

end