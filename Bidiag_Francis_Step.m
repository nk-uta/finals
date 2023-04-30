function Bi = Bidiag_Francis_Step(Bi)

m = size(Bi, 1);
tau00 = Bi(1, 1) ^ 2;
tau10 = Bi(1, 2) * Bi(1, 1);
tau_m1_m1 = Bi(m-1, m) ^ 2 + Bi(m, m) ^ 2;

column = [tau00 - tau_m1_m1
		tau10];
G = Givens_rotation( column);
G' * column;
Bi(1:2, 1:2) = Bi(1:2, 1:2) * G;

for i = 1:size(Bi, 1) - 1
		row = [Bi(i, i),
					Bi(i + 1, i)];
		G_hat = Givens_rotation( row);
		Bi(i:i+1,:) = G_hat' * Bi(i:i+1,:);
		
		if i < m - 1
		column = [Bi(i, i + 1)
					Bi(i, i + 2)];
		G_tilde = Givens_rotation( column);
		Bi(:, i+1:i+2) = Bi(:, i+1:i+2) * G_tilde;
		end
end	
		
%{		
column2 = [Bi(1, 1)
			Bi(2, 1)];
G_hat = Givens_rotation( column2);
Bi(1:2,:) = G_hat' * Bi(1:2,:);


column3 = [Bi(1, 2)
			Bi(1, 3)];
G_tilde = Givens_rotation( column3);
Bi(:, 2:3) = Bi(:, 2:3) * G_tilde;


column4 = [Bi(2, 2)
			Bi(3, 2)];
G_hat = Givens_rotation( column4);
Bi(2:3,:) = G_hat' * Bi(2:3,:);


column5 = [Bi(2, 3)
			Bi(2, 4)];
G_tilde = Givens_rotation( column5);
Bi(:, 3:4) = Bi(:, 3:4) * G_tilde;


column6 = [Bi(3, 3)
			Bi(4, 3)];
G_hat = Givens_rotation( column6);
Bi(3:4,:) = G_hat' * Bi(3:4,:);


column7 = [Bi(3, 4)
			Bi(3, 5)];
G_tilde = Givens_rotation( column7);
Bi(:, 4:5) = Bi(:, 4:5) * G_tilde;

column8 = [Bi(4, 4)
			Bi(5, 4)];
G_hat = Givens_rotation( column8);
Bi(4:5,:) = G_hat' * Bi(4:5,:);
%}

Bi;

end