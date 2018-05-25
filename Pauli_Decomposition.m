function a=Pauli_Decomposition(H)

sigma{1}=[0 1; 1 0];
sigma{2}=[0 -1i; 1i 0];
sigma{3}=[1 0 ; 0 -1];
sigma{4}=eye(2);

for i=1:4
    for j=1:4
        a(i,j)=(1/4)*trace(kron(sigma{i},sigma{j})*H);
    end
end

