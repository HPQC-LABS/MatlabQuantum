function a=Pauli_Decomposition(H)

sigma{1}=[0 1; 1 0];
sigma{2}=[0 -1i; 1i 0];
sigma{3}=[1 0 ; 0 -1];
sigma{4}=eye(2);

% for i=1:4
%     for j=1:4
%         a(i,j)=(1/4)*trace(kron(sigma{i},sigma{j})*H);
%     end
% end

%% alternative more general way:

for i=1:4
    for j=1:4
        temp=kron(sigma{i},sigma{j});
        basis(:,sub2ind([4,4],i,j))=temp(:);
    end
end

a=basis\H(:);

b=zeros(4);
for i=1:4
    for j=1:4
        b=b+a(sub2ind([4,4],i,j))*kron(sigma{i},sigma{j});
    end
end

disp(H-b)
