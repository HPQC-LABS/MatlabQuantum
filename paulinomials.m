function [a, error]=paulinomials(H)
%%
m{1}=[0 1 0;...
      1 0 0;...
      0 0 0];

m{2}=[0 -1i 0;...
      1i 0 0;...
      0 0 0];

m{3}=[1  0 0;...
      0 -1 0;...
      0 0 0];

m{4}=[0 0 1;...
      0 0 0;...
      1 0 0];

m{5}=[0  0 -1i;...
      0  0  0;...
      1i 0  0];

m{6}=[0 0 0;...
      0 0 1;...
      0 1 0];

m{7}=[0  0  0;...
      0  0 -1i;...
      0 1i  0];

m{8}=[1 0 0;...
      0 1 0;...
      0 0 -2]/sqrt(3);

m{9}=eye(3);

for i=1:length(m)
    for j=1:length(m)
        temp=kron(m{i},m{j});
        basis(:,sub2ind([length(m),length(m)],i,j))=temp(:);
        %c(i,j)=(1/4)*trace(kron(m{i},m{j})*H); % this doesn't work because the basis is not orthogonal.
    end
end

a=basis\H(:);

b=zeros(9);d=b;
for i=1:9
    for j=1:9
        b=b+a(sub2ind([length(m),length(m)],i,j))*kron(m{i},m{j});
%        d=d+c(i,j)*kron(m{i},m{j});
    end
end

error=H-b;
%disp(H-d)
