function [c, error]=paulinomials(H)

sigma{1}=[0 1; 1 0];
sigma{2}=[0 -1i; 1i 0];
sigma{3}=[1 0 ; 0 -1];
sigma{4}=eye(2);

m{1}=[0 1 0;1 0 0;0 0 0];
m{2}=[0 -1i 0;1i 0 0;0 0 0];
m{3}=[1  0 0;0 -1 0;0 0 0];
m{4}=[0 0 1;0 0 0;1 0 0];
m{5}=[0  0 -1i;0  0  0; 1i 0  0];
m{6}=[0 0 0;0 0 1;0 1 0];
m{7}=[0  0  0;0  0 -1i;0 1i  0];
m{8}=[1 0 0;0 1 0;0 0 -2]/sqrt(3);
m{9}=eye(3);

switch length(H)
    case 2; for i=1:4; c(i)=(1/2)*trace(sigma{i}*H);end
    case 3
        for i=1:length(H)^2
            basis(:,i)=m{i}(:);
        end
        c=basis\H(:);
    case 4
        for i=1:length(H)
            for j=1:length(H)
                c(i,j)=(1/4)*trace(kron(sigma{i},sigma{j})*H);
            end
        end
    case 6
        for i=1:4
            for j=1:9
                temp=kron(sigma{i},m{j});
                basis(:,sub2ind([4,9],i,j))=temp(:);
            end
        end
        c=basis\H(:);
    case 9
        for i=1:length(H)
            for j=1:length(H)
                temp=kron(m{i},m{j});
                basis(:,sub2ind([length(m),length(m)],i,j))=temp(:);
            end
        end
        c=basis\H(:);
end

error=0;
if length(H)==6
    t=zeros(length(H));
    for i=1:4
        for j=1:9
            t=t+c(sub2ind([4,9],i,j))*kron(sigma{i},m{j});
            error=H-t;
        end
    end
end