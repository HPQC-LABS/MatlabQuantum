b0=[1 0];
b1=[0 1];
t0=[1 0 0];
t1=[0 1 0];
t2=[0 0 1];

t1b1=kron(t1,b1);
t2b1=kron(t2,b1);

b0t1=kron(b0,t1);
b1t0=kron(b1,t0);
b1t1=kron(b1,t1);
b1t2=kron(b1,t2);
b0t2=kron(b0,t2);

sigma{1}=[0 1; 1 0];
sigma{2}=[0 -1i; 1i 0];
sigma{3}=[1 0 ; 0 -1];
sigma{4}=eye(2);

x=sigma{1};
y=sigma{2};
z=sigma{3};

x1=kron(x,eye(2));
x2=kron(eye(2),x);
y1=kron(y,eye(2));
y2=kron(eye(2),y);
z1=kron(z,eye(2));
z2=kron(eye(2),z);

m{1}=[0 1 0;1 0 0;0 0 0];
m{2}=[0 -1i 0;1i 0 0;0 0 0];
m{3}=[1  0 0;0 -1 0;0 0 0];
m{4}=[0 0 1;0 0 0;1 0 0];
m{5}=[0  0 -1i;0  0  0; 1i 0  0];
m{6}=[0 0 0;0 0 1;0 1 0];
m{7}=[0  0  0;0  0 -1i;0 1i  0];
m{8}=[1 0 0;0 1 0;0 0 -2]/sqrt(3);
m{9}=eye(3);
