clear
clc

load('data_M.mat')
load('zero_point');
load('number.mat');

b=0.1;
c(1)=0.1;
c(2)=0.2;
%c(3)=0.2;

delta(1)=pi/2;
delta(2)=pi/2;

C(1)=2*b;
C(2)=3*c(1)*exp(sqrt(-1)*delta(1));
C(3)=4*c(2)*exp(sqrt(-1)*delta(2));
l=10^3;

A=ones(l);
J=zeros(l,l);
M1_1=M1_1_1.*tril(A,-1)*C(1)+M1_1_1.* triu(A,1)*conj(C(1))+M1_1_2.*tril(A,-1)*conj(C(1))+M1_1_2.* triu(A,1)*C(1);
M3_1=M3_1_1.*tril(A,-1)*C(2)*conj(C(1))+M3_1_1.* triu(A,1)*C(1)*conj(C(2))+M3_1_2.*tril(A,-1)*C(1)*conj(C(2))+M3_1_2.* triu(A,1)*C(2)*conj(C(1));
M5_1=M5_1_1.*tril(A,-1)*C(3)*conj(C(2))+M5_1_1.* triu(A,1)*C(2)*conj(C(3))+M5_1_2.*tril(A,-1)*C(2)*conj(C(3))+M5_1_2.* triu(A,1)*C(3)*conj(C(2));

M2_2=M2_2_1.*tril(A,-1)*C(2)+M2_2_1.* triu(A,1)*conj(C(2))+M2_2_2.*tril(A,-1)*conj(C(2))+M2_2_2.* triu(A,1)*C(2);
M4_2=M4_2_1.*tril(A,-1)*C(3)*conj(C(1))+M4_2_1.* triu(A,1)*C(1)*conj(C(3))+M4_2_2.*tril(A,-1)*C(1)*conj(C(3))+M4_2_2.* triu(A,1)*C(3)*conj(C(1));

M3_3=M3_3_1.*tril(A,-1)*C(3)+M3_3_1.* triu(A,1)*conj(C(3))+M3_3_2.*tril(A,-1)*conj(C(3))+M3_3_2.* triu(A,1)*C(3);

J=1/(1+2*b^2+3*c(1)^2+4*c(2)^2).*(eye(l)+4*b^2*M2_0+9*c(1)^2*M4_0+16*c(2)^2*M6_0+M1_1+M3_1+M5_1+M2_2+M4_2+M3_3);

for i=1:l
    U(i,i)=1/zero_point(i);
end
J1=U*J*U;
disp(1)
[D T]=eig(J1);
disp(2) 
for i=1:l
    eigv_statical(i)=sqrt(1/real(T(i,i)));
    eigs_statical(:,i)=U*D(:,i);
end
%eigv_statical=sort(eigv_statical);
save([pwd,'/eigs_statical.mat'],'eigs_statical');
save([pwd,'/eigv_statical.mat'],'eigv_statical');