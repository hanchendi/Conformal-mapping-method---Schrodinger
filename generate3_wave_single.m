clear
clc

load('eigs_statical.mat');
load('zero_point.mat');
load('normalized.mat');
load('number.mat');
mode=5;
M=100;
N=300;

b=0.1;
c(1)=0.1;
c(2)=0.2;
delta(1)=pi/2;
delta(2)=pi/2;
l=1000;
rr=linspace(0,1,N);
rr=rr';
delta_theta=2*pi/N;
theta=delta_theta:delta_theta:2*pi;

xx1=rr*cos(theta);
yy1=rr*sin(theta);
z=xx1+sqrt(-1)*yy1;
w=(z+b*z.^2+c(1)*exp(sqrt(-1)*delta(1))*z.^3+c(2)*exp(sqrt(-1)*delta(2))*z.^4)/sqrt(1+2*b^2+3*c(1)^2+4*c(2)^2);
xx2=real(w);
yy2=imag(w);
psi=zeros(N,N);

A(1:l)=abs(eigs_statical(1:l,mode));
t=1;
for i=1:M
    aa=find(A==max(A));
    if length(aa)~=1
        for j=1:length(aa)
            a(t)=aa(j);
            A(a(t))=0;
            t=t+1;
        end
    else
        a(t)=aa;
        A(a(t))=0;
        t=t+1;
    end
end

for i=1:t-1
    psi=psi+eigs_statical(a(i),mode)*normalized(a(i))*(besselj(number(a(i),1),rr.*zero_point(a(i)))*exp(sqrt(-1)*number(a(i),1)*(theta)));
end

xx2(:,N+1)=xx2(:,1);
yy2(:,N+1)=yy2(:,1);
psi(:,N+1)=psi(:,1);
figure(1)
pcolor(xx2,yy2,abs(psi).^2);
shading flat
axis off
axis equal
phi_C=abs(psi).^2;

save([pwd,'/Conformal_',num2str(mode),'.mat'],'phi_C');