clear
clc

load('number.mat');
load('zero_point.mat');

l=1000;
normalized=zeros(1,l);
for i=1:l
    b1=number(i,1);
    b2=number(i,2);
    b3=zero_point(i,1);
    if b1>=0
        normalized(i)=sqrt(pi)*(besselj(b1+1,b3));
    else
        normalized(i)=sqrt(pi)*(besselj(b1-1,b3));
    end
    normalized(i)=1/normalized(i);
end

M2_0=zeros(l,l);
M4_0=zeros(l,l);
M6_0=zeros(l,l);

M1_1_1=zeros(l,l);
M1_1_2=zeros(l,l);
M3_1_1=zeros(l,l);
M3_1_2=zeros(l,l);
M5_1_1=zeros(l,l);
M5_1_2=zeros(l,l);

M2_2_1=zeros(l,l);
M2_2_2=zeros(l,l);
M4_2_1=zeros(l,l);
M4_2_2=zeros(l,l);

M3_3_1=zeros(l,l);
M3_3_2=zeros(l,l);

for i=1:l
    for j=1:i
     
        if number(i,1)==number(j,1)
            M2_0(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j))).*r.^3,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
            M4_0(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j))).*r.^5,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
            M6_0(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j))).*r.^7,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
        end
        
        if number(i,1)-number(j,1)==1
            M1_1_1(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^2,0, 1,'RelTol',1e-8,'AbsTol',1e-12);

            
            M3_1_1(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^4,0, 1,'RelTol',1e-8,'AbsTol',1e-12);

            
            M5_1_1(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^6,0, 1,'RelTol',1e-8,'AbsTol',1e-12);

        end
        
        if number(i,1)-number(j,1)==-1
            
            M1_1_2(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^2,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
            M3_1_2(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^4,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
            M5_1_2(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^6,0, 1,'RelTol',1e-8,'AbsTol',1e-12);

        end
        
        if number(i,1)-number(j,1)==2
            
            M2_2_1(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^3,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
            M4_2_1(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^5,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
        end
        
        if number(i,1)-number(j,1)==-2
            
            M2_2_2(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^3,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
            M4_2_2(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^5,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
        end
        
        if number(i,1)-number(j,1)==3
            
            M3_3_1(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^4,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
        end
        
        if number(i,1)-number(j,1)==-3
            
            M3_3_2(i,j)=2*pi*normalized(i)*normalized(j)*quadgk(@(r)...
            (besselj(number(i,1), r*zero_point(i)).*besselj(number(j,1), r*zero_point(j)))...
            .*r.^4,0, 1,'RelTol',1e-8,'AbsTol',1e-12);
            
        end
    end
    disp(i)
end

for i=1:l
    for j=i+1:l
        
        M2_0(i,j)=M2_0(j,i);
        M4_0(i,j)=M4_0(j,i);
        M6_0(i,j)=M6_0(j,i);
        
        M1_1_1(i,j)=M1_1_1(j,i);
        M1_1_2(i,j)=M1_1_2(j,i);
        M3_1_1(i,j)=M3_1_1(j,i);
        M3_1_2(i,j)=M3_1_2(j,i);
        M5_1_1(i,j)=M5_1_1(j,i);
        M5_1_2(i,j)=M5_1_2(j,i);
        
        M2_2_1(i,j)=M2_2_1(j,i);
        M2_2_2(i,j)=M2_2_2(j,i);
        M4_2_1(i,j)=M4_2_1(j,i);
        M4_2_2(i,j)=M4_2_2(j,i);
        
        M3_3_1(i,j)=M3_3_1(j,i);
        M3_3_2(i,j)=M3_3_2(j,i);
    end
end

save([pwd,'/data_M.mat'],'M2_0','M4_0','M6_0','M1_1_1','M1_1_2','M3_1_1',...
    'M3_1_2','M5_1_1','M5_1_2','M2_2_1','M2_2_2','M4_2_1','M4_2_2','M3_3_1','M3_3_2');