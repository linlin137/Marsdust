%
% All rights are retained by the authors and Tsinghua University and University of Stuttgart.
% Please contact linhl24@mails.tsinghua.edu.cn for licensing inquiries.
% 
% Authors: Hanlin Lin
% Contact: linhl24@mails.tsinghua.edu.cn
% 

close all;

clear;

l=2990;
dt=10;
N=256;

rxy1=zeros(1,l/dt+1);
rxy2=zeros(1,l/dt+1);
rxy3=zeros(1,l/dt+1);
rxy4=zeros(1,l/dt+1);
rxy5=zeros(1,l/dt+1);
rxy6=zeros(1,l/dt+1);
rxy12=zeros(1,l/dt+1);
rxy13=zeros(1,l/dt+1);

asp1=zeros(1,l/dt+1);
asp2=zeros(1,l/dt+1);
asp3=zeros(1,l/dt+1);
asp4=zeros(1,l/dt+1);
asp5=zeros(1,l/dt+1);
asp6=zeros(1,l/dt+1);
asp12=zeros(1,l/dt+1);
asp13=zeros(1,l/dt+1);

t=0:dt:l;

path="rare_25_N256TIME600.0seed0.2976Re0.000055Nq0R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

    z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp1(k/dt+1)=2*max(rr)/height/2;
    rxy1(k/dt+1)=max(rr);
    clear rr zz;
end



path="rare_25_N256TIME600.0seed0.2976Re0.000055Nq1000R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

      z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp2(k/dt+1)=2*max(rr)/height/2;
    rxy2(k/dt+1)=max(rr);
    clear rr zz;
end

path="rare_25_N256TIME600.0seed0.2976Re0.000055Nq3000R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

    z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp3(k/dt+1)=2*max(rr)/height/2;
    rxy3(k/dt+1)=max(rr);
    clear rr zz;
end

path="rare_25_N256TIME600.0seed0.2976Re0.000055Nq5000R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

    z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp4(k/dt+1)=2*max(rr)/height/2;
    rxy4(k/dt+1)=max(rr);
    clear rr zz;
end
path="rare_25_N256TIME600.0seed0.2976Re0.000055Nq10000R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

    z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp5(k/dt+1)=2*max(rr)/height/2;
    rxy5(k/dt+1)=max(rr);
    clear rr zz;
end

path="rare_25_N256TIME600.0seed0.2976Re0.000055Nq100000R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

    z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp6(k/dt+1)=2*max(rr)/height/2;
    rxy6(k/dt+1)=max(rr);
    clear rr zz;
end

path="rare_25_N256TIME600.0seed0.3976Re0.000055Nq0R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

    z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp12(k/dt+1)=2*max(rr)/height/2;
    rxy12(k/dt+1)=max(rr);
    clear rr zz;
end

path="rare_25_N256TIME600.0seed0.4976Re0.000055Nq0R00.008.txt";
data2=load(path);

for k=0:dt:l
    x=data2(k*N+1:1:(k+1)*N,1);
    y=data2(k*N+1:1:(k+1)*N,2);
    z=data2(k*N+1:1:(k+1)*N,3);

    z2=sort(z,'descend');
    height=(z2(1)-z2(N/2));
    
    count2=0;
    for j=1:1:N
        if(z(j)> z2(1)-2.3*height)
            count2=count2+1;
            rr(count2)=sqrt(x(j)*x(j)+y(j)*y(j));
            zz(count2)=z(j);
        end
    end
    asp13(k/dt+1)=2*max(rr)/height/2;
    rxy13(k/dt+1)=max(rr);
    clear rr zz;
end

datafilename="2ijnew.mat";
save(datafilename,"N","dt","l","rxy1","rxy2","rxy3","rxy4","rxy5","rxy6","rxy12","rxy13","asp1","asp2","asp3","asp4","asp5","asp6","asp12","asp13");