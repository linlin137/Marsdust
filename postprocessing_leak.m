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

runout1=zeros(1,l/dt+1);
runout2=zeros(1,l/dt+1);
runout3=zeros(1,l/dt+1);
runout4=zeros(1,l/dt+1);
runout5=zeros(1,l/dt+1);
runout6=zeros(1,l/dt+1);
runout12=zeros(1,l/dt+1);
runout13=zeros(1,l/dt+1);
runout22=zeros(1,l/dt+1);
runout23=zeros(1,l/dt+1);
runout14=zeros(1,l/dt+1);
runout24=zeros(1,l/dt+1);

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
        end
    end
    runout1(k/dt+1)=N-count2;
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
        end
    end
    runout2(k/dt+1)=N-count2;
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
        end
    end
    runout3(k/dt+1)=N-count2;
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
        end
    end
    runout4(k/dt+1)=N-count2;
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
        end
    end
    runout5(k/dt+1)=N-count2;
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
        end
    end
    runout6(k/dt+1)=N-count2;
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
        end
    end
    runout12(k/dt+1)=N-count2;
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
        end
    end
    runout13(k/dt+1)=N-count2;
end




datafilename="2g1.mat";
save(datafilename,"N","dt","l","runout1","runout2","runout3","runout4","runout5","runout6","runout12","runout13");