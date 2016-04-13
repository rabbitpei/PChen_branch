clear;
clc;
close all;

e=exp(1);
for i=1:7
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.5;

p(1)=-0.07;
p(2)=-0.06;
p(3)=-0.05;
p(4)=-0.04;
p(5)=-0.03;
p(6)=0.006;
for i=7:27
    p(i)=(i-5)/50;
end
for i=1:27
    pr(i)=p(28-i);
end

pp(1)=0.2;
pp(2)=0.1;
pp(3)=0.005;
pp(4)=-0.1;

sample_num=10;
D= [-1 1 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 ;...
    1 0 1 0 0 0 0 ;...
    1 0 -1 -1 0 0 1 ;...
    0 0 0 0 -1 0 1 ;...
    0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 -1 -1];
TT=100;

my_edges=zeros(10,4,TT);
CC=zeros(7,4,sample_num);
nodein=zeros(2,3);
edgein=zeros(2,3);
for bt=1:TT
    for l=1:4
        q(l)=0.96^(1/abs(pp(l)));
        
        E=[-2/5*q(l) 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 ;...
            0 0 0 0 -7/5 0 0 ;...
            0 0 0 0 0 -8/5 0 ;...
            0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        pos=0;
        node_in1=0;
        node_in2=0;
        edge_in1=0;
        edge_in2=0;
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:7
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,2000);
        end
        TC=reshape(CC(:,l,:),7,sample_num);
        
        for e=1:sample_num
            vec1=[TC(1,1:e-1),TC(1,e+1:sample_num)];
            vec2=[TC(2,1:e-1),TC(2,e+1:sample_num)];
            vec3=[TC(5,1:e-1),TC(5,e+1:sample_num)];
            vec4=[TC(6,1:e-1),TC(6,e+1:sample_num)];
            aedge(1,e)=abs(corr(vec1',vec2'));
            aedge(2,e)=abs(corr(vec3',vec4'));
        end
        
        if l>=2
            for k=1:sample_num
                if TC(1,k)>mu1-sigma1&&TC(1,k)<mu1+sigma1
                    node_in1=node_in1+1;
                end
                if TC(5,k)>mu2-sigma2&&TC(5,k)<mu2+sigma2
                    node_in2=node_in2+1;
                end
                
                if aedge(1,k)>mu3-sigma3&&aedge(1,k)<mu3+sigma3
                    edge_in1=edge_in1+1;
                end
                if aedge(2,k)>mu4-sigma4&&aedge(2,k)<mu4+sigma4
                    edge_in2=edge_in2+1;
                end
            end
            nodein(1,l-1)=nodein(1,l-1)+node_in1/TT/sample_num;
            nodein(2,l-1)=nodein(2,l-1)+node_in2/TT/sample_num;
            edgein(1,l-1)=edgein(1,l-1)+edge_in1/TT/sample_num;
            edgein(2,l-1)=edgein(2,l-1)+edge_in2/TT/sample_num;
        end
        
        [mu1,sigma1]=normfit(TC(1,:));
        [mu2,sigma2]=normfit(TC(5,:));
        [mu3,sigma3]=normfit(aedge(1,:));
        [mu4,sigma4]=normfit(aedge(2,:));
        l
    end
    
    bt
end

nodein
edgein








