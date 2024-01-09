%% Addjust stiffness matrix K to consider support

function [Kwsupp,T]=ModKforSupp(K,supp,Tinput)

    Kwsupp=K;
    A=size(supp);
    SuppSize=A(1);
    A=size(K);
    N=A(1);
    T=Tinput;
    
    for i=1:SuppSize
        TempNodeNum=supp(i,1);
        if supp(i,2)==1
            Kvv=K(TempNodeNum*3-2,TempNodeNum*3-2);
            Kwsupp(TempNodeNum*3-2,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3-2)=zeros(N,1);
            Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=Kvv;
            % Kwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=1;
            T(TempNodeNum*3-2)=0;
        end
        if supp(i,3)==1
            Kvv=K(TempNodeNum*3-1,TempNodeNum*3-1);
            Kwsupp(TempNodeNum*3-1,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3-1)=zeros(N,1);
            Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=Kvv;
            % Kwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=1;
            T(TempNodeNum*3-1)=0;
        end
        if supp(i,4)==1
            Kvv=K(TempNodeNum*3,TempNodeNum*3);
            Kwsupp(TempNodeNum*3,1:N)=zeros(1,N);
            Kwsupp(1:N,TempNodeNum*3)=zeros(N,1);
            Kwsupp(TempNodeNum*3,TempNodeNum*3)=Kvv;
            % Kwsupp(TempNodeNum*3,TempNodeNum*3)=1;
            T(TempNodeNum*3)=0;
        end    
    end   
    
end