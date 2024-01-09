%% Addjust mass matrix M to consider support

function [Mwsupp]=ModMforSupp(M,supp)

    Mwsupp=M;
    A=size(supp);
    SuppSize=A(1);
    A=size(M);
    N=A(1);
    
    for i=1:SuppSize
        TempNodeNum=supp(i,1);
        if supp(i,2)==1

            Mwsupp(TempNodeNum*3-2,1:N)=zeros(1,N);
            Mwsupp(1:N,TempNodeNum*3-2)=zeros(N,1);
            Mwsupp(TempNodeNum*3-2,TempNodeNum*3-2)=0;

        end
        if supp(i,3)==1

            Mwsupp(TempNodeNum*3-1,1:N)=zeros(1,N);
            Mwsupp(1:N,TempNodeNum*3-1)=zeros(N,1);
            Mwsupp(TempNodeNum*3-1,TempNodeNum*3-1)=0;

        end
        if supp(i,4)==1

            Mwsupp(TempNodeNum*3,1:N)=zeros(1,N);
            Mwsupp(1:N,TempNodeNum*3)=zeros(N,1);
            Mwsupp(TempNodeNum*3,TempNodeNum*3)=0;

        end    
    end   
    
end