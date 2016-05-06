function [NMI PWF Modu Ncut]=evalClustering(W,K,Z,realZ)
%% W adjacent matrix
%% K number of cluters
%% Z returned community label 1:K
%% realZ real label 1:K
%% Modu modularity
%% Ncut
%% NMI normalized mutual information
%% PWF pairwise F measure

Modu=[];
Ncut=[];
NMI=[];
PWF=[];


%% Compute F measure and NMI
if isempty(realZ)==0
    [R  P NMI]=infoeval(Z,realZ);    
    n=length(Z);
    pair_relation_esti=repmat(Z, 1,n)-repmat(Z', n, 1);
    pair_relation_real=repmat(realZ, 1,n)-repmat(realZ', n, 1);
    
    pair_relation_esti=pair_relation_esti+diag(ones(1,n));
    pair_relation_real=pair_relation_real+diag(ones(1,n));
    
    Number_Pairs_Actually_in_Same_Cluster=length(find(pair_relation_real==0));
    Number_Pairs_Predicted_in_Same_Cluster=length(find(pair_relation_esti==0));
    ind_esti=find(pair_relation_esti==0);
    ind_real=find(pair_relation_real==0);
    ind_corr=intersect(ind_esti, ind_real);
    if length(ind_esti)>0
        Pre=length(ind_corr)/length(ind_esti);
    else
        Pre=0;
    end
    if length(ind_real)>0
       Rec=length(ind_corr)/length(ind_real);
    else
        Rec=0;
    end
    if Pre+Rec>0
        PWF=(2*Pre*Rec)./(Pre+Rec);
    else
        PWF=0;
    end
end


% compute Ncut and Modu
if isempty(W)==0
    n=size(W,1);
    if size(Z,2)==1
        I=[1:n]';
        J=Z;
        Z=sparse(I,J, ones(n,1), n,K);
    end
    Modu=0;
    AllSimi=sum(sum(W));

    Ncut=0;
    delta_Comm=zeros(n,n);
    for k=1:K
        I=find(Z(:,k));
        delta_Comm(I,I)=1;
        J=find(Z(:,k)==0);
        CW=W(I,I);
        BW=W(I,:);
        DW=W(I,J);
%         InCSimi=sum(sum(CW));
%         BeCSimi=sum(sum(BW));
%         Modu=Modu+(InCSimi/AllSimi-(BeCSimi/AllSimi)^2);

        volC=sum(sum(BW));
        cutC=sum(sum(DW));
        Ncut=Ncut+cutC/volC;
    end
     delta_Comm=delta_Comm-diag(diag(delta_Comm));
     din=sum(W,1);
     dout=sum(W,2);
     m=sum(din)/2;
     Modu=(1/(2*m)).*sum(sum((W-(dout*din)./(2*m)).*delta_Comm));
     Modu=full(Modu);
end




