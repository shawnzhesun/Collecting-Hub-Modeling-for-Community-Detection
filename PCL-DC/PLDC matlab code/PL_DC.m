function net=PL_DC(W, X, K, net, swit)
% this function implements the PL model combined with Discriminative Content Model
% input: W directional adjacent matrix
%        K number of clusters
%        X attributes of nodes if available, otherwise empty []
%        net structure
%            net.verbosity display level
%            net.opt optimization option opt{1}.tol opt{1}.maxit for EM,
%            opt{2}.tol opt{2}.maxit for M step, opt{3}.tol opt{3}.maxit
%            for solving w_k at each iteration in M-step maximization
%            net.state for rand state
%            net.G initial membership of freee form (n*K)
%            net.Y inital membership of specic form (n*K)
%            net.B initial assignments of popularity (n*1)
%            net.A initial assignments of productivity (n*1)
%            net.R inital popularity for each node in each community (n*K)
%            net.H inital productivity for each node in each community (n*K)
%            net.lamda  the regularizer parameter 
%            net.model store the final feature weights for each community (d*K),
%            i.e. w in the paper 
%        swit: the switch to choose 1:PL-1(popularity link model) 2:PL-2(productivity link model)
% output: net with the above fields
% 5-27-2009 Tianbao Yang version final
% usage net=PL_DC(W, X, K, net, swit)
format long


if nargin<5
    help PL_DC
end
if swit~=1 && swit~=2
    help PL_DC
end
    

n=size(W,1);
W=W-diag(diag(W));
W=W./sum(sum(W));



if isempty(net.state)
    state=1;
    net.state=state;
else
    state=net.state;
end
rand('state',state);

if isempty(net.verbosity)
    verbosity=1;
    net.verbosity=verbosity;
else
    verbosity=net.verbosity;
end


opt=net.opt;
if  isempty(opt{1}.tol)
    tol=1e-5;
    maxit=1000;
    stop='fulliteration';
else
    tol=opt{1}.tol;
    maxit=[];
    stop='converge';
end

if  isempty(opt{1}.maxit)==0
    maxit=opt{1}.maxit;
    stop='fulliteration';
end

if isempty(net.G)
    G=rand(n,K);
    G=G./repmat(sum(G,2), 1,K);
else
    G=net.G;
end

if isempty(net.B)
    B=rand(n,1);
    B=B./sum(B);
else
    B=net.B;
end

if isempty(net.A)
    A=rand(n,1);
    A=A./sum(A);
else
    A=net.A;
end


if isempty(net.R)
    R=repmat(B,1,K).*G; % the community-dependent popularity
else
    R=net.R;
end

if isempty(net.H)
    H=repmat(A,1,K).*G; % the community-dependent productivity
else
    H=net.H;
end


if isempty(X)==0
    d=size(X,2);
    
    if isempty(net.model)
        model=zeros(d, K);
    else
        model=net.model;
    end
    
    if isempty(net.lamda)
        lamda=0;
        net.lamda=lamda;
    else
        lamda=net.lamda;
    end
    
else
    model=[];
    lamda=[];
end




Y=G;
%% construct data sturcture
[I J v]=find(W);
S=[I J v];


linkN=size(S,1);
din=sum(W,1);  % indegree
dout=sum(W,2); % outdegree

if swit==1
    A=dout;
    H=repmat(A,1,K).*G;
elseif swit==2
    B=din';
    R=repmat(B,1,K).*G;
end


% start iteration
obj_curt=objective1(S, Y, R, H, swit);
obj_last=obj_curt-1;
it=0;
while  1
    if obj_curt>=obj_last &&  (obj_curt-obj_last)<tol
        net.stop='converge';
        break;
    end
    if (strcmp(stop, 'fulliteration') && it>=maxit)
        net.stop='fulliteration';
        break;
    end
    %% E-step compute tau or eta
    if swit==1
        tau=sum(R,1);
    elseif swit==2
        eta=sum(H,1);
    end
    
    % compute q
    q=zeros(size(S,1), 2+K);
    if swit==1
        q(:,1:2)=S(:,1:2);
        q(:,3:K+2)=G(S(:,1), :).*R(S(:,2),:)./repmat(tau, size(S,1),1);
        q(:,3:K+2)=q(:,3:K+2)./repmat(sum(q(:,3:K+2),2), 1, K);
    elseif swit==2
        q(:, 1:2)=S(:,1:2);
        q(:, 3:K+2)=H(S(:,1), :).*G(S(:,2),:)./repmat(eta, size(S,1),1);
        q(:,3:K+2)=q(:,3:K+2)./repmat(sum(q(:, 3:K+2),2), 1, K);
    end
    
    % compute n_in, n_out, n_io
    S_Mat=[repmat(S(:,3), 1, K)];
    STQ=S_Mat.*q(:,3:K+2);
    %for i=1:n
    %    inlink_i=find(S(:,2)==i);
    %    n_in(i,:)=sum(STQ(inlink_i,:), 1);
    %    outlink_i=find(S(:,1)==i);
    %    n_out(i,:)=sum(STQ(outlink_i,:),1);
    %end
    n_in=zeros(n,K);
    n_out=zeros(n,K);
    for k=1:K
       n_in(1:max(S(:,2)),k)=accumarray(S(:,2), STQ(:,k));
       n_out(1:max(S(:,1)),k)=accumarray(S(:,1), STQ(:,k));
   end
    
    n_io=n_in+n_out;
    m=sum(n_out,1);
    if swit==1
        m_tau=m./tau;
    elseif swit==2
        m_eta=m./eta;
    end
    
    %% M-step
    ind_in0=find(din==0);
    ind_in1=find(din>0);
    ind_out0=find(dout==0);
    ind_out1=find(dout>0);
    
    %% node class 1
    ind_in0_out0=intersect(ind_in0, ind_out0);
    B(ind_in0_out0,:)=0;
    A(ind_in0_out0,:)=0;
    G(ind_in0_out0,:)=1/K;
    
    %% node class 2
    ind_in0_out1=intersect(ind_in0, ind_out1);
    if swit==1 % PL-1
        B(ind_in0_out1,:)=0;
        G(ind_in0_out1,:)=n_out(ind_in0_out1,:)./repmat(sum(n_out(ind_in0_out1,:),2),1,K);
    elseif swit==2
        A(ind_in0_out1,:)=sum(n_io(ind_in0_out1,:)./repmat(m_eta, length(ind_in0_out1),1),2);
        G(ind_in0_out1,:)=n_io(ind_in0_out1,:)./(repmat(m_eta, length(ind_in0_out1), 1).*repmat(A(ind_in0_out1), 1, K));
        G(ind_in0_out1,:)=G(ind_in0_out1,:)./repmat(sum(G(ind_in0_out1,:),2),1, K);
    end
    %% node class 3
    ind_in1_out0=intersect(ind_in1, ind_out0);
    if swit==1
        B(ind_in1_out0,:)=sum(n_io(ind_in1_out0,:)./repmat(m_tau, length(ind_in1_out0),1),2);
        G(ind_in1_out0,:)=n_io(ind_in1_out0,:)./(repmat(m_tau, length(ind_in1_out0), 1).*repmat(B(ind_in1_out0), 1, K));
        G(ind_in1_out0,:)=G(ind_in1_out0,:)./repmat(sum(G(ind_in1_out0,:),2),1,K);
    elseif swit==2
        A(ind_in1_out0,:)=0;
        G(ind_in1_out0,:)=n_in(ind_in1_out0,:)./repmat(sum(n_in(ind_in1_out0,:),2),1,K);
    end
    %% node class 4
    ind_in1_out1=intersect(ind_in1, ind_out1);
    if swit==1
        for l=1:length(ind_in1_out1)
            i=ind_in1_out1(l);
            
            % approach 1
            arg{1}=n_io(i,:)./m_tau;
            arg{2}=sum(n_out(i,:))./m_tau;
            arg{3}=1;
            [p0,err,k,y]=newton_solve_equation('fun','dfun',0,1e-10,1e-10,1000, arg{:});
            B(i)=p0;
            G(i,:)=n_io(i,:)./(B(i).*m_tau+sum(n_out(i,:)));
            G(i,:)=G(i,:)./sum(G(i,:));
        end
    elseif swit==2
        for l=1:length(ind_in1_out1)
            i=ind_in1_out1(l);
            % approach 1
            arg{1}=n_io(i,:)./m_eta;
            arg{2}=sum(n_in(i,:))./m_eta;
            arg{3}=1;
            [p0,err,k,y]=newton_solve_equation('fun','dfun',0,1e-10,1e-10,1000, arg{:});
            A(i)=p0;
            G(i,:)=n_io(i,:)./(A(i).*m_eta+sum(n_in(i,:)));
            G(i,:)=G(i,:)./sum(G(i,:));
        end
    end
    if swit==1
        B=B./sum(B);
    end
    if swit==2
        A=A./sum(A);
    end
    
    %% discriminative content model
    if isempty(X)==0
        if length(opt)==1
            if it==1
                disp('using BFGS to solve the optimization problem in M step');
            end
            st=LOGm_ml(struct('lambda',lamda,'intercept',false),G, X);
            model=st.para;
            Y=exp(X*model);
            Y=Y./repmat(sum(Y,2),1,K);
        elseif length(opt)==2
            if it==1
                disp('using Nesterov method to solve the optimization problem in M step');
            end
            [Y model]=LogisticModel_Nesterov(X, G, lamda, opt{2}, 1);
        elseif length(opt)==3
            if it==1
                disp('using Newton method to solve the optimization problem in M step');
            end
            [Y model]=LogisticModel_Newton(model, X, G, lamda, {opt{2:3}}, verbosity);
        end
    else
        Y=G;
    end
    
    %% for next iteration
    R=Y.*repmat(B,1,K);
    H=Y.*repmat(A,1,K);
    obj_last=obj_curt;
    obj_curt=objective1(S, Y, R, H, swit);
    it=it+1;
    obj(it)=obj_curt;
    if verbosity>=1
        disp(sprintf('1..........................................iteration=%d,  first Objective=%.20f', it, obj_curt));
    elseif verbosity==0 && mod(it, 100)==0
        disp(sprintf('1..........................................iteration=%d, first Objective=%.20f', it, obj_curt));
    end
end

net.G=G;
net.Y=Y;
net.B=B;
net.R=R;
net.A=A;
net.H=H;
net.obj=obj;
net.model=model;
net.lamda=lamda;

if swit==1
    net.algorithm='PL1';
elseif swit==2
    net.algorithm='PL2';
end





%%
function obj=objective1(S, G, R, H, swit)
n=size(G,1);
K=size(G,2);
tau=sum(R,1);
eta=sum(H,1);


obj=0;

if swit==1
    for m=1:size(S,1)
        i=S(m,1);
        j=S(m,2);
        sij=S(m,3);
        v1=zeros(1,K);
        ind=find(R(j,:)>0);
        v1(ind)=R(j,ind)./tau(ind);
        logterm=(G(i,ind)*v1(ind)');
        obj=obj+sij*log(logterm);
    end
end

if swit==2
    for m=1:size(S,1)
        i=S(m,1);
        j=S(m,2);
        sij=S(m,3);
        v1=zeros(1,K);
        ind=find(H(i,:)>0);
        v1(ind)=H(i,ind)./eta(ind);
        logterm=(v1(ind)*G(j,ind)');
        obj=obj+sij*log(logterm);
    end
end


function [Y model]=LogisticModel_Newton(model, X, G, lamda, opt, verbosity)
if  isempty(opt{1}.tol)
    tol=1e-5;
    maxit=1;
    stop='fulliteration';
else
    tol=opt{1}.tol;
    maxit=[];
    stop='converge';
end

if  isempty(opt{1}.maxit)==0
    maxit=opt{1}.maxit;
    stop='fulliteration';
end

[n K]=size(G);
[d K]=size(model);
it=0;
A=X*model;
obj_curt=sum(sum(G.*A)) - sum( log( sum(exp(A), 2) ) )- (lamda/2)*( sum(sum(model.*model) ) );
obj_last=-Inf;
% if verbosity>=2
%     disp(sprintf('2...............iteration=%d,  Logistic objective =%.10f', it, obj_curt));
% end
while  1
    if obj_curt>=obj_last && (obj_curt-obj_last)<tol
        stop='converge';
        break;
    end
    if (strcmp(stop, 'fulliteration') && it>=maxit)
        stop='fulliteration';
        break;
    end
    obj_last=obj_curt;
    tau= sum(exp(X*model), 2);
    ogparam{1}=G;
    ogparam{2}=repmat(1./tau, 1, K);
    ogparam{3}=X;
    ogparam{4}=lamda;
    [model A]=newton_raphson_model(model, A, ogparam, opt{2}, verbosity);
    obj_curt=sum(sum(G.*A)) - sum( log( sum(exp(A), 2) ) )- (lamda/2)*( sum(sum(model.*model) ) );
    it=it+1;
    if verbosity>=2 && mod(it, 10)==0
        disp(sprintf('2...............iteration=%d,  Logistic objective =%.10f', it, obj_curt));
    end
    if verbosity>=3
        disp(sprintf('2...............iteration=%d,  Logistic objective =%.10f', it, obj_curt));
    end
end

maxA=repmat(max(A,[],2), 1, K);
Y=exp(A);
Y=Y./repmat(sum(Y,2), 1, K);


function [m Ey err]=newton_raphson_model(m, Ey,  param, opt, verbosity)
if  isempty(opt.tol)
    tol=1e-5;
    maxit=200;
else
    tol=opt.tol;
    maxit=200;
end

if  isempty(opt.maxit)==0
    maxit=opt.maxit;
end

p=param{1};
q=param{2};
X=param{3};
lamda=param{4};
[d K]=size(m);
lamda=repmat(lamda, 1, K);

Xp=X'*p;

obj_new=sum(sum(q.*exp(Ey))) + (1/2).*sum(lamda.* sum(m.*m, 1) ) - sum(sum(p.*Ey));

W=q.*exp(Ey);

dmp=Xp-X'*W;

dm=repmat(lamda, d, 1).*m - dmp;

obj_old=-Inf;
it=0;
if verbosity>=3&& mod(it, 10)==0
    disp(sprintf('3......iteration=%d,  newton mean Objective=%.10f', it, obj_new));
end
while obj_new > obj_old || obj_old-obj_new>tol              % begin Newton's iterations
    obj_old = obj_new;
    m_old=m;
    for k=1:K
        Bk=X'*diag(W(:,k))*X;
        Lk = chol(Bk+lamda(k)*eye(d));                            % L'*L=B=eye(n)+sW*K*sW
        bk = Bk*m(:,k)+dmp(:,k);
        m_init(:,k) = Lk\(Lk'\bk);
    end
    deltam=m_init-m_old;
    
    m=m_init;
    Ey=X*m_init;
    obj_new=sum(sum(q.*exp(Ey))) + (1/2).*sum(lamda.* sum(m.*m, 1) ) - sum(sum(p.*Ey));
    
    t=1;
    alpha=0.2;
    beta=0.1;
    while obj_new > obj_old + alpha*t*(sum(sum(dm.*deltam)))
        t=t*beta;
        m=(1-t).*m_old+ t.* m_init;
        Ey=X*m;
        obj_new=sum(sum(q.*exp(Ey))) + (1/2).*sum(lamda.* sum(m.*m, 1) ) - sum(sum(p.*Ey));
    end
    obj_new=sum(sum(q.*exp(Ey))) + (1/2).*sum(lamda.* sum(m.*m, 1) ) - sum(sum(p.*Ey));
    
    W=q.*exp(Ey);
    
    dmp=Xp-X'*W;
    
    dm=repmat(lamda, d, 1).*m - dmp;
    
    it=it+1;
    if verbosity>=3
        disp(sprintf('3......iteration=%d,  newton mean Objective=%.10f', it, obj_new));
    end
    if isempty(maxit)==0 &&it>maxit
        break;
    end
end
err=obj_new-obj_old;
maxd=max(max(abs(dm)));
if verbosity>=3
    disp(sprintf('3......iteration=%d,  err=%.10f, grad=%.10f', it, err,  maxd));
end


%%
function [Y model]=LogisticModel_Nesterov(X, G, lamda, opt, verbosity)
% this function implements the logist regression model using Nesterov
% method to optimize
if  isempty(opt.tol)
    tol=1e-5;
    maxit=1000;
    stop='fulliteration';
else
    tol=opt.tol;
    maxit=[];
    stop='converge';
end

if  isempty(opt.maxit)==0
    maxit=opt.maxit;
    stop='fulliteration';
end

n=size(G,1);
K=size(G,2);
d=size(X,2);
if size(X,1)~=n
    error('error: number of examples in X');
end

it=1;
model_Yt=zeros(K,d);
model_Yt_1=zeros(K,d);
t_1=1;
t_2=0;
Lt_1=10;
obj_Yt=objective2(X, G, model_Yt, lamda);
obj_Yt_1=obj_Yt-1;
while 1
    if (obj_Yt>=obj_Yt_1 && (obj_Yt-obj_Yt_1)<=tol)
        break;
    elseif ( strcmp(stop,'fulliteration') && it>=maxit)
        break;
    end
    alphat=(t_2-1)/t_1;
    model_Xt=(1+alphat).*model_Yt-alphat.*model_Yt_1;
    [obj_Xt tempY]=objective2(X, G, model_Xt, lamda);
    grad_Xt=gradient2(X, G, tempY, model_Xt, lamda);
    L=Lt_1;
    model_Xt_L=model_Xt+(1/L).*grad_Xt;
    obj_Xt_L=objective2(X, G, model_Xt_L, lamda);
    norm_grad_Xt=(norm(grad_Xt,'fro'))^2;
    while obj_Xt_L<obj_Xt+(1/(2*L))*norm_grad_Xt
        L=2*L;
        model_Xt_L=model_Xt+(1/L).*grad_Xt;
        obj_Xt_L=objective2(X, G, model_Xt_L, lamda);
    end
    Lt_1=L;
    model_Yt_1=model_Yt;
    model_Yt=model_Xt_L;
    obj_Yt_1=obj_Yt;
    obj_Yt=obj_Xt_L;
    t_2=t_1;
    t_1=(1/2)*(1+sqrt(1+4*t_2*t_2));
    it=it+1;
    if verbosity>=2
        disp(sprintf('2...........iteration=%d,  first Objective=%.10f', it, obj_Yt));
    end
end

model=model_Yt;
A=X*model';
maxA=max(A,[],2);
A=A-repmat(maxA, 1, K);
Y=exp(A);
Y=Y./repmat(sum(Y,2),1,K);

%%
function [obj Y]=objective2(X, G, model, lamda)
% this function computes the objective in logistic regression
n=size(G,1);
K=size(G,2);
obj=0;
A=X*model';
maxA=max(A,[],2);
A=A-repmat(maxA, 1, K);
Y=exp(A);
Y=Y./repmat(sum(Y,2),1,K);

%ind=find(G>0);
obj=G(:)'*log(Y(:))-(lamda/2)* (model(:)'*model(:));


function grad=gradient2(X, G, Y, model, lamda)
% this function computes the gradient in logistic regression
K=size(model,1);
d=size(model,2);

grad=zeros(K,d);
grad=(G-Y)'*X-lamda.*model;


%%
function [g b]=iterateGb(Ain, Aout, alpha, g,b, ep, N,stop, verbosity)

A=Ain+Aout;
obj_curt=objective3(Ain, Aout, alpha,g,b);
obj_last=obj_curt-1;
it=0;
while 1
    if obj_curt>=obj_last && (obj_curt-obj_last)<ep
        break;
    end
    if (strcmp(stop,'fulliteration') && it>=N)
        break;
    end
    %g=A./(b.*alpha+sum(Aout));
    %g=g./sum(g);
    arg{1}=A; arg{2}=b.*alpha; arg{3}=1;
    lamda=newton_solve_equation('fun','dfun',0,1e-10,1e-10,10000,arg);
    g=A./(b.*alpha+lamda);
    g=g./sum(g);
    b=sum(Ain)/(alpha*g');
    obj_last=obj_curt;
    obj_curt=objective3(Ain, Aout, alpha, g, b);
    it=it+1;
    if verbosity>=2
        disp(sprintf('2................iteration=%d,  Logistic Objective=%.10f', it, obj_curt));
    end
end

%%
function obj=objective3(Ain, Aout, alpha, g, b)
obj=0;
A=Ain+Aout;

ind=find(A);
obj=obj+A(ind)*log(g(ind))';

obj=obj+sum(Ain)*b;

obj=obj-b*alpha*g';


