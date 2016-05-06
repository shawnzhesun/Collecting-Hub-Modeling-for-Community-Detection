function net=test_PL_DC()

% only link 
load emptynet;
load toy
net=PL_DC(W,[],2,net,1);

%link plus content
% load emptynet;
% opt=net.opt{1};
% net.opt=[];
% net.opt{1}=opt;
% net.lamda=0.1;
% %% using BFGS: I am using package realsed by others, I can not offer it 
% net=PL_DC(W,X,2,net,1);

%% using nesterov method solve the optimization in M-step
load emptynet;
opt=net.opt{1};
net.opt=[];
net.opt{1}=opt;
net.opt{2}=opt;
net.lamda=0.1;
net=PL_DC(W,X,2,net,1);


%% using newton method solve the optimization in M-step
load emptynet;
net.lamda=0.1;
%% using Nesterov
net=PL_DC(W,X,2,net,1);
