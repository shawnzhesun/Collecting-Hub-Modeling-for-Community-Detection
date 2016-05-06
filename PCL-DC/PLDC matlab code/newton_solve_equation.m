function [p0,err,k,y]=newton_solve_equation(f,df,p0,delta,epsilon,maxN, varargin)

%Input - f is the object function input as a string 'f'
%      - df is the derivative of f input as a string 'df'
%      - p0 is the initial approximation to a zero of f
%        - delta is the tolerance for p0
%        - epsilon is the tolerance for the function values y
%        - max1 is the maximum number of iterations
%Output - p0 is the Newton-Raphson approximation to the zero
%         - err is the error estimate for p0
%         - k is the number of iterations
%         - y is the function value f(p0)
f=str2func(f);
df=str2func(df);
y=feval(f,p0,varargin{:});
for k=1:maxN    
    p1=p0-y./feval(df,p0,varargin{:});    
    err=abs(p1-p0);
    relerr=2*err/(abs(p1)+delta);
    p0=p1;
    y=feval(f,p0,varargin{:});
    %if (err<delta)|(relerr<delta)|(abs(y)<epsilon),break,end
    if (abs(y)<epsilon),break,end
end


function obj=fun(x, A,C,d)
% this funciton computes f(x)=\sum_k A_k/(C_k+x)-d
obj=sum(A./(C+x))-d;
 
function gra=dfun(x, A,C,d)
% this function computes f'(x) defnied above
gra=-sum(A./(C+x).^2);


function obj=fun2(x, arg)
% this function computes \sum_k a_kx/(b_kx+c)-d
a=arg{1};
b=arg{2};
c=arg{3};
d=arg{4};
obj=sum(a.*x/(b.*x+c))-d;


function gra=dfun2(x, arg)
% this function computes  the gradient of \sum_k a_kx/(b_kx+c)-d
a=arg{1};
b=arg{2};
c=arg{3};
d=arg{4};
gra=sum(a.*c/((b.*x+c).^2))


function obj=fun_ab(x, A, B, C, d)
% this function computes the objective in terms of a and b in APPL 
obj=sum(A./(B.*x(1)+C.*x(2)+d))-1;

function grad=dfun_ab(x,A, B, C, d)
% this funciton computes the gradient in terms of a and b in APPL
grad1=-sum((A.*B)./(B.*x(1)+C.*x(2)+d).^2);
grad2=-sum((A.*C)./(B.*x(1)+C.*x(2)+d).^2);
grad=[grad1; grad2];
