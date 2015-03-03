function [k,OptDual] = AbsCdN(V,type,OptionalArg)
% ABSCDN Calculate the absolute condition number
% ABSCDN(V,type,OptionalArg) returns the absolute \ell_\infty-condition number of the system formed by the columns of V in 
% type 0: the Euclidean space \ell_2^n
% type 1: the sequence space \ell_p^n, in which case OptionalArg holds the value of p 
% type 2: the space P_m of algebraic polynomials of degree at most m
% type 3: the space T_m of trigonometric polynomials of degree at most m
% type 4: the polygonal space whose unit ball possesses the columns of OptionalArg as vertices 
% type 5: the polygonal space whose dual unit ball possesses the columns of OptionalArg as vertices
% [k,OptDual]=ABSCDN(V,type,OptionalArg) also returns the dual system OptDual which yields the minimal relative condition number


if nargin<2
  disp('Error in AbsCdN: the function takes at least two arguments')
  return
elseif ( (nargin==2) & (type==1) )
  disp('Error in AbsCdN: for \ell_p^n spaces [type 1], the value of p has to be specified')
  return
elseif ( (nargin==2) & (type==4) )
  disp('Error in AbsCdN: for polygonal spaces [type 4], the set of vertices of the unit ball has to be specified')
  return
elseif ( (nargin==2) & (type==5) )
  disp('Error in AbsCdN: for polygonal spaces [type 5], the set of vertices of the dual unit ball has to be specified')
  return
end 

[n,N]=size(V);

% determination of s2=min \|T'\|
% minimize_U max_j \|u_j\|_1 subject to V*U'=I using fminimax 
% we set x=[c_1;...;c_n] where U'=[c_1,...,c_n], so that the duality constraint VU'=I becomes V c_i = e_i, i=1:n.  

% linear equality constraints -- replace VU'=I by Wx=y
W=V;
for i=2:n
  W=[V, zeros(n,(i-1)*N);zeros((i-1)*n,N),W];
end
y=zeros(n*n,1);
for i=1:n
  y((i-1)*n+i)=1;
end

% choice of the initial vector x0
U0=CanDual(V);
x0=Mx2Vc(U0');

% definition of the objective function
function f = objfun(x)
U=Vc2Mx(x,N,n)';
if type==0
  for j=1:N
    uj=U(:,j);
    f(j)=norm(uj);
  end
elseif type==1
  p=OptionalArg;
  if p==1
    q=inf;
  else
    q=p/(p-1);
  end
  for j=1:N
    uj=U(:,j);
    f(j)=norm(uj,q);
  end
elseif type==2
  for j=1:N
    uj=U(:,j);
    [nj,Valsj,s]=DualNormAlg(uj);
    for l=1:s
      f(s*(j-1)+l)=Valsj(l);
    end
  end
elseif type==3
  for j=1:N
    uj=U(:,j);
    [nj,Valsj,s]=DualNormTrig(uj);
    for l=1:s
      f(s*(j-1)+l)=Valsj(l);
    end
  end
elseif type==4
  Vertices=OptionalArg;
  for j=1:N
    uj=U(:,j);
    [nj,Valsj,s]=DualNormPolygon(Vertices,uj);
    for l=1:s
      f(s*(j-1)+l)=Valsj(l);
    end
  end
elseif type==5
  Vertices=OptionalArg;
  for j=1:N
    uj=U(:,j);
    [nj,Valsj,s]=DualNormPolygon2(Vertices,uj);
    for l=1:s
      f(s*(j-1)+l)=Valsj(l);
    end
  end
end  
end

% the minimization itself
% options=optimset('GradObj','on');
% options=optimset('MaxFunEvals',500000);
[xx,fval,s2]=fminimax(@objfun,x0,[],[],W,y,[],[],[]);


% definition of s1 = \|T\|
function f = TheNorm(x)
if type==0
  f = norm(x);
elseif type==1 
  p=OptionalArg;
  f = norm(x,p);
elseif type==2
  f = NormAlg(x);
elseif type==3
  f = NormTrig(x);
elseif type==4
  Vertices=OptionalArg;
  f = NormPolygon(Vertices,x);
elseif type==5
  Vertices=OptionalArg;
  f = NormPolygon2(Vertices,x);
end
end
VEps=V*Signs1(N);
for h=1:2^(N-1)
  NormVEps(h) = TheNorm(VEps(:,h));
end
s1 = max(NormVEps);


% outputs
k=s1*s2;
OptDual=Vc2Mx(xx,N,n)';


end