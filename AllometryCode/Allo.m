function [k,OptSyst,OptDual] = Allo(n,N,type,OptionalArg)
% ALLO Calculates allometry constants
% ALLO(n,N,type,OptionalArg) returns the value of the Nth allometry constant of the n-dimensional space
% type 0: the Euclidean space \ell_2^n
% type 1: the sequence space \ell_p^n, in which case OptionalArg holds the value of p 
% type 2: the space P_m of algebraic polynomials of degree at most m
% type 3: the space T_m of trigonometric polynomials of degree at most m
% type 4: the polygonal space whose unit ball possesses the columns of OptionalArg as vertices 
% type 5: the polygonal space whose dual unit ball possesses the columns of OptionalArg as vertices
% [k,OptSyst,OptDual]=ALLO(n,N,type,OptionalArg) also returns the best conditioned N-system OptSyst and its optimal dual system OptDual



if nargin<3
  disp('Error in Allo: the function takes at least three arguments')
  return
elseif ( (nargin==3) & (type==1) )
  disp('Error in Allo: for \ell_p^n spaces [type 1], the value of p has to be specified')
  return
elseif ( (nargin==3) & (type==4) )
  disp('Error in Allo: for polygonal spaces [type 4], the set of vertices of the unit ball has to be specified')
  return
elseif ( (nargin==3) & (type==5) )
  disp('Error in Allo: for polygonal spaces [type 5], the set of vertices of the dual unit ball has to be specified')
  return
elseif nargin==3
  OptionalArg=0;
end 

% objective function
function f = objfun(x)
V=Vc2Mx(x(1:n*N),n,N);
AllEps=Signs1(N);
if type==0
  for h=1:2^(N-1)
    f(h) = norm(V*AllEps(:,h));
  end
elseif type==1
  p=OptionalArg;
  for h=1:2^(N-1)
    f(h) = norm(V*AllEps(:,h),p);
  end
elseif type==2
  for h=1:2^(N-1)
    [nh,Valsh,s]=NormAlg(V*AllEps(:,h));
    for l=1:s
      f(s*(h-1)+l)=Valsh(l);
    end
  end
elseif type==3
  for h=1:2^(N-1)
    [nh,Valsh,s]=NormTrig(V*AllEps(:,h));
    for l=1:s
      f(s*(h-1)+l)=Valsh(l);
    end
  end
elseif type==4
  Vertices=OptionalArg;
  for h=1:2^(N-1)
    [nh,Valsh,s]=NormPolygon(Vertices,V*AllEps(:,h));
    for l=1:s
      f(s*(h-1)+l)=Valsh(l);
    end
  end
elseif type==5
  Vertices=OptionalArg;
  for h=1:2^(N-1)
    [nh,Valsh,s]=NormPolygon2(Vertices,V*AllEps(:,h));
    for l=1:s
      f(s*(h-1)+l)=Valsh(l);
    end
  end
end 
end

% definition of the dual norm
function g = TheDualNorm(y)
if type==0
  g = norm(y);
elseif type==1
  p=OptionalArg;
  if p == 1
    q = inf;
  else
    q = p/(p-1);
  end
  g = norm(y,q);
elseif type==2
  g = DualNormAlg(y);
elseif type==3
  g = DualNormTrig(y); 
elseif type==4
  Vertices=OptionalArg;
  g = DualNormPolygon(Vertices,y);
elseif type==5
  Vertices=OptionalArg;
  g = DualNormPolygon2(Vertices,y);
end
end

% nonlinear equality constraints: duality VU'=I, normalization \|lambda_j\|=1
function [f,feq] = eqconstraints(x)
V=Vc2Mx(x(1:n*N),n,N);
U=Vc2Mx(x(n*N+1:2*n*N),n,N);
f=[];
for i=1:n
  for k=1:n
    feq((i-1)*n+k)=dot(U(i,:),V(k,:));
  end
  feq((i-1)*n+i)=feq((i-1)*n+i)-1;
end
for j=1:N
  Uj=U(:,j);
  feq(n*n+j)=TheDualNorm(Uj)-1;
end
end

% choice of the initial vector and minimization
NbRdTries=15;
V=randn(n,N);
[kk,U]=AbsCdN(V,type,OptionalArg);
for m=1:NbRdTries-1
  V1=rand(n,N);
  [kk1,U1]=AbsCdN(V1,type,OptionalArg);
  if kk1<kk
    kk=kk1;
    V=V1;
    U=U1;
  end
end
x0=[Mx2Vc(V);Mx2Vc(U)];
%options=optimset('GradObj','on');
options=optimset('MaxFunEvals',500000);
[x,fval,k]=fminimax(@objfun,x0,[],[],[],[],[],[],@eqconstraints,options);
 

% outputs
OptSyst=Vc2Mx(x(1:n*N),n,N);
OptDual=Vc2Mx(x(n*N+1:2*n*N),n,N);

end