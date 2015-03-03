function k = CdN(V,U,type,OptionalArg)
% CDN Calculate the condition number
% CDN(V,U,type,OptionalArg) returns the \ell_\infty-condition number of the system formed by the columns of V relatively to the dual system formed by the columns of U, in 
% type 0: the Euclidean space \ell_2^n
% type 1: the sequence space \ell_p^n, in which case OptionalArg holds the value of p 
% type 2: the space P_n of algebraic polynomials of degree at most n
% type 3: the space T_n of trigonometric polynomials of degree at most n
% type 4: the polygonal space whose unit ball possesses the columns of OptionalArg as vertices 
% type 5: the polygonal space whose dual unit ball possesses the columns of OptionalArg as vertices

if nargin<3
  disp('Error in CdN: the function takes at least three arguments')
  return
elseif ( (nargin==3) & (type==1) )
  disp('Error in CdN: for \ell_p^n spaces [type 1], the value of p has to be specified')
  return
elseif ( (nargin==3) & (type==4) )
  disp('Error in CdN: for polygonal spaces [type 4], the set of vertices of the unit ball has to be specified')
  return
elseif ( (nargin==3) & (type==5) )
  disp('Error in CdN: for polygonal spaces [type 5], the set of vertices of the dual unit ball has to be specified')
  return
end

ZERO=1.e-4;
[n,N] = size(V);
if ( norm(V*U'-eye(n),inf)>ZERO )
  disp('Error in CdN: the given systems are not dual')
  k=inf;
  return
end


function f = TheNorm(x)
if type == 0
  f = norm(x);
elseif type == 1 
  p=OptionalArg;
  f = norm(x,p);
elseif type == 2
  f = NormAlg(x);
elseif type == 3
  f = NormTrig(x);
elseif type == 4
  Vertices=OptionalArg;
  f = NormPolygon(Vertices,x);
elseif type == 5
  Vertices=OptionalArg;
  f = NormPolygon2(Vertices,x);
end
end

function g = TheDualNorm(y)
if type == 0
  g = norm(y);
elseif type == 1
  p=OptionalArg;
  if p == 1
    q = inf;
  else
    q = p/(p-1);
  end
  g = norm(y,q);
elseif type == 2
  g = DualNormAlg(y);
elseif type == 3
  g = DualNormTrig(y);
elseif type == 4
  Vertices=OptionalArg;
  g = DualNormPolygon(Vertices,y);
elseif type == 5
  Vertices=OptionalArg;
  g = DualNormPolygon2(Vertices,y);
end
end


VEps=V*Signs1(N);
for h=1:2^(N-1)
  NormVEps(h)=TheNorm(VEps(:,h));
end
s1=max(NormVEps);

for j=1:N
  NormColU(j)=TheDualNorm(U(:,j));
end
s2=max(NormColU);

k=s1*s2;

end 