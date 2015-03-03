function [k,OptSyst,OptDual] = SymAllo(n,N,type)
% ALLO Calculates `symmetric' allometry constants
% ALLO(n,N,type,OptionalArg) returns the value of the Nth `symmetric' allometry constant of an n-dimensional space, only in the case
% type 2: the space P_m of algebraic polynomials of degree at most m
% [k,OptSyst,OptDual]=ALLO(n,N,type,OptionalArg) also returns the best conditioned symmetric N-system OptSyst and its optimal dual system OptDual



if nargin~=3
  disp('Error in SymAllo: the function takes exactly three arguments')
  return
elseif type~=2
  disp('Error in SymAllo: the third argument, type, must be equal to 2')
  return
end 

% objective function
function f = objfun(x)
V=Vc2Mx(x(1:n*N),n,N);
AllEps=Signs1(N);
  for h=1:2^(N-1)
    [nh,Valsh,s]=NormAlg(V*AllEps(:,h));
    for l=1:s
      f(s*(h-1)+l)=Valsh(l);
    end
  end
end

% nonlinear equality constraints: duality VU'=I, normalization \|lambda_j\|=1, and symmetry relation (-1)^{n-i} V_{i,N+1-j} = V_{i,j}
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
NN=floor((N+1)/2);
for j=1:NN
  Uj=U(:,j);
  feq(n*n+j)=DualNormAlg(Uj)-1;
end
for j=1:NN
  for i=1:n
    feq(n^2+NN+(j-1)*n+i)=(-1)^(n-i)*V(i,N+1-j)-V(i,j);
  end
end
end

% choice of the initial vector and minimization
NbRdTries=15;
V=randn(n,N);
[kk,U]=AbsCdN(V,type);
for m=1:NbRdTries-1
  V1=rand(n,N);
  [kk1,U1]=AbsCdN(V1,type);
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