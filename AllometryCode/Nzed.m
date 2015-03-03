function [Vstar,Ustar] = Nzed(V,U)
% NZED Calculate the best normalization of a frame relatively to a dual frame
% [Vstar,Ustar]=NZED(V,U) returns matrices Vstar and Ustar with columns \|u_j\| v_j and u_j / \|u_j\|, where v_j and u_j are the columns of V and U
% [Vstar,Ustar]=NZED(V) performs the same task with U=CanDual(V) 


[n,N]=size(V);

if nargin==1
  U=CanDual(V);
elseif nargin==2
  ZERO=1.e-9;
  if (norm(V*U'-eye(n),inf)>ZERO)
    disp('Error in Nzed: the given frames are not dual');
    return
  end
end

for j=1:N
  d(j)=norm(U(:,j));
end
D=diag(d);
Vstar=V*D;
Ustar=U/D;

end