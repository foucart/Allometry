function [Norm,Vals,s] = DualNormPolygon2(V,lambda)
% DUALNORMPOLYGON2 Calculate the dual norm in a polygonal space
% DUALNORMPOLYGON2(V,lambda) returns the norm, relatively to the two-dimensional space whose dual unit ball possesses the columns of the matrix V as vertices, of the linear functional whose values on [1;0] and [0;1] are given in lambda 
% [Norm,Vals,s]=DUALNORMPOLYGON2(V,lambda) also returns in Vals the values of the linear functional at the vertices of the unit ball, and in s the size of Vals

[n,m]=size(V);
if n ~= 2
  disp('Error in DualNormPolygon2: the vertices must be vectors of size 2')
  return
elseif (m/2 - floor(m/2)) ~= 0
  disp('Error in DualNormPolygon2: the number of vertices must be even') 
  return
else
  for i=1:m-1
    u(i,:) = 1/(V(1,i)*V(2,i+1)-V(1,i+1)*V(2,i))*[V(2,i+1)-V(2,i),V(1,i)-V(1,i+1)];
  end
    u(m,:) = 1/(V(1,m)*V(2,1)-V(1,1)*V(2,m))*[V(2,1)-V(2,m),V(1,m)-V(1,1)];
  Vals=u*lambda;  
  Norm=max(Vals);
  s=m;
end

end