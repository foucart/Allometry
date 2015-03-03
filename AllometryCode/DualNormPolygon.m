function [Norm,Vals,s] = DualNormPolygon(V,lambda)
% DUALNORMPOLYGON Calculate the dual norm in a polygonal space
% DUALNORMPOLYGON(V,lambda) returns the norm, relatively to the two-dimensional space whose unit ball possesses the columns of the matrix V as vertices, of the linear functional whose values on [1;0] and [0;1] are given in lambda 
% [Norm,Vals,s]=DUALNORMPOLYGON(V,lambda) also returns in Vals the values of the linear functional at the vertices of the unit ball, and in s the size of Vals

[n,m]=size(V);
if n ~= 2
  disp('Error in DualNormPolygon: the vertices must be vectors of size 2')
  return
elseif (m/2 - floor(m/2)) ~= 0
  disp('Error in DualNormPolygon: the number of vertices must be even') 
  return
else
  Vals=V'*lambda;
  Norm=max(Vals);
  s=m;
end

end