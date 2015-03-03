function [Norm,Vals,s] = NormPolygon2(V,x)
% NORMPOLYGON2 Calculate the norm in a polygonal space
% NORMPOLYGON2(V,x) returns the norm of the vector x in the two-dimensional space whose dual unit ball possesses the columns of the matrix V as vertices
% [Norm,Vals,s]=NORMPOLYGON2(V,x) also returns in Vals the values of the vertices of the dual unit ball applied to x, and in s the size of Vals

[n,m]=size(V);
if n ~= 2
  disp('Error in NormPolygon2: the vertices must be vectors of size 2')
  return
elseif (m/2 - floor(m/2)) ~= 0
  disp('Error in NormPolygon2: the number of vertices must be even') 
  return
else
  Vals=V'*x;
  Norm=max(Vals);
  s=m;
end

end