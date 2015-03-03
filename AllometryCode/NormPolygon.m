function [Norm,Vals,s] = NormPolygon(V,x)
% NORMPOLYGON Calculate the norm in a polygonal space
% NORMPOLYGON(V,x) returns the norm of the vector x in the two-dimensional space whose unit ball possesses the columns of the matrix V as vertices
% [Norm,Vals,s]=NORMPOLYGON(V,x) also returns in Vals the values of which Norm is the maximum, and in s the size of Vals

[n,m]=size(V);
if n ~= 2
  disp('Error in NormPolygon: the vertices must be vectors of size 2')
  return
elseif (m/2 - floor(m/2)) ~= 0
  disp('Error in NormPolygon: the number of vertices must be even') 
  return
else
  for i=1:m-1
    u(i,:) = 1/(V(1,i)*V(2,i+1)-V(1,i+1)*V(2,i))*[V(2,i+1)-V(2,i),V(1,i)-V(1,i+1)];
  end
    u(m,:) = 1/(V(1,m)*V(2,1)-V(1,1)*V(2,m))*[V(2,1)-V(2,m),V(1,m)-V(1,1)];
  Vals=u*x;  
  Norm=max(Vals);
  s=m;
end

end 