function x = Mx2Vc(V)
% MX2VC Produce a vector in which the columns of a matrix are stacked on top of each other
% MX2VC(V) returns the vector [v1;...;v_n] constructed from the matrix V=[v_1,...,v_N]

[n,N]=size(V);
for j=1:N
  xx((j-1)*n+1:j*n)=V(:,j);
end
x=xx';

end