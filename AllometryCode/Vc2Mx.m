function V = Vc2Mx(x,n,N)
% VC2MX Produce a matrix whose columns are blocks of a vector
% VC2NX(x,n,N) returns the matrix [v_1,...,v_N] constructed from the vector x=[v_1;...;v_N]

if size(x,1)~=n*N
  disp('Error in Vc2Mx: the size of the input vector is not of the form n*N');
else 
  for j=1:N
    V(:,j)=x((j-1)*n+1:j*n);
  end
end

end