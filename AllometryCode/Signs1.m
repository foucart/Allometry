function Eps = Signs1(N)
% SIGNS1 Produce vectors with entries +1 or -1
% SIGNS1(N) returns in an N X 2^(N-1) matrix all the the 2^(N-1) vectors Eps of size N for which Eps(1)=+1 and Eps(i)=+1 or -1, i>=2

Eps=[];
for i=1:N
  Eps=[ones(1,2^(i-1)); Eps, -Eps];
end

end