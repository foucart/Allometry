function Eps = Signs(N)
% SIGNS Produce vectors with entries +1 or -1
% SIGNS(N) returns in an N X 2^N matrix all the the 2^N vectors Eps of size N for which Eps(i)=+1 or -1, i>=1 

Eps=[Signs1(N),-Signs1(N)];

end






