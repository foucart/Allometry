function U = CanDual(V)
% CANDUAL Calculate the canonical dual frame
% CANDUAL(V) returns a matrix whose columns hold the canonical dual frame of the frame formed by the columns of V

S=V*V';
U=S\V;  

end