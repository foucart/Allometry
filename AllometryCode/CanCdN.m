function k = CanCdN(V)
% CANCDN Calculate the canonical condition number
% CANCDN(V) returns the condition number of the frame formed by the columns of V relatively to its canonical dual frame

U=CanDual(V);
k=CdN(V,U,0);

end