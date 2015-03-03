function [Norm,Vals,s,ExtPts] = NormAlg(p)
% NORMALG Calculate the max-norm on [-1,1] of a real algebraic polynomial
% NORMALG(p) returns the max-norm on [-1,1] of a polynomial whose coefficients are stored in decreasing order in the vector p
% [Norm,Vals,s,ExtPts]=NORMALG(p) also returns in ExtPts the size(p) points where p achieves local extrema (including endpoints), in Vals the positive and negative values of the polynomial at these points, and in s the size of the vector Vals

deg=size(p,1)-1;
CritPts=roots(polyder(p));
m=size(CritPts,1);
if m<deg-1
  CritPts=[zeros(deg-1-m,1);CritPts];
end  
ExtPts=[-1;real(CritPts);1];
for i=2:deg
  Pti=ExtPts(i);
  if Pti > 1
    ExtPts(i) = 1;
  elseif Pti < -1
    ExtPts(i) = -1;
  end
end
ExtPts=sort(ExtPts);
Vals=polyval(p,ExtPts);
Vals=[Vals;-Vals];
Norm=max(Vals);
s=2*deg+2;

end