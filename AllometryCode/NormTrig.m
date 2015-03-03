function [Norm,Vals,s] = NormTrig(q)
% NORMTRIG Calculate the max-norm of a real trigonometric polynomial
% NORMTRIG(q) returns the max-norm of the trigonometric polynomial whose coefficients on the basis (cos(nx),sin(nx),...,cos(x),sin(x),1) are stored in the vector q
% [Norm,Vals,s]=NORMTRIG(q) also returns in Vals the positive and negative values of the polynomial at the points where it might achieve local extrema, and in s the size of the vector Vals

deg=(size(q,1)-1)/2;
if ( (deg~=1) & (deg~=2) )
  disp('Error in NormTrig: the code to compute the norm of a trigonometric polynomial of degree >2 is not yet written')
  return
elseif deg==1
  a=q(1);
  b=q(2);
  c=q(3);
  Vals=[sqrt(a*a+b*b)+c;sqrt(a*a+b*b)-c];
  Norm=max(Vals);
  s=2;
elseif deg==2
  a=q(1);
  b=q(2);
  c=q(3);
  d=q(4);
  e=q(5);
  C(1)=16*(a^2+b^2);
  C(2)=8*(b*d+a*c);
  C(3)=c^2+d^2-16*(a^2+b^2);
  C(4)=-8*a*c-4*b*d;
  C(5)=4*b^2-c^2;
  u=roots(C);
  m=size(u,1);
  if m<4
    u=[zeros(4-m,1);u];
  end
  u=[-1;real(u);1];
  for i=2:5
    ui=u(i);
    if ui>1
      u(i)=1;
    elseif ui<-1
      u(i)=-1;
    end 
  end
  u=sort(u);
  Vals=[a*(2*u.^2-1)+2*b*u.*sqrt(1-u.^2)+c*u+d*sqrt(1-u.^2)+e;...
        a*(2*u.^2-1)-2*b*u.*sqrt(1-u.^2)+c*u-d*sqrt(1-u.^2)+e];
  Vals=[Vals;-Vals];
  Norm=max(Vals);
  s=24;
end

end