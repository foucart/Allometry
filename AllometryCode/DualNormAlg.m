function [Norm,Vals,s] = DualNormAlg(lambda)
% DUALNORMALG Calculate the norm (subordinate to the max-norm on [-1,1]) of a linear functional on a space of real algebraic polynomials
% DUALNORMALG(lambda) returns the norm of the linear functional whose values on the monomials are stored in decreasing order in the vector lambda
% [Norm,Vals,s]=DUALNORMALG(lambda) also returns in Vals the positive and negative values of the linear functional on the extreme polynomials where it might achieve its norm, and in s the size of Vals 
 

deg=size(lambda,1)-1;
if deg>2
  disp('Error in DualNormAlg: the code to compute the dual norm on P_n is not yet written for n>2')
  return
elseif deg==2
  Alpha=lambda(1);
  Beta=lambda(2);
  Gamma=lambda(3);
  F0=2*Alpha-Gamma;
  Fp1=(Alpha-2*Beta-Gamma)/2;
  Fn1=(Alpha+2*Beta-Gamma)/2;
  tp=(Alpha+Beta)/(Gamma+Beta);
  tn=(Alpha-Beta)/(-Gamma+Beta);
  if ( (tp>0) & (tp<1) )
    Ftp=(2*Alpha-4*tp*Beta+(tp^2-2*tp-1)*Gamma)/(1+tp)^2;
  else
    Ftp=0;
  end
  if ( (tn<0) & (tn>-1) )
    Ftn=(2*Alpha-4*tn*Beta+(tn^2+2*tn-1)*Gamma)/(1-tn)^2; 
  else
    Ftn=0;
  end
  Vals=[Gamma;F0;Fp1;Fn1;Ftp;Ftn];
  Vals=[Vals;-Vals];
  Norm=max(Vals);
  s=12;
elseif deg==1 
  Alpha=lambda(1);
  Beta=lambda(2);
  Vals=[Alpha;Beta];
  Vals=[Vals;-Vals];
  Norm=max(Vals);
  s=4;
end 

end