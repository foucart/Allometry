function [Norm,Vals,s] = DualNormTrig(lambda)
% DUALNORMTRIG Calculate the norm (subordinate to the max-norm) of a linear functional on a space of real trigonometric polynomials
% DUALNORMTRIG(lambda) returns the norm of the linear functional whose values on the basis (cos(nx),sin(nx),...,cos(x),sin(x),1) are stored in the vector lambda
% [Norm,Vals,s]=DUALNORMTRIG(lambda) also returns in Vals the positive and negative values of the linear functional on the extreme trigonometric polynomials where it might achieve its norm, and in s the size of Vals 


deg=(size(lambda,1)-1)/2;
if ( (deg~=1) & (deg~=2) )
  disp('Error in DualNormTrig: the code to compute the dual norm on T_n is not yet written for n>2')
  return
elseif deg==1
  Alpha=lambda(1);
  Beta=lambda(2);
  Gamma=lambda(3);
  Vals=[Gamma;-Gamma;sqrt(Alpha*Alpha+Beta*Beta)];
  Norm=max(Vals);
  s=3;
elseif deg==2
  Alpha=lambda(1);
  Beta=lambda(2);
  Gamma=lambda(3);
  Delta=lambda(4);
  Eps=lambda(5);
  C(1)=-4*( (Alpha*Eps+Delta^2-Gamma^2)^2+(Beta*Eps-2*Alpha*Delta)^2 );
  C(2)=-4*( (Alpha*Eps+Delta^2-Gamma^2)*(Beta*Delta+Alpha*Gamma-Gamma*Eps) + (Beta*Eps-2*Alpha*Delta)*(Beta*Gamma-Alpha*Delta-Delta*Eps) );
  C(3)=4*( (Alpha*Eps+Delta^2-Gamma^2)^2-(Beta*Eps-2*Alpha*Delta)*(2*Gamma*Delta-Beta*Eps) ) - ((Beta*Delta+Alpha*Gamma-Gamma*Eps)^2+(Beta*Gamma-Alpha*Delta-Delta*Eps)^2);
  C(4)=4*( (Alpha*Eps+Delta^2-Gamma^2)*(Beta*Delta+Alpha*Gamma-Gamma*Eps) ) - 2*(Beta*Gamma-Alpha*Delta-Delta*Eps)*(2*Gamma*Delta-Beta*Eps);
  C(5)=(Beta*Delta+Alpha*Gamma-Gamma*Eps)^2-(2*Gamma*Delta-Beta*Eps)^2;
  u=roots(C);
  m=size(u,1);
  if m<4
    u=[zeros(4-m,1);u];
  end
  u=sort([-1;real(u);1]);
  for i=1:6
    ui=u(i);
    if ui>i
      u(i)=1;
      D1(i)=0;
      D2(i)=0;
    elseif ui<-1
      u(i)=-1;
      D1(i)=0;
      D2(i)=0;
    else
      Den1=Gamma*sqrt(1-ui*ui)-Delta*ui;
      ScdDen1=Gamma*ui+Delta*sqrt(1-ui*ui)+Eps;
      if Den1~=0
        D1i=.5*(2*Alpha*ui*sqrt(1-ui*ui) - Beta*(2*ui*ui-1))/Den1;
      elseif ScdDen1~=0
        D1i=.5*(Alpha*(2*ui*ui-1)+Beta*2*ui*sqrt(1-ui*ui)+2*Gamma*ui+2*Delta*sqrt(1-ui*ui)+Eps)/ScdDen1;
      else
        D1i=0;
      end
      D1(i)=min(max(0,D1i),1);
      Den2=Gamma*sqrt(1-ui*ui)+Delta*ui;
      ScdDen2=Gamma*ui-Delta*sqrt(1-ui*ui)+Eps;
      if Den2~=0
        D2i=.5*(2*Alpha*ui*sqrt(1-ui*ui) + Beta*(2*ui*ui-1))/Den2;
      elseif ScdDen2~=0
        D2i=.5*(Alpha*(2*ui*ui-1)-Beta*2*ui*sqrt(1-ui*ui)+2*Gamma*ui-2*Delta*sqrt(1-ui*ui)+Eps)/ScdDen2;
      else
        D2i=0;
      end
      D2(i)=min(max(0,D2i),1);
    end 
  end
  Fat0=sqrt(Alpha*Alpha);
  [x,Fat1]=NormTrig([Alpha/4;Beta/4;-Gamma;-Delta;-Eps/4]);
  D1=D1';
  D2=D2';
  FatD1=(Alpha*(2*u.^2-1)+Beta*2*u.*sqrt(1-u.^2)-4*D1.*(Gamma*u+Delta*sqrt(1-u.^2)) + Eps*D1.*(D1-2)) ./ (1+D1).^2;
  FatD2=(Alpha*(2*u.^2-1)-Beta*2*u.*sqrt(1-u.^2)-4*D2.*(Gamma*u-Delta*sqrt(1-u.^2)) + Eps*D2.*(D2-2)) ./ (1+D2).^2;
  Vals=[Eps;FatD1;FatD2];
  Vals=[Vals;-Vals];
  Vals=[Fat0;Fat1;Vals];
  Norm=max(Vals);
  s=51;
end

end