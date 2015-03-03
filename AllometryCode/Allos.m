function [k,OptSyst,OptDual,Comment,K]=Allos(n,N,type,OptionalArg,NbMzions)
% ALLOS Execute Allo several times
% ALLOS(n,N,type,OptionalArg) returns the value of the Nth allometry constant of the n-dimensional space
% type 0: the Euclidean space \ell_2^n
% type 1: the sequence space \ell_p^n, in which case OptionalArg holds the value of p 
% type 2: the space P_m of algebraic polynomials of degree at most m
% type 3: the space T_m of trigonometric polynomials of degree at most m
% type 4: the polygonal space whose unit ball possesses the columns of OptionalArg as vertices 
% type 5: the polygonal space whose dual unit ball possesses the columns of OptionalArg as vertices
% after having performed the minimization 10 times
% [k,OptSyst,OptDual,Comment,K]=ALLOS(n,N,type,OptionalArg) also returns best conditioned N-system OptSyst and its optimal dual system OptDual, a remark Comment on the reliability of the result, and the table K of the 10 values obtained


if nargin<3
  disp('Error in Allos: the function takes at least three arguments')
  return
elseif ( (nargin==3) & (type==1) )
  disp('Error in Allos: for \ell_p^n spaces [type 1], the value of p has to be specified')
  return
elseif ( (nargin==3) & (type==4) )
  disp('Error in Allos: for polygonal spaces [type 4], the set of vertices of the unit ball has to be specified')
  return
elseif ( (nargin==3) & (type==5) )
  disp('Error in Allos: for polygonal spaces [type 5], the set of vertices of the dual unit ball has to be specified')
  return
elseif nargin==3
  OptionalArg=0;
  NbMzions=10;
elseif nargin==4
  NbMzions=10;
end

k=n;
NbOcc=1;
for l=1:NbMzions
  [kkk,OS,OD]=Allo(n,N,type,OptionalArg);
%  K(l)=kk;
  kk=CdN(OS,OD,type,OptionalArg);
  K(l)=kk;
  ZERO=1.e-4;
  if kk<k-ZERO
    k=kk;
    OptSyst=OS;
    OptDual=OD;
    NbOcc=1;
  elseif abs(k-kk)<=ZERO
    NbOcc=NbOcc+1;
    OptSyst=OS;
    OptDual=OD;
  end
end
if NbOcc==NbMzions
  Comment=['The result is very reliable: with ' int2str(NbMzions) ' minimizations, the value of the allometry constant was always the same'];
  for i=size(Comment,2)+1:140
    Comment(i)=' '; 
  end
elseif NbOcc>=2/3*NbMzions
  Comment=['The result is rather reliable: with ' int2str(NbMzions) ' minimizations, the smallest value for the allometry constant was obtained ' int2str(NbOcc) ' times'];
  for i=size(Comment,2)+1:140
    Comment(i)=' '; 
  end
else
  Comment=['The result is not so reliable: with ' int2str(NbMzions) ' minimizations, the smallest value for the allometry constant was obtained ' int2str(NbOcc) ' times'];
  for i=size(Comment,2)+1:140
    Comment(i)=' '; 
  end
end

end