function V = Reg2nGone(n)
% REG2NGONE Construct a regular 2n-gone
% REG2NGONE(n) stores the coordinates of the vertices of the regular 2n-gone in a 2x2n matrix

V = [cos((1:2*n)*pi/n);sin((1:2*n)*pi/n)];

end