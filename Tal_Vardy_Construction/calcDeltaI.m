function c = calcDeltaI( a, b, ap, bp )
% calculate the capacity change from merging y and yp
%   
y1 = calcCapacity(a,b);
y2 = calcCapacity(ap,bp);
as = a+ap;
bs = b+bp;
y3 = calcCapacity(as,bs);

c = y1 + y2 - y3;
end

