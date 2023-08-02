function [ang] = GetPhase2(a,b)

ang = atan(a/b);

c = a;
if c>0 && b >0
    ang = ang;
elseif c>0 && b<0
    ang = ang + pi;
elseif c<0 && b<0
    ang = ang + pi;
elseif c<0 && b>0
    ang = ang + 2*pi;
end






















