function xx = smooth5(x)

xx(1) = x(1);
xx(2) = x(2);

for j= 3:max(size(x))-2
xx(j) = 1/5*(x(j-2)+x(j-1)+x(j)+x(j+1)+x(j+2));
end

xx(max(size(x))-1) = x(max(size(x))-1);
xx(max(size(x))) = x(max(size(x)));


xx = xx';

