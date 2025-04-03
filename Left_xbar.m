function x = Left_xbar(x1,x2,y1,y2)
% centroid of the curvature area to the left end

if y1==0 & y2==0
    x   = 0.5*(x1+x2);
else
    Ax1 = (x2-x1)*y1*(x1+x2)/2;
    Ax2 = 0.5 * (x2-x1)*(y2-y1)*(x1+2/3*(x2-x1));
    A   = 0.5*(y1+y2)*(x2-x1);
    if A ~= 0
        x   = (Ax1+Ax2)/A;
    else
        x = 0.5*(x1+x2);
    end
    
end
