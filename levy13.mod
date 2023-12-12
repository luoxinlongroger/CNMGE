param pi := 3.1415926;
var x1 := 2;
var x2 := 2;



minimize obj:  sin(3*pi*x1)^2+ ((x1-1)^2)*(1+sin(3*pi*x2)^2) + ((x2-1)^2)*(1+sin(3*pi*x2)^2);
# 
#