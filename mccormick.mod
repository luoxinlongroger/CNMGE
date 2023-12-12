param n := 2;
var x1:=2,>=-1.5,<=4;
var x2:=2,>=-3,<=4;


minimize obj: sin(x1+x2)+(x1-x2)^2-1.5*x1+2.5*x2+1;
