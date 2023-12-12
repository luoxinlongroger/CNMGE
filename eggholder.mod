param pi := 3.1415926;
var x1 := 2, >=-512, <=512;
var x2 := 2, >=-512, <=512;



minimize obj: -(x2+47) * sin(sqrt(abs(x2+x1/2+47)))-x1 * sin(sqrt(abs(x1-(x2+47))));
#
