param n := 4;
param b{i in 1..n};
let b[1] := 8;
let b[2] := 18;
let b[3] := 44;
let b[4] := 114;
var x{i in 1..n}>=0,<=n ;
for{i in 1..n} 
{let x[i] := 2};


minimize obj: sum{i in 1..n}( (sum{j in 1..n} x[j]^i) -b[i]  )^2;
