$ontext
l. Variant of Golden and Wong (1981) with compact subtour elimination constraints
$offtext

sets
  i     set of nodes    / i1*i6 /
  d(i)  depot node      / i1 /
  k     set of available postmen or vehicles / k1*k3 / ;

alias(i,j);

sets
  e(i,j) set of arc
  r(i,j) set of arcs to be serviced ;

parameters
  c(i,j) costs of traversing arc i-j
  / i1.i2 1, i1.i3 6, i2.i3 2, i2.i4 2, i2.i6 5, i3.i4 1, i3.i5 2, i4.i5 1, i4.i6 2 /

  demand_q(i,j) demand on arc i-j or j-i
  / i1.i2 4, i2.i6 3, i4.i5 5 / ;
  
*parameters
*  c(i,j) costs of traversing arc i-j
*  / i1.i2 1.3, i1.i3 6.1, i2.i3 2, i2.i4 2.1, i2.i6 5, i3.i4 1, i3.i5 2, i4.i5 1, i4.i6 2 /
*
*  demand_q(i,j) demand on arc i-j or j-i
*  / i1.i2 4, i2.i6 3, i4.i5 5 / ;
*
  c(i,j)$(ord(i)>ord(j)) = c(j,i) ;
  demand_q(i,j)$(ord(i)>ord(j)) = demand_q(j,i) ;

  e(i,j) = no ;
  e(i,j)$c(i,j) = yes ;
  r(i,j) = no ;
  r(i,j)$demand_q(i,j) = yes ;
*  e(i,j) = e(i,j) + e(j,i);
*  c(i,j) = c(i,j) + c(j,i);
*  r(i,j) = r(i,j) + r(j,i);

display e,r,c,demand_q;
scalar
  W capacity per postman / 8/
  W2 capacity per postman / 38/ ;

variables
  z objective variable ;

positive variables
  f(i,j,k) flow variable ;

binary variables
  x(i,j,k) i-j is traversed by postman k
  l(i,j,k) i-j is serviced by postman k ;

equations
  objective_1          objective function
  constraint_2(i,k)    route continuity
  constraint_3(i,j)    arcs with positive demand are serviced exactly once
  constraint_4(i,j,k)  arcs i-j are serviced by postman k only if he covers it
  constraint_5(k)      capacities of vehicles are not violated
  constraint_5b(k)
  constraint_6a(i,k)   eliminate illegal subtours - part 1
  constraint_6b(i,j,k) eliminate illegal subtours - part 2 ;

  objective_1..                   z =e= sum((i,j,k)$e(i,j), c(i,j)*x(i,j,k)) ;
  constraint_2(i,k)..             sum(j$e(j,i), x(j,i,k)) =e= sum(j$e(i,j), x(i,j,k)) ;
  constraint_3(i,j)$r(i,j)..      sum(k, l(i,j,k)+l(j,i,k)) =g= 1 ;
  constraint_4(i,j,k)$e(i,j)..    x(i,j,k) =g= l(i,j,k) ;
  constraint_5(k)..               sum((i,j)$r(i,j), demand_q(i,j)*l(i,j,k)) =l= W ;
  constraint_5b(k)..              sum((i,j)$e(i,j), c(i,j)*x(i,j,k)) =l= W2 ;
  constraint_6a(i,k)$(not d(i)).. sum(j$e(i,j), f(i,j,k)) - sum(j$e(j,i), f(j,i,k)) =e= sum(j$e(i,j), l(i,j,k)) ;
  constraint_6b(i,j,k)$e(i,j)..   f(i,j,k) =l= sqr(card(i))*x(i,j,k) ;


model carp / all /;

option MIP = cplex ;

solve carp minimizing z using MIP ;

display x.l, l.l ;
