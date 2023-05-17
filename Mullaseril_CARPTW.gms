$ontext

Model of Mullaseril (1997) with linearized time-window constraints
As of: 2023-03-19

$offtext

option resLim = 200 ;

* $call g********e instance_excel.xlsx o=instance_excel.gdx index=Parameter!A1:D7
        
* ****** instance_excel.gdx

scalars
   n_nod number of nodes
   n_veh number of cars
   c_veh capacity of cars
   n_arc number of required arcs
;

$load n_nod, n_veh, c_veh, n_arc

$eval nn n_nod+1
$eval na n_arc+1
$eval nv n_veh

sets
   i                  set of nodes  / i1*i%nn% /
   arc                index of arc  / 1*%na% /
   depot_arc(arc)     dummy depot arc
   required_arc(arc)  arcs that need to be served and dummy depot arc
   v                  vehicles      / k1*k%nv% /
   twb                time window boundaries / a, b /
   c_captions         column captions of data table / ord_i, ord_j, cost, demand /
;

alias(g, h, i, j) ;

set
   is_arc(arc,i,j) yes if i is initial and j is terminal node of arc arc

parameters
   data(arc,c_captions) imported data table
   ord_i(arc)           number of initial node i of arc arc
   ord_j(arc)           number of terminal node j of arc arc
   cost_a(arc)          cost of traversing arc arc
   demand_a(arc)        demand of arc
;

$load data

   ord_i(arc)$(ord(arc)<card(arc))    = data(arc,'ord_i') ;
   ord_j(arc)$(ord(arc)<card(arc))    = data(arc,'ord_j') ;
   cost_a(arc)$(ord(arc)<card(arc))   = data(arc,'cost') ;
   demand_a(arc)$(ord(arc)<card(arc)) = data(arc,'demand') ;
   ord_i(arc)$(ord(arc)=card(arc)) = %nn% ;
   ord_j(arc)$(ord(arc)=card(arc)) = 1 ;
   cost_a(arc)$(ord(arc)=card(arc)) = 0 ;
   demand_a(arc)$(ord(arc)=card(arc)) = 0 ;
   depot_arc(arc) = no ;
   depot_arc(arc)$(ord(arc)=card(arc)) = yes ;

   required_arc(arc) = no ;
   required_arc(arc)$(demand_a(arc) or depot_arc(arc)) = yes ;

   is_arc(arc,i,j) = no ;
   is_arc(arc,i,j)$(ord(i)=ord_i(arc) and ord(j)=ord_j(arc)) = yes ;

parameters
   c(i,j)          cost of arc i-j
;

   c(i,j) = sum(arc$(ord(i)=ord_i(arc) and ord(j)=ord_j(arc)), cost_a(arc)) ;
   c('i%nn%',i)$(ord(i)>1) = c('i1',i) ;
   c(i,'i%nn%') = c(i,'i1') ;
   c('i%nn%','i%nn%') = 0 ;
   
display c, ord_i, i,j, arc, demand_a;

** Assume symmetric network
   c(i,j)$(c(i,j)=0) = c(j,i) ;

$eval ne card(i)*(card(i)-1)/2

scalar bigM ;

   bigM = sum((i,j), c(i,j)) ;

* Die Zeitfenster muessten aus der Excel-Tabelle mit eingelesen werden
table   t_window(i,j,twb) time window for service of arc i-j or j-i ;


   t_window(i,j,'a')$(sum(arc, is_arc(arc,i,j))) = 0 ;
   t_window(i,j,'b')$(sum(arc, is_arc(arc,i,j))) = BigM ;


*******************************************************************************************************************

* Floyd-Warshall algorithm

parameters
   distance(i,j) shortest path length from node i to node j ;

   distance(i,j) = inf ;
   distance(i,i) = 0 ;
   distance(i,j)$sum(arc, is_arc(arc,i,j)+is_arc(arc,j,i)) = c(i,j) ;

   loop(i,
         loop((h,j)$((distance(h,i)<inf) and (distance(i,j)<inf)),
                 distance(h,j) = min(distance(h,j), distance(h,i) + distance(i,j)) ;
         );
   );

alias(arc, arcPrime) ;

parameters
   time_lag_arcs(arc,arcPrime) 2-dimensional time lag
   distance_arcs(arc,arcPrime) distance between arcs
   t_window_arcs(arc,twb)
;

   t_window_arcs(arc,twb) = sum((i,j)$is_arc(arc,i,j), t_window(i,j,twb)) ;
   time_lag_arcs(arc,arcPrime) = sum((g,h,i,j)$(is_arc(arc,g,h) and is_arc(arcPrime,i,j)), c(g,h)+distance(h,i)) ;
   distance_arcs(arc,arcPrime) = sum((g,h,i,j)$(is_arc(arc,g,h) and is_arc(arcPrime,i,j)), distance(h,i)) ;



*############################################################################################################################

variables
   z objective variable ;

positive variables
   start(arc,v) start time of service on arc by vehicle v ;

binary variables
   x(arc,arcPrime,v) arc is traversed by vehicle v
   y(arc,v)          demand is delivered to arc by vehicle v ;

equations
   objective                                objective function
   eq_3_20_node_balance(arc,v)              each traversed arc possesses one predecessor and one successor
   eq_3_21_service_once(arc)                service each arc only once
*   eq_3_21_service_once2(arc,i,arcPrime,h)
   eq_3_22_capacity_constraint(v)           capacity constraint for vehicle v
   eq_3_23_service_cover(arc,v)             serviced arcs must be traversed
   eq_3_28_time_constraint(arc,arcPrime,v)  time window constraints
;

   objective..                                                z =e= sum((arc,arcPrime,v)$(required_arc(arc) and required_arc(arcPrime)), distance_arcs(arc,arcPrime)*x(arc,arcPrime,v)) ;
   eq_3_20_node_balance(arc,v)$required_arc(arc)..            sum(arcPrime$required_arc(arcPrime), x(arcPrime,arc,v)) =e= sum(arcPrime$required_arc(arcPrime), x(arc,arcPrime,v)) ;
   eq_3_21_service_once(arc)$(required_arc(arc) and not depot_arc(arc))..
*                                                             sum(v, y(arc,i,j,v) + y(arcPrime,j,i,v) ) =e= 1 ;
                                                              sum(v, y(arc,v)) =e= 1 ;
   eq_3_22_capacity_constraint(v)..                           sum(arc$required_arc(arc), demand_a(arc)*y(arc,v)) =l= c_veh ;
   eq_3_23_service_cover(arc,v)$required_arc(arc)..           sum(arcPrime$required_arc(arcPrime), x(arcPrime,arc,v)) =g= y(arc,v) ;

   eq_3_28_time_constraint(arc,arcPrime,v)$(required_arc(arc) and required_arc(arcPrime) and not depot_arc(arcPrime))..
                                                              start(arc,v) + time_lag_arcs(arc,arcPrime) - start(arcPrime,v) =l= bigM*(1-x(arc,arcPrime,v)) ;

   start.lo(arc,v) = t_window_arcs(arc,'a') ;
   start.up(arc,v) = t_window_arcs(arc,'b') ;

   x.fx(arc,arc,v) = 0 ;

model carp / all / ;

option mip = cplex ;

carp.optfile = 1 ;

solve carp minimizing z using mip ;

display x.l, y.l ;
