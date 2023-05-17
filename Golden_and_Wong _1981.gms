$ontext
Model of Golden and Wong (1981) with compact subtour elimination constraints
$offText

option resLim  = 200;

* $call g********e instance_excel.xlsx o=instance_excel.gdx index=Parameter!A1:D7

* ******* instance_excel.gdx

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
;

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
   q(i,j)          demand on arc i-j
;

   c(i,j) = sum(arc$(ord(i)=ord_i(arc) and ord(j)=ord_j(arc)), cost_a(arc)) ;
   c('i%nn%',i)$(ord(i)>1) = c('i1',i) ;
   c(i,'i%nn%') = c(i,'i1') ;
   c('i%nn%','i%nn%') = 0 ;

display c, i,j, arc ;

** Assume symmetric network
   c(i,j)$(c(i,j)=0) = c(j,i) ;

  loop(arc,
    q(i,j)$(ord(i)=ord_i(arc) and ord(j)=ord_j(arc) and demand_a(arc)) = demand_a(arc);
    );

alias(v,k);

q(i,j) = q(i,j) + q(j,i);

display q,c;


VARIABLES
z              objective Variable ;
Positive Variable
f(i,j,k)  flow variable;

Binary Variable
visit(i,j,k)     i-j is traversed by postman p
Service(i,j,k)   i-j is serviced by postman p
;   

EQUATIONS
obj the objective function of the problem

e_2_flow_balance(i,k)               Route continuity
e_3_flow_balance(i,j)
e_4_flow_balance(i,j,k)             Illegal and legal subtours
e_5_Kapazitaet_Bedingung(k)         Vehicle capacity is not violated when servicing the arc.
e_6_1_subtour_break(i,k)
e_6_2_subtour_break(i,j,k)
e_6_3_lowerbound(i,j,k)
;

obj                                     .. Z    =   E   =   sum(i, sum(j$c(i,j),sum(k, c(i,j)*visit(i,j,k))));
e_2_flow_balance(i,k)                   ..sum(j$c(j,i),visit(j,i,k))=E=sum(j$c(i,j),visit(i,j,k))  ;
e_3_flow_balance(i,j)$c(i,j)            ..sum(k, Service(i,j,k) + Service(j,i,k))=e=ceil(q(i,j)/c_veh);
e_4_flow_balance(i,j,k)$c(i,j)          .. visit(i,j,k) =g= Service(i,j,k);
e_5_Kapazitaet_Bedingung(k)             .. sum(i,sum(j$q(i,j), q(i,j)*Service(i,j,k)))=L=c_veh;

e_6_1_subtour_break(i,k)$(ord(i)<>1)    .. sum(j$c(i,j),f(i,j,k)) - sum(j$c(j,i),f(j,i,k)) =E=sum(j$c(i,j),Service(i,j,k));
e_6_2_subtour_break(i,j,k)$c(i,j)       .. f(i,j,k)=L=card(i)*card(i)*visit(i,j,k);
e_6_3_lowerbound(i,j,k)$c(i,j)          .. f(i,j,k) =g= 0;

MODEL Tours  /all/;

SOLVE Tours min Z using MIP;

option MIP = cplex

display visit.l, service.l;
*---end input file--- 
