%%% Visualize the less conservative method for input constraint inner
%%% approximation. Ref: Nonlinear Model Predictive Control using Feedback 
%%% Linearization and Local Inner Convex Constraint Approximations.
clear all 
clc

import casadi.*

%% Define system dynamics
ode = @(x,u) 1.8 * x + (0.2 * x.^4 + 0.875) * u;
ode_fl =  @(x,u_tilde) u_tilde;
delta_t = 0.4;

% x_SX = SX.sym('x',1);
% u_SX = SX.sym('u',1);

% sys_int = RK4(ode,x_SX,u_SX,delta_t);
% A_dis_SX = jacobian(sys_int,x_SX);
% B_dis_SX = jacobian(sys_int, u_SX);
% % 
% A_dis = full(DM(A_dis_SX));
% B_dis = full(DM(B_dis_SX));
% % 
% disp(A_dis)
% disp(B_dis)

%% Plot constraints
plot_interval = [-2 2 -15 15];
f_constraint1 = @(x,u_tilde) -0.4 * x.^4 + 1.8 * x - 1.75 - u_tilde;
f_constraint2 = @(x,u_tilde) u_tilde -(0.4 * x.^4 + 1.8 * x + 1.75);

x_min = plot_interval(1);
x_max = plot_interval(2);

% hold on; fcontour(@(x,y) 3*x+4*y) ; hold off


syms f1(x) f2(x)
f1(x) =  -0.4 * x.^4 + 1.8 * x - 1.75;
f2(x) =  0.4 * x.^4 + 1.8 * x + 1.75;
Df1 = diff(f1(x),x);
Df2 = diff(f2(x),x);

Df1_inline = inline(Df1, 'x');
Df2_inline = inline(Df2, 'x');

x_eval = [0];
slope1 = Df1_inline(x_eval);
slope2 = Df2_inline(x_eval);

f1_0 = double(f1(x_eval) - slope1 * (x_eval));
f2_0 = double(f2(x_eval) - slope2 * (x_eval));

x_interval = x_min:x_max;
f1_interval = slope1 * (x_interval - x_eval) + f1_0;
f2_interval = slope2 * (x_interval - x_eval) + f2_0;

g_ub_fn = @(x) slope1 * (x - x_eval) + f1_0;
g_lb_fn = @(x) slope2 * (x - x_eval) + f2_0;


u_xmin_ub = double(f2(x_min));
u_xmax_ub = double(f2(x_max));
u_xmin_lb = double(f1(x_min));
u_xmax_lb = double(f1(x_max));

u_xeval_ub = double(f2(x_eval));
u_xeval_lb = double(f1(x_eval));

out_appr_fn_ub1 = @(x) u_xmin_ub + (u_xmin_ub - u_xeval_ub) / (x_min - x_eval) * (x - x_min);
out_appr_fn_ub2 = @(x) u_xmax_ub + (u_xmax_ub - u_xeval_ub) / (x_max - x_eval) * (x - x_max);
out_appr_fn_lb1 = @(x) u_xmin_lb + (u_xmin_lb - u_xeval_lb) / (x_min - x_eval) * (x - x_min);
out_appr_fn_lb2 = @(x) u_xmax_lb + (u_xmax_lb - u_xeval_lb) / (x_max - x_eval) * (x - x_max);

x_l_interval = x_min:x_eval;
x_u_interval = x_eval:x_max;
plot(x_l_interval, out_appr_fn_ub1(x_l_interval), Color='g');hold on;
plot(x_u_interval, out_appr_fn_ub2(x_u_interval), Color='g');
plot(x_l_interval, out_appr_fn_lb1(x_l_interval), Color='g');
plot(x_u_interval, out_appr_fn_lb2(x_u_interval), Color='g');

fimplicit(f_constraint1,'r', plot_interval);
fimplicit(f_constraint2,'r', plot_interval); 


A_g = [slope1,-1;
       -slope2,1;
       1,0;
       -1,0];
b_g = [x_eval * slope1 - f1_0;
       -x_eval * slope2 + f2_0;
       2;
       2];
P_g = polytope(A_g,b_g);
plot(P_g);
% plot(x_interval,f1_interval);
% plot(x_interval,f2_interval);


%% Outer approximation of the reachable set
ode_fl =  @(x,u_tilde) u_tilde;
delta_t = 0.4;

x_SX = SX.sym('x',1);
u_SX = SX.sym('u',1);

sys_fl_int = RK4(ode_fl,x_SX,u_SX,delta_t);
A_fl_dis_SX = jacobian(sys_fl_int,x_SX);
B_fl_dis_SX = jacobian(sys_fl_int, u_SX);
% % 
A_dis_fl = full(DM(A_fl_dis_SX));
B_dis_fl = full(DM(B_fl_dis_SX));
% % 
disp(A_dis_fl);
disp(B_dis_fl);

pi_constraint_lb = @(x) -0.4 * x.^4 + 1.8 * x - 1.75;
pi_constraint_ub = @(x) 0.4 * x.^4 + 1.8 * x + 1.75;

x_init = 1.9;

O_0_lb = pi_constraint_lb(x_init);
O_0_ub = pi_constraint_ub(x_init);

plot([x_init, x_init],[O_0_lb, O_0_ub]);

X1_ub = A_dis_fl * x_init + B_dis_fl * O_0_ub;
X1_lb = A_dis_fl * x_init + B_dis_fl * O_0_lb;

x_interval_1 = X1_lb :0.01: x_max;

O_1_ub = out_appr_fn_ub2(x_interval_1);
O_1_lb = out_appr_fn_lb2(x_interval_1);

G1 = g_ub_fn(X1_lb);
G2 = g_lb_fn(X1_lb);

slope1_g1 = 1.0;
slope2_g1 = 2.0;
A_g1 = [slope1_g1,-1;
       -slope2_g1,1;
       1,0;
       -1,0];
b_g1 = [X1_lb * slope1_g1 - G1;
       -X1_lb * slope2_g1 + G2;
       x_max;
       -X1_lb];
P_g1 = polytope(A_g1,b_g1);
plot(P_g1)


X2_ub = A_dis_fl * x_interval_1 + B_dis_fl * O_1_ub;
X2_lb = B_dis_fl * x_interval_1 + B_dis_fl * O_1_lb;

x2_min = min(min(X2_ub), min(X2_lb));
x2_max = max(max(X2_ub), max(X2_lb));

