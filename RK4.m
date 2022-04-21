function [x_next] =  RK4(f,x,u,delta_t)

k1 = f(x,u);
k2 = f(x + delta_t / 2 * k1, u);
k3 = f(x + delta_t / 2 * k2, u);
k4 = f(x +delta_t * k3, u);

x_next = x + delta_t / 6 * (k1 + 2 * k2 +2 * k3 +k4);