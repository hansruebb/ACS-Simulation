u = symunit;

%res = acs_simulation();

%approx values at t+2
R = 287.05 * (u.J/ (u.kg * u.K));
m_dot = 0.005 * (u.kg / u.s);
t_t = 210 * (u.K);
p_t = 4 * (u.MPa);
gamma = 1.401;

%Calculate recommended throat size
a_star = ((m_dot * sqrt(t_t)) / p_t ) * ( 1 / (sqrt(gamma / R) * power(((gamma + 1) / 2), -(gamma + 1) / (2 * (gamma - 1))))); 
r_t = sqrt(a_star / pi);
[r_t_val, r_t_u] = separateUnits(r_t);
disp("The recommended throat radius is: ")
disp(unitConvert(vpa(r_t_val) * simplify(r_t_u), u.mm))
disp("")

%Calculate thrust for given throat size & expansion ratio
expansion_ratio = 4;
throat_radius = 0.4 * (u.mm);
p_0 = 1 * (u.bar);
syms M_e;
exit_mach = vpasolve(expansion_ratio == ((gamma + 1) / 2)^(-(gamma + 1) / (2 * (gamma - 1))) * ((1 + (((gamma - 1) / 2) * M_e^2))^((gamma + 1) / (2 * (gamma - 1))) / M_e), M_e);
exit_temp = (1 + ((gamma - 1) / 2) * exit_mach^2)^(-1) * t_t;
exit_pressure = (1 + ((gamma - 1) / 2) * exit_mach^2)^(-gamma / (gamma - 1)) * p_t;
v_e = exit_mach * sqrt(gamma * R * exit_temp);
[v_e_val, v_e_u] = separateUnits(v_e);
exit_velocity = unitConvert(vpa(v_e_val) * simplify(v_e_u), u.m / u.s);
f = m_dot * exit_velocity + (exit_pressure - p_0) * pi * throat_radius^2 * expansion_ratio;
[f_val, f_u] = separateUnits(f);
thrust = unitConvert(vpa(f_val) * simplify(f_u), u.N);
disp("The estimated thrust is: ")
disp(thrust)
disp("")