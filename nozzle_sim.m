u = symunit;

%res = acs_simulation();


R = 287.05 * (u.J/ (u.kg * u.K));

%approx values at t+2
m_dot = -0.005 * (u.kg / u.s);
t_t = 210 * (u.K);
p_t = 4 * (u.MPa);
gamma = 1.401;


a_star = ((m_dot * sqrt(t_t)) / p_t ) * ( 1 / (sqrt(gamma / R) * power(((gamma + 1) / 2), -(gamma + 1) / (2 * (gamma - 1))))); 
r_t = sqrt(a_star / pi);
[r_t_val, r_t_u] = separateUnits(r_t);
disp(unitConvert(vpa(r_t_val) * simplify(r_t_u), u.mm))

%r_t = unitConvert(eval(sqrt(a_star / pi)), u.mm);

%disp(r_t)


