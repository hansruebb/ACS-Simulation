u = symunit;

res = sim("acs_simulation");

%User setting
throat_radius = 0.4 * u.mm;
expansion_ratio = 4;
outside_pressure = 1 * (u.bar);

%Values from simulation: (EXECPT GAS CONSTANT AND SPECIFIC HEAT; SHOULD BE CHANGED)
gas_constant = simscape.Value(ones(55, 1) * 287.05, symunit2str(u.J/ (u.kg * u.K), 'Simscape'));
specific_heat_ratio = simscape.Value(ones(55, 1) * 1.401);
mass_flow = simscape.Value(simlog_acs_simulation.Pressure_Regulator.mdot_B.series.values, simlog_acs_simulation.Pressure_Regulator.mdot_B.series.unit);
init_temp = simscape.Value(simlog_acs_simulation.Pressure_Regulator.B.T.series.values, simlog_acs_simulation.Pressure_Regulator.B.T.series.unit);
init_pressure = simscape.Value(simlog_acs_simulation.Pressure_Regulator.B.p.series.values, simlog_acs_simulation.Pressure_Regulator.B.p.series.unit);
time_steps = simlog_acs_simulation.Pressure_Regulator.mdot_B.series.time;

r = gas_constant;
gamma = specific_heat_ratio;
m_dot = -mass_flow;
t_t = init_temp;
p_t = init_pressure;
r_t = simscape.Value(ones(55, 1) * double(getValue(throat_radius)), symunit2str(getUnit(throat_radius)));
ratio = simscape.Value(ones(55, 1) * double(expansion_ratio));
p_0 = simscape.Value(ones(55, 1) * double(getValue(outside_pressure)), symunit2str(getUnit(outside_pressure)));

%Calculate recommended throat size
a_star_recommended = ((m_dot .* sqrt(t_t)) ./ p_t ) .* ( 1 ./ (sqrt(gamma ./ r) .* power(((gamma + 1) ./ 2), -(gamma + 1) ./ (2 .* (gamma - 1))))); 
r_t_recommended = sqrt(a_star_recommended / pi);
%disp(convert(r_t_recommended, 'mm'))

%Calculate thrust for given throat size & expansion ratio
%r_t = r_t_recommended;
syms M_e;
g = double(gamma(1));
exit_mach = vpasolve(expansion_ratio == ((g + 1) / 2)^(-(g + 1) / (2 * (g - 1))) * ((1 + (((g - 1) / 2) * M_e^2))^((g + 1) / (2 * (g - 1))) / M_e), M_e);
exit_mach = simscape.Value(ones(55, 1) * double(exit_mach));
exit_temp = (1 + ((gamma - 1) / 2) .* exit_mach.^2).^(-ones(55, 1)) .* t_t;
exit_pressure = (1 + ((gamma - 1) / 2) .* exit_mach.^2).^(-gamma ./ (gamma - 1)) .* p_t;
v_e = exit_mach .* sqrt(gamma .* r .* exit_temp);
f = m_dot .* v_e + (exit_pressure - p_0) * pi .* r_t.^2 .* ratio;
thrust = convert(f, 'N');

%Plot thrust and output some additional stats
figure
plot(time_steps, value(thrust))
title("Thrust over time")
xlabel("Time (s)")
ylabel("Thrust (N)")

F = griddedInterpolant(time_steps, value(thrust));
f_t = @(t) (F(t));
total_impuls = integral(f_t, time_steps(1), time_steps(end));
disp(['The total impuls is: ', num2str(total_impuls), ' Ns'])

M = griddedInterpolant(time_steps, value(convert(m_dot, 'kg/s')));
m_t = @(t) (M(t));
total_mass = integral(m_t, time_steps(1), time_steps(end));
specific_impuls = total_impuls / (9.8066 * total_mass);
disp(['The specific impuls is: ', num2str(specific_impuls), ' s'])


function v = getValue(expr)
    [v, ~] = separateUnits(expr);
end

function u = getUnit(expr)
    [~, u] = separateUnits(expr);
end