u = symunit;

for r = 0.35:0.05:0.75
    for e = 5:0.2:7
        throat_radius = r * u.mm;
        expansion_ratio = e;
        
        throat_area_raw_mm = double(getValue(unitConvert(throat_radius, u.mm)^2 * pi));
        [tp, isp] = run_sim(throat_radius, expansion_ratio);
        disp(['The specific impuls for radius ', num2str(r), 'mm and expansion ratio ', num2str(e) , ' is: ', num2str(isp), ' s'])
        disp(['The total impuls for radius ', num2str(r), 'mm and expansion ratio ', num2str(e) , ' is: ', num2str(tp), ' Ns'])
    end 
end

function [total_impuls, specific_impuls] = run_sim(throat_radius, expansion_ratio)

    u = symunit;
    
    %Outside pressure, gas constant and specific heat ratio are hard-coded
    %here
    outside_pressure = 1 * (u.bar);
    
    simlog_acs_simulation = sim("acs_simulation").simlog_acs_simulation;

    time_steps = simlog_acs_simulation.Pressure_Regulator.mdot_B.series.time;
    gas_constant = simscape.Value(ones(length(time_steps), 1) * 287.05, symunit2str(u.J/ (u.kg * u.K), 'Simscape'));
    specific_heat_ratio = simscape.Value(ones(length(time_steps), 1) * 1.401);
    mass_flow = simscape.Value(simlog_acs_simulation.Pressure_Regulator.mdot_B.series.values, simlog_acs_simulation.Pressure_Regulator.mdot_B.series.unit);
    init_temp = simscape.Value(simlog_acs_simulation.Pressure_Regulator.B.T.series.values, simlog_acs_simulation.Pressure_Regulator.B.T.series.unit);
    init_pressure = simscape.Value(simlog_acs_simulation.Pressure_Regulator.B.p.series.values, simlog_acs_simulation.Pressure_Regulator.B.p.series.unit);
    
    r = gas_constant;
    gamma = specific_heat_ratio;
    m_dot = -mass_flow;
    t_t = init_temp;
    p_t = init_pressure;
    r_t = simscape.Value(ones(length(time_steps), 1) * double(getValue(throat_radius)), symunit2str(getUnit(throat_radius)));
    ratio = simscape.Value(ones(length(time_steps), 1) * double(expansion_ratio));
    p_0 = simscape.Value(ones(length(time_steps), 1) * double(getValue(outside_pressure)), symunit2str(getUnit(outside_pressure)));
    
    %Calculate recommended throat size
    %a_star_recommended = ((m_dot .* sqrt(t_t)) ./ p_t ) .* ( 1 ./ (sqrt(gamma ./ r) .* power(((gamma + 1) ./ 2), -(gamma + 1) ./ (2 .* (gamma - 1))))); 
    %r_t_recommended = sqrt(a_star_recommended / pi);
    %disp(convert(r_t_recommended, 'mm'))
    %r_t = r_t_recommended;
    
    %Calculate thrust for given throat size & expansion ratio
    syms M_e;
    %TODO Calculate time dependent
    g = double(gamma(1));
    exit_mach = vpasolve(expansion_ratio == ((g + 1) / 2)^(-(g + 1) / (2 * (g - 1))) * ((1 + (((g - 1) / 2) * M_e^2))^((g + 1) / (2 * (g - 1))) / M_e), M_e);
    exit_mach = simscape.Value(ones(length(time_steps), 1) * double(exit_mach));
    t_e = (1 + ((gamma - 1) / 2) .* exit_mach.^2).^(-ones(length(time_steps), 1)) .* t_t;
    p_e = (1 + ((gamma - 1) / 2) .* exit_mach.^2).^(-gamma ./ (gamma - 1)) .* p_t;
    v_e = exit_mach .* sqrt(gamma .* r .* t_e);
    f = m_dot .* v_e + (p_e - p_0) * pi .* r_t.^2 .* ratio;
    thrust = convert(f, 'N');
    
    %Plot thrust and output some additional stats
    %figure
    %plot(time_steps, value(thrust))
    %title("Thrust over time")
    %xlabel("Time (s)")
    %ylabel("Thrust (N)")
    
    F = griddedInterpolant(time_steps, value(thrust));
    f_t = @(t) (F(t));
    total_impuls = integral(f_t, time_steps(1), time_steps(end));
    %disp(['The total impuls is: ', num2str(total_impuls), ' Ns'])
    
    M = griddedInterpolant(time_steps, value(convert(m_dot, 'kg/s')));
    m_t = @(t) (M(t));
    total_mass = integral(m_t, time_steps(1), time_steps(end));
    specific_impuls = total_impuls / (9.8066 * total_mass);
    %disp(['The specific impuls is: ', num2str(specific_impuls), ' s'])

end

function v = getValue(expr)
    [v, ~] = separateUnits(expr);
end

function u = getUnit(expr)
    [~, u] = separateUnits(expr);
end