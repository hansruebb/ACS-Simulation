u = symunit;

load('constants.mat')

r_vec = 0.02:0.02:0.4; %0.05:0.05:1.0;
e_vec = 1.5:0.5:8;
res = zeros(size(r_vec, 2), size(e_vec, 2));

for r_idx = 1:size(r_vec, 2)
    for e_idx = 1:size(e_vec, 2)
        throat_radius = r_vec(r_idx) * u.mm;
        expansion_ratio = e_vec(e_idx);
        
        THROAT_AREA_SQMM = double(getValue(unitConvert(throat_radius, u.mm)^2 * pi));
        [tp, isp] = run_sim(throat_radius, expansion_ratio);
        res(r_idx, e_idx) = isp;
        disp(['The specific impuls for radius ', num2str(r_vec(r_idx)), 'mm and expansion ratio ', num2str(e_vec(e_idx)) , ' is: ', num2str(isp), ' s'])
        disp(['The total impuls for radius ', num2str(r_vec(r_idx)), 'mm and expansion ratio ', num2str(e_vec(e_idx)) , ' is: ', num2str(tp), ' Ns'])
    end 
end

surf(r_vec, e_vec, res');
xlabel("Throat radius (mm)");
ylabel("Expansion ratio");
zlabel("Specific impuls (s)");

function [total_impuls, specific_impuls] = run_sim(throat_radius, expansion_ratio)

    u = symunit;

    simlog_acs_simulation = sim("acs_simulation").simlog_acs_simulation;

    time_steps = simlog_acs_simulation.Nozzle.mdot_A.series.time;                                                                   % Simulation time steps
    r = simscape.Value(ones(length(time_steps), 1) * 287.05, symunit2str(u.J/ (u.kg * u.K), 'Simscape'));                           % Gas constant (hard-coded for now)
    gamma = simscape.Value(ones(length(time_steps), 1) * 1.401);                                                                    % Specific heat ratio (hard-coded for now)
    m_dot = simscape.Value(simlog_acs_simulation.Nozzle.mdot_A.series.values, simlog_acs_simulation.Nozzle.mdot_A.series.unit);     % Input mass flow
    t_t = simscape.Value(simlog_acs_simulation.Nozzle.A.T.series.values, simlog_acs_simulation.Nozzle.A.T.series.unit);             % Input temperature
    p_t = simscape.Value(simlog_acs_simulation.Nozzle.A.p.series.values, simlog_acs_simulation.Nozzle.A.p.series.unit);             % Input pressure
    ratio = simscape.Value(ones(length(time_steps), 1) * double(expansion_ratio));                                                  % Expansion ratio
    p_0 = simscape.Value(ones(length(time_steps), 1) * 1.0, symunit2str(u.bar, 'Simscape'));                                        % Atmospheric pressure (hard-coded for now)
    r_t = simscape.Value(ones(length(time_steps), 1) * double(getValue(throat_radius)), symunit2str(getUnit(throat_radius)));       % Throat radius
    
    %Calculate recommended throat size
    %a_star_recommended = ((m_dot .* sqrt(t_t)) ./ p_t ) .* ( 1 ./ (sqrt(gamma ./ r) .* power(((gamma + 1) ./ 2), -(gamma + 1) ./ (2 .* (gamma - 1))))); 
    %r_t_recommended = sqrt(a_star_recommended / pi);
    %disp(convert(r_t_recommended, 'mm'))
    %r_t = r_t_recommended;
    
    %Calculate thrust for given throat size & expansion ratio
    syms M_e;
    %TODO Calculate time dependent
    g = double(gamma(1));
    exit_mach = vpasolve(expansion_ratio == ((g + 1) / 2)^(-(g + 1) / (2 * (g - 1))) * ((1 + (((g - 1) / 2) * M_e^2))^((g + 1) / (2 * (g - 1))) / M_e), M_e, 3.0);
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

    %disp(['The exit mach is: ', num2str(exit_mach(1).value)])
    
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