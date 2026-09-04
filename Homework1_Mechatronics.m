% Joint Locations
A = [1.4 0.485 0];
B = [1.67 0.99 0];
C = [0.255 1.035 0];
D = [0.285 0.055 0];
E = [0.195 2.54 0];
F = [-0.98 2.57 0];
G = [0.05 0.2 0];
%Link Vectors
AB = B - A;
BC = C - B;
CD = D - C;
CDE = E - D;
EF = F - E;
FG = G - F;
%Link lengths
lAB = norm(AB);
lBC = norm(BC);
lCD = norm(CD);
lCDE = norm(CDE);
lEF = norm(EF);
lFG = norm(FG);
% Center of Mass of each link
S1 = (A + B)/2;
S2 = (B + C)/2;
S3 = (D + E)/2;
S4 = (E + F)/2;
S5 = (F + G)/2;
% Weights of each link
WAB = [0 -1 0];
WBC = [0 -1 0];
WCD = [0 -1 0];
WCDE = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];
%Force of each link
syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin
F_A = [FAx FAy 0];
F_B = [FBx FBy 0];
F_C = [FCx FCy 0];
F_D = [FDx FDy 0];
F_E = [FEx FEy 0];
F_F = [FFx FFy 0];
F_G = [FGx FGy 0];
T_in = [0 0 Tin];

% Angular Velocities Calculations

% Loop ABCDA

syms wBC wCDE

omega_AB = [0 0 1]; % WHY IS THIS 1?
omega_BC = [0 0 wBC];
omega_CDE = [0 0 wCDE];

eqn11 = cross(omega_AB, B - A) + cross(omega_BC, C - B) + cross(omega_CDE, D - C) == [0 0 0];

loop1Solution = solve(eqn11, [wBC wCDE]);

angularVelocity_BC = double(loop1Solution.wBC);
angularVelocity_CDE = double(loop1Solution.wCDE);

% Second Loop GFEDG

syms wGF wFE 

omega_CDE = [0 0 angularVelocity_CDE];

omega_BC = [0 0 angularVelocity_BC];

omega_GF = [0 0 wGF];

omega_FE = [0 0 wFE];

eqn12 = cross(omega_GF, F - G) + cross(omega_FE, E - F) + cross(omega_CDE, D - E) == [0 0 0];

loop2Solution = solve(eqn12, [wGF wFE]);

angularVelocity_GF = double(loop2Solution.wGF);

angularVelocity_FE = double(loop2Solution.wFE);

% Circle Intersections Technique

initial_theta = atan2(B(2)-A(2), B(1) - A(1));

if (initial_theta < 0)
    inputAngle = 2 * pi + initial_theta;
else
    inputAngle = initial_theta;
end

for theta =  1:1:360

    % New position of Joint B
    B_new = A + [lAB*cos(inputAngle + deg2rad(theta)) lAB*sin(inputAngle + deg2rad(theta)) 0];

    % New position of Joint C
    % with B_new as center, BC as radius
    % with D as center and DC as radius

    [C_x, C_y] = circcirc(B_new(1), B_new(2), lBC, D(1), D(2), lCD);

    if any(isnan(C_x))
        fprintf('Mechanism reaches its limit at theta = %d degrees\n', theta-1);
        break
    end

    circIntersectC_x = any(isnan(vpa(C_x)));
    circIntersectC_y = any(isnan(vpa(C_y)));

    if circIntersectC_y == 0 && circIntersectC_x == 0

        C_1 = [C_x(1) C_y(1) 0];
        C_2 = [C_x(2) C_y(2) 0];
        
        % Distance to determine which one of the solutions is correct
        dist1 = norm(C_1 - C);
        dist2 = norm(C_2 - C);
        if dist1 < dist2
            C_new = C_1;
        else
            C_new = C_2;
        end
    else
        fprintf('New position C cannot be determined at angle: %d degree', theta);
    end

    % New position of E
    % E is co-linear with C so we can use that to calculate it

    E_new = D + lCDE/lCD * (C_new - D);

    % New position of Joint F
    % with E_new as center, EF as radius
    % with G as center and GF as radius

    [Fx, Fy] = circcirc(E_new(1), E_new(2), lEF, G(1), G(2), lFG);

    % Checking if there is a NaN

    circIntersectF_x = any(isnan(vpa(Fx)));
    circIntersectF_y = any(isnan(vpa(Fy)));

    if circIntersectF_y == 0 && circIntersectF_x == 0
        F_1 = [Fx(1) Fy(1) 0];
        F_2 = [Fx(2) Fy(2) 0];

        % Determine which one of the solutions is correct for F
        distF1 = norm(F_1 - F);
        distF2 = norm(F_2 - F);
        if distF1 < distF2
            F_new = F_1;
        else
            F_new = F_2;
        end
    end

    % Store values for plotting

    new_B_x(theta) = B_new(1);
    new_B_y(theta) = B_new(2);
    new_C_x(theta) = C_new(1);
    new_C_y(theta) = C_new(2);
    new_E_x(theta) = E_new(1);
    new_E_y(theta) = E_new(2);
    new_F_x(theta) = F_new(1);
    new_F_y(theta) = F_new(2);

    B = B_new;
    C = C_new;
    E = E_new;
    F = F_new;



end

figure;            
hold on;
grid on;
axis equal;

plot(new_B_x, new_B_y, 'b-', 'LineWidth', 1.5);
plot(new_C_x, new_C_y, 'r-', 'LineWidth', 1.5);
plot(new_E_x, new_E_y, 'g-', 'LineWidth', 1.5);
plot(new_F_x, new_F_y, 'm-', 'LineWidth', 1.5);

xlabel('X Position');
ylabel('Y Position');
title('Joint Trajectories');
legend('B', 'C', 'E', 'F');

hold off;
