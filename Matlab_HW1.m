
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

WAB = [0 -237.8 0];
WBC = [0 -575.1 0];
WCDE = [0 -992.1 0];
WEF = [0 -479.0 0];
WFG = [0 -1042.4 0];


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


%% Static Solutions

ArtifactWeight = [0 -200 0];


%Equilibrium of each link

%AB
eqn1 = F_A + F_B + WAB == 0;
eqn2 = cross((A-S1), F_A) + cross(B-S1, F_B) + T_in == 0;
%BC
eqn3 = -F_B + F_C + WBC == 0;
eqn4 = cross((B-S2), -F_B) + cross((C-S2), F_C) == 0;
%CDE
eqn5 = -F_C + F_D + F_E + WCDE == 0;
eqn6 = cross((C-S3), -F_C) + cross((D-S3), F_D) + cross((E-S3), F_E) == 0;
%FE
eqn7 = F_F - F_E + WEF + ArtifactWeight == 0;
eqn8 = cross((F-S4), F_F) + cross((E-S4), -F_E) == 0;
%FG
eqn9 = -F_F + F_G + WFG == 0; 
eqn10 = cross((F-S5), -F_F) + cross((G-S5), F_G) == 0;

eqnMatrix = [eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn9 eqn10];
StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

ForceAx = double(StaticSolution.FAx)
ForceAy = double(StaticSolution.FAy)
ForceBx = double(StaticSolution.FBx)
ForceBy = double(StaticSolution.FBy)
ForceCx = double(StaticSolution.FCx)
ForceCy = double(StaticSolution.FCy)
ForceDx = double(StaticSolution.FDx)
ForceDy = double(StaticSolution.FDy)
ForceEx = double(StaticSolution.FEx)
ForceEy = double(StaticSolution.FEy)
ForceFx = double(StaticSolution.FFx)
ForceFy = double(StaticSolution.FFy)
ForceGx = double(StaticSolution.FGx)
ForceGy = double(StaticSolution.FGy)
InputTorque = double(StaticSolution.Tin)

%% Newton's Second Law

% Angular Velocities Calculations

% Loop ABCDA
syms wBC wCDE
omega_AB = [0 0 ((12500*(2*pi))/(9*3600)) ]; % WHY IS THIS 1?
omega_BC = [0 0 wBC];
omega_CDE = [0 0 wCDE];

eqn11 = cross(omega_AB, B - A) + cross(omega_BC, C - B) + cross(omega_CDE, D - C) == [0 0 0];
loop1Solution = solve(eqn11, [wBC wCDE]);
angularvelocityBC = double(loop1Solution.wBC);
angularvelocityCDE = double(loop1Solution.wCDE);

% Second Loop GFEDG
syms wGF wFE
omega_CDE = [0 0 angularvelocityCDE];
omega_BC = [0 0 angularvelocityBC];
omega_GF = [0 0 wGF];
omega_FE = [0 0 wFE];
eqn12 = cross(omega_GF, F - G) + cross(omega_FE, E - F) + cross(omega_CDE, D - E);
loop2Solution = solve(eqn12, [wGF wFE]);
angularvelocityBC = double(loop2Solution.wGF);
angularvelocityCDE = double(loop2Solution.wFE);

% Velocity at Joint A = 0
% Velocity at Joint B
vB = cross(omega_AB, B - A);
% Velocity at Joint C
% vC = vC_A = vC_B + vB_A
vC_B = cross(omega_BC, C - B);
vC = vC_B + vB;
% Velocity at Joint D = 0
% Velocity at Joint E
% vE = vE_D = vE_A because D and A are both grounded
vE = cross(omega_CDE, E -D);
% Velocity at Joint F
% vF = vF_G = vF_A because G and A are both grounded
vF = cross(omega_GF, F - G);
% Velocity at Joint G = 0
% Velocities of Center of Masses
% Velocity of S1
% vS1 = vS1_A
vS1 = cross(omega_AB, S1 - A);
% Velocity of S2
% vS2 = vS2_A = vS2_B + vB
vS2 = cross(omega_AB, S2 - A) + vB;
% Velocity of S3
% vS3 = vS3_D
vS3 = cross(omega_CDE, S3 - D);
% Velocity of S4
% vS4 = vS4_D = vS4_E + vE
vS4 = cross(omega_CDE, S4 - E) + vE;
% Velocity of S5
% vS5 = vS5_G
vS5 = cross(omega_GF, S5 - G);


%Angular Accelerations Calculations

syms aBC aCDE
alpha_AB = [0 0 0]; % why is this 0?
alphaBC = [0 0 aBC];
alphaCDE = [0 0 aCDE];
a_B_A = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));
a_C_B = cross(alphaBC, C-B) + cross(omega_BC, cross(omega_BC, C-B));
a_D_C = cross(alphaCDE, D-C) + cross(omega_CDE, cross(omega_CDE, D-C));
eqn13 = a_B_A + a_C_B + a_D_C == 0;
loop1AccSolution = solve(eqn13, [aBC aCDE]);
alpha_BC = double(loop1AccSolution.aBC);
alpha_CDE = double(loop1AccSolution.aCDE);

% Loop 2 velocity (solve for omega_GF, omega_FE) -- needed before loop 2 accel
syms wGF wFE
omega_GF = [0 0 wGF];
omega_FE = [0 0 wFE];
eqn12 = cross(omega_GF, F-G) + cross(omega_FE, E-F) + cross(omega_CDE, D-E) == 0;
loop2VelSolution = solve(eqn12, [wGF wFE]);
omega_GF = [0 0 double(loop2VelSolution.wGF)];
omega_FE = [0 0 double(loop2VelSolution.wFE)];

% Loop 2 acceleration (solve for alphaGF, alphaFE)
syms aGF aFE
alphaGF = [0 0 aGF];
alphaFE = [0 0 aFE];
alphaCDE_vector = [0 0 alpha_CDE];
a_F_G = cross(alphaGF, F-G) + cross(omega_GF, cross(omega_GF, F-G));
a_E_F = cross(alphaFE, E-F) + cross(omega_FE, cross(omega_FE, E-F));
a_D_E = cross(alphaCDE_vector, D-E) + cross(omega_CDE, cross(omega_CDE, D-E));
eqn14 = a_F_G + a_E_F + a_D_E == 0;
loop2AccSolution = solve(eqn14, [aGF aFE]);
alpha_GF = double(loop2AccSolution.aGF);
alpha_FE = double(loop2AccSolution.aFE);


% Accelerations Calculations
alphaBC_vector = [0 0 alpha_BC];
alphaCDE_vector = [0 0 alpha_CDE];
alphaFE_vector = [0 0 alpha_FE];
alphaGF_vector = [0 0 alpha_GF];


%Joints - reference is second letter, i.e. a_BA is acceleration of joint B in reference to a
a_BA = cross(alpha_AB, B-A) + cross(omega_AB, (cross(omega_AB, B-A)));
a_CB = cross(alphaBC_vector, C-B) + cross(omega_BC, (cross(omega_BC, C-B)));
a_CA = a_CB + a_BA;
a_ED = cross(alphaCDE_vector, E-D) + cross(omega_CDE, (cross(omega_CDE, E-D)));
a_FG = cross(alphaGF_vector, F-G) + cross(omega_GF, cross(omega_GF, F-G));

%Centers of Mass, each respective to a grounded pin

a_s1_A = cross(alpha_AB, S1-A) + cross(omega_AB, cross(omega_AB, S1-A));
a_s2_A = a_BA + cross(alphaBC_vector, S2-B) + cross(omega_BC, cross(omega_BC, S2-B));
a_s3_D = cross(alphaCDE_vector, S3-D) + cross(omega_CDE, cross(omega_CDE, S3-D));
a_s4_D = a_ED + cross(alphaFE_vector, S4-E) + cross(omega_FE, cross(omega_FE, S4-E));
a_s5_G = cross(alphaGF_vector, S5-G) + cross(omega_GF, cross(omega_GF, S5-G));

% Masses and Moments of Inertia

MassAB = 23.78;
MassBC = 57.51;
MassCDE = 99.21;
MassEF = 47.90;
MassGF = 104.24;


% Izz --> ASK IF THIS IS THE CORRECT MOMENT OF INERTIA
J_AB = 90.62;
J_BC = 156.87;
J_CDE = 226.28;
J_EF = 325.88;
J_FG = 281.78;

%Forces

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NInputTorque = [0 0 NTin];


% Link AB
Neqn1 = NForceA + NForceB + WAB == MassAB * a_s1_A;
Neqn2 = cross(A-S1, NForceA) + cross(B-S1, NForceB) + NInputTorque == J_AB*alpha_AB;
% Link BC
Neqn3 = -NForceB + NForceC + WBC == MassBC * a_s2_A;
Neqn4 = cross(B-S2, -NForceB) + cross(C-S2, NForceC) == J_BC * alphaBC_vector;
% Link CDE
Neqn5 = -NForceC + NForceD + NForceE + WCDE == MassCDE * a_s3_D;
Neqn6 = cross(C-S3, -NForceC) + cross(D-S3, NForceD) + cross(E-S3, NForceE) == J_CDE * alphaCDE_vector;
% Link EF
Neqn7 = -NForceE + NForceF + WEF + ArtifactWeight == MassEF * a_s4_D;
Neqn8 = cross(E-S4, -NForceE) + cross(F-S4, NForceF) == J_EF * alphaFE_vector;
% Link FG
Neqn9 = -NForceF + NForceG + WFG == MassGF * a_s5_G;
Neqn10 = cross( F-S5, -NForceF) + cross(G-S5, NForceG) == J_FG * alphaGF_vector;


%Solution
NEqnMatrix = [Neqn1 Neqn2 Neqn3 Neqn4 Neqn5 Neqn6 Neqn7 Neqn8 Neqn9 Neqn10];
DynamicSolution = solve(NEqnMatrix, [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin]);

%Values
NForce_Ax = double(DynamicSolution.NFAx)
NForce_Ay = double(DynamicSolution.NFAy)
NForce_Bx = double(DynamicSolution.NFBx)
NForce_By = double(DynamicSolution.NFBy)
NForce_Cx = double(DynamicSolution.NFCx)
NForce_Cy = double(DynamicSolution.NFCy)
NForce_Dx = double(DynamicSolution.NFDx)
NForce_Dy = double(DynamicSolution.NFDy)
NForce_Ex = double(DynamicSolution.NFEx)
NForce_Ey = double(DynamicSolution.NFEy)
NForce_Fx = double(DynamicSolution.NFFx)
NForce_Fy = double(DynamicSolution.NFFy)
NForce_Gx = double(DynamicSolution.NFGx)
NForce_Gy = double(DynamicSolution.NFGy)
NInputTorque = double(DynamicSolution.NTin)

%% Comparisons

deltaAx = NForce_Ax - ForceAx;
deltaAy = NForce_Ay - ForceAy;
deltaBx = NForce_Bx - ForceBx;
deltaBy = NForce_By - ForceBy;
deltaCx = NForce_Cx - ForceCx;
deltaCy = NForce_Cy - ForceCy;
deltaDx = NForce_Dx - ForceDx;
deltaDy = NForce_Dy - ForceDy;
deltaFx = NForce_Fx - ForceFx;
deltaFy = NForce_Fy - ForceFy;
deltaGx = NForce_Gx - ForceGx;
deltaGy = NForce_Gy - ForceGy;
deltaTin = NInputTorque - InputTorque;

%% Circle Intersections Technique

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