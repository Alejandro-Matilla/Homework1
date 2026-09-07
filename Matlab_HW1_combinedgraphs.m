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


ForceAx = double(StaticSolution.FAx);
ForceAy = double(StaticSolution.FAy);
ForceBx = double(StaticSolution.FBx);
ForceBy = double(StaticSolution.FBy);
ForceCx = double(StaticSolution.FCx);
ForceCy = double(StaticSolution.FCy);
ForceDx = double(StaticSolution.FDx);
ForceDy = double(StaticSolution.FDy);
ForceEx = double(StaticSolution.FEx);
ForceEy = double(StaticSolution.FEy);
ForceFx = double(StaticSolution.FFx);
ForceFy = double(StaticSolution.FFy);
ForceGx = double(StaticSolution.FGx);
ForceGy = double(StaticSolution.FGy);
InputTorque = double(StaticSolution.Tin);

%% Newton's Second Law

% Angular Velocities Calculations

% Loop ABCDA
syms wBC wCDE
omega_AB = [0 0 ((12500*(2*pi))/(9*3600)) ]; 
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
NForce_Ax = double(DynamicSolution.NFAx);
NForce_Ay = double(DynamicSolution.NFAy);
NForce_Bx = double(DynamicSolution.NFBx);
NForce_By = double(DynamicSolution.NFBy);
NForce_Cx = double(DynamicSolution.NFCx);
NForce_Cy = double(DynamicSolution.NFCy);
NForce_Dx = double(DynamicSolution.NFDx);
NForce_Dy = double(DynamicSolution.NFDy);
NForce_Ex = double(DynamicSolution.NFEx);
NForce_Ey = double(DynamicSolution.NFEy);
NForce_Fx = double(DynamicSolution.NFFx);
NForce_Fy = double(DynamicSolution.NFFy);
NForce_Gx = double(DynamicSolution.NFGx);
NForce_Gy = double(DynamicSolution.NFGy);
NInputTorque = double(DynamicSolution.NTin);

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

%Static Preallocation
N = 360;
ForceAx_l_all = zeros(1,N); 
ForceAy_l_all = zeros(1,N);
ForceBx_l_all = zeros(1,N);
ForceBy_l_all = zeros(1,N);
ForceCx_l_all = zeros(1,N); 
ForceCy_l_all = zeros(1,N);
ForceDx_l_all = zeros(1,N); 
ForceDy_l_all = zeros(1,N);
ForceEx_l_all = zeros(1,N);
ForceEy_l_all = zeros(1,N);
ForceFx_l_all = zeros(1,N); 
ForceFy_l_all = zeros(1,N);
ForceGx_l_all = zeros(1,N); 
ForceGy_l_all = zeros(1,N);
InputTorque_l_all = zeros(1,N);

%syms for Static Force Solution
syms FAx_l FAy_l FBx_l FBy_l FCx_l FCy_l FDx_l FDy_l FEx_l FEy_l FFx_l FFy_l FGx_l FGy_l Tin_l

% static forces vector preallocation

F_A_loop = [FAx_l FAy_l 0];
F_B_loop = [FBx_l FBy_l 0];
F_C_loop = [FCx_l FCy_l 0];
F_D_loop = [FDx_l FDy_l 0];
F_E_loop = [FEx_l FEy_l 0];
F_F_loop = [FFx_l FFy_l 0];
F_G_loop = [FGx_l FGy_l 0];
T_in_loop = [0 0 Tin_l];

% Dynamic Solution Preallocation

NForceAx_l_all = zeros(1,N); 
NForceAy_l_all = zeros(1,N);
NForceBx_l_all = zeros(1,N);
NForceBy_l_all = zeros(1,N);
NForceCx_l_all = zeros(1,N); 
NForceCy_l_all = zeros(1,N);
NForceDx_l_all = zeros(1,N); 
NForceDy_l_all = zeros(1,N);
NForceEx_l_all = zeros(1,N);
NForceEy_l_all = zeros(1,N);
NForceFx_l_all = zeros(1,N); 
NForceFy_l_all = zeros(1,N);
NForceGx_l_all = zeros(1,N); 
NForceGy_l_all = zeros(1,N);
NInputTorque = zeros(1,N);


%syms for Dynamic Solution

syms NFAx_l NFAy_l NFBx_l NFBy_l NFCx_l NFCy_l NFDx_l NFDy_l NFEx_l NFEy_l NFFx_l NFFy_l NFGx_l NFGy_l NTin_l

% dynamic forces vector preallocation

NF_A_loop = [NFAx_l NFAy_l 0];
NF_B_loop = [NFBx_l NFBy_l 0];
NF_C_loop = [NFCx_l NFCy_l 0];
NF_D_loop = [NFDx_l NFDy_l 0];
NF_E_loop = [NFEx_l NFEy_l 0];
NF_F_loop = [NFFx_l NFFy_l 0];
NF_G_loop = [NFGx_l NFGy_l 0];
NT_in_loop = [0 0 NTin_l];

%starting theta position

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

% New COMs
S1_new = (A + B_new)/2;
S2_new = (B_new + C_new)/2;
S3_new = (D + E_new)/2;
S4_new = (E_new + F_new)/2;
S5_new = (F_new + G)/2;

%New Static Solution

%Equilibrium of each link
%AB
Loopeqn1 = F_A_loop + F_B_loop + WAB == 0;
Loopeqn2 = cross((A-S1_new), F_A_loop) + cross(B_new-S1_new, F_B_loop) + T_in_loop == 0;
%BC
Loopeqn3 = -F_B_loop + F_C_loop + WBC == 0;
Loopeqn4 = cross((B_new-S2_new), -F_B_loop) + cross((C_new-S2_new), F_C_loop) == 0;
%CDE
Loopeqn5 = -F_C_loop + F_D_loop + F_E_loop + WCDE == 0;
Loopeqn6 = cross((C_new-S3_new), -F_C_loop) + cross((D-S3_new), F_D_loop) + cross((E_new-S3_new), F_E_loop) == 0;
%FE
Loopeqn7 = F_F_loop - F_E_loop + WEF + ArtifactWeight == 0;
Loopeqn8 = cross((F_new-S4_new), F_F_loop) + cross((E_new-S4_new), -F_E_loop) == 0;
%FG
Loopeqn9 = -F_F_loop + F_G_loop + WFG == 0;
Loopeqn10 = cross((F_new-S5_new), -F_F_loop) + cross((G-S5_new), F_G_loop) == 0;


LoopeqnMatrix = [Loopeqn1 Loopeqn2 Loopeqn3 Loopeqn4 Loopeqn5 Loopeqn6 Loopeqn7 Loopeqn8 Loopeqn9 Loopeqn10];
LoopStaticSolution = solve(LoopeqnMatrix, [FAx_l FAy_l FBx_l FBy_l FCx_l FCy_l FDx_l FDy_l FEx_l FEy_l FFx_l FFy_l FGx_l FGy_l Tin_l]);


%get each force value for each theta (1-360), store for later plotting
ForceAx_l_all(theta) = LoopStaticSolution.FAx_l;
ForceAy_l_all(theta) = LoopStaticSolution.FAy_l;
ForceBx_l_all(theta) = LoopStaticSolution.FBx_l;
ForceBy_l_all(theta) = LoopStaticSolution.FBy_l;
ForceCx_l_all(theta) = LoopStaticSolution.FCx_l;
ForceCy_l_all(theta) = LoopStaticSolution.FCy_l;
ForceDx_l_all(theta) = LoopStaticSolution.FDx_l;
ForceDy_l_all(theta) = LoopStaticSolution.FDy_l;
ForceEx_l_all(theta) = LoopStaticSolution.FEx_l;
ForceEy_l_all(theta) = LoopStaticSolution.FEy_l;
ForceFx_l_all(theta) = LoopStaticSolution.FFx_l;
ForceFy_l_all(theta) = LoopStaticSolution.FFy_l;
ForceGx_l_all(theta) = LoopStaticSolution.FGx_l;
ForceGy_l_all(theta) = LoopStaticSolution.FGy_l;
InputTorque_l_all(theta) = LoopStaticSolution.Tin_l;

% Calculating angular velocities for links at every theta
    
    syms wBC_new wCDE_new

    omega_BC_new = [0 0 wBC_new];
    omega_CDE_new = [0 0 wCDE_new];

    %omega_AB is constant so no need to recalculate

    eqn18 = cross(omega_AB, B_new - A) + cross(omega_BC_new, C_new - B_new) + cross(omega_CDE_new, D - C_new) == [0 0 0];

    loop1Solution = solve(eqn18, [wBC_new wCDE_new]);

    angularVelocity_BC_new = double(loop1Solution.wBC_new);
    angularVelocity_CDE_new = double(loop1Solution.wCDE_new);

    % Second Loop GFEDG

    syms wGF_new wFE_new

    omega_CDE_new = [0 0 angularVelocity_CDE_new];

    omega_BC_new = [0 0 angularVelocity_BC_new];

    omega_GF_new = [0 0 wGF_new];

    omega_FE_new = [0 0 wFE_new];

    eqn19 = cross(omega_GF_new, F_new - G) + cross(omega_FE_new, E_new - F_new) + cross(omega_CDE_new, D - E_new) == [0 0 0];

    loop2Solution = solve(eqn19, [wGF_new wFE_new]);

    angularVelocity_GF_new = double(loop2Solution.wGF_new);

    angularVelocity_FE_new = double(loop2Solution.wFE_new);

    omega_GF_new = [0 0 angularVelocity_GF_new];
    omega_FE_new = [0 0 angularVelocity_FE_new];

    % Calculating Velocities for joints at each angle of theta
    vB_new = cross(omega_AB, B_new - A);

    vC_B_new = cross(omega_BC_new, C_new - B_new);
    vC_new = vC_B_new + vB_new;
    vE_new = cross(omega_CDE_new, E_new - D);
    vF_new = cross(omega_GF_new, F_new - G);

    %Angular Accelerations Calculations at each angle

    syms aBC_new aCDE_new
    alphaBC_new = [0 0 aBC_new];
    alphaCDE_new = [0 0 aCDE_new];
    a_B_A_new = cross(alpha_AB, B_new-A) + cross(omega_AB, cross(omega_AB, B_new-A));
    a_C_B_new = cross(alphaBC_new, C_new-B_new) + cross(omega_BC_new, cross(omega_BC_new, C_new-B_new));
    a_D_C_new = cross(alphaCDE_new, D-C_new) + cross(omega_CDE_new, cross(omega_CDE_new, D-C_new));
    eqn20 = a_B_A_new + a_C_B_new + a_D_C_new == 0;
    loop1AccSolution = solve(eqn20, [aBC_new aCDE_new]);
    alpha_BC_new = double(loop1AccSolution.aBC_new);
    alpha_CDE_new = double(loop1AccSolution.aCDE_new);

    % Loop 2 acceleration (solve for alphaGF, alphaFE)
    syms aGF_new aFE_new
    alphaGF_new = [0 0 aGF_new];
    alphaFE_new = [0 0 aFE_new];
    alphaCDE_vector_new = [0 0 alpha_CDE_new];
    a_F_G_new = cross(alphaGF_new, F_new - G) + cross(omega_GF_new, cross(omega_GF_new, F_new - G));
    a_E_F_new = cross(alphaFE_new, E_new - F_new) + cross(omega_FE_new, cross(omega_FE_new, E_new - F_new));
    a_D_E_new = cross(alphaCDE_vector_new, D - E_new) + cross(omega_CDE_new, cross(omega_CDE_new, D - E_new));
    eqn21 = a_F_G_new + a_E_F_new + a_D_E_new == 0;
    loop2AccSolution = solve(eqn21, [aGF_new aFE_new]);
    alpha_GF_new = double(loop2AccSolution.aGF_new);
    alpha_FE_new = double(loop2AccSolution.aFE_new);


    % Accelerations Calculations
    alphaBC_vector_new = [0 0 alpha_BC_new];
    alphaCDE_vector_new = [0 0 alpha_CDE_new];
    alphaFE_vector_new = [0 0 alpha_FE_new];
    alphaGF_vector_new = [0 0 alpha_GF_new];



    %Joints - 
    a_BA_new = cross(alpha_AB, B_new - A) + cross(omega_AB, (cross(omega_AB, B_new - A)));
    a_CB_new = cross(alphaBC_vector_new, C_new - B_new) + cross(omega_BC_new, (cross(omega_BC_new, C_new - B_new)));
    a_CA_new = a_CB_new + a_BA_new;
    a_ED_new = cross(alphaCDE_vector_new, E_new - D) + cross(omega_CDE_new, (cross(omega_CDE_new, E_new - D)));
    a_FG_new = cross(alphaGF_vector_new, F_new - G) + cross(omega_GF_new, cross(omega_GF_new, F_new - G));

    % Joint, Centers of Mass
    a_s1_A_l = cross(alpha_AB, S1_new-A) + cross(omega_AB, cross(omega_AB, S1_new-A));
    a_s2_A_l = a_BA_new + cross(alphaBC_vector_new, S2_new-B_new) + cross(omega_BC_new, cross(omega_BC_new, S2_new-B_new));
    a_s3_D_l = cross(alphaCDE_vector_new, S3_new-D) + cross(omega_CDE_new, cross(omega_CDE_new, S3_new-D));
    a_s4_D_l = a_ED_new + cross(alphaFE_vector_new, S4_new-E_new) + cross(omega_FE_new, cross(omega_FE_new, S4_new-E_new));
    a_s5_G_l = cross(alphaGF_vector_new, S5_new-G) + cross(omega_GF_new, cross(omega_GF_new, S5_new-G));

    % Store values for plotting

    new_B_x(theta) = B_new(1);
    new_B_y(theta) = B_new(2);
    new_C_x(theta) = C_new(1);
    new_C_y(theta) = C_new(2);
    new_E_x(theta) = E_new(1);
    new_E_y(theta) = E_new(2);
    new_F_x(theta) = F_new(1);
    new_F_y(theta) = F_new(2);

    % Store angular velocities for plotting
    new_omega_BC(theta) = angularVelocity_BC_new;
    new_omega_CDE(theta) = angularVelocity_CDE_new;
    new_omega_GF(theta) = angularVelocity_GF_new;
    new_omega_FE(theta) = angularVelocity_FE_new;

    % Store linear velocities for plotting
    new_vB(theta) = norm(vB_new);
    new_vC(theta) = norm(vC_new);
    new_vE(theta) = norm(vE_new);
    new_vF(theta) = norm(vF_new);

    % Store angular accelerations for plotting
    new_alpha_BC(theta) = alpha_BC_new;
    new_alpha_CDE(theta) = alpha_CDE_new;
    new_alpha_GF(theta) = alpha_GF_new;
    new_alpha_FE(theta) = alpha_FE_new;

    % Store accelerations of joints for plotting
    new_aB(theta) = norm(a_BA_new);
    new_aC(theta) = norm(a_CA_new);
    new_aE(theta) = norm(a_ED_new);
    new_aF(theta) = norm(a_FG_new);


%Dynamic Forces Equations 
%Link AB
LoopNeqn1 = NF_A_loop + NF_B_loop + WAB == MassAB * a_s1_A_l;
LoopNeqn2 = cross(A-S1_new, NF_A_loop) + cross(B_new - S1_new, NF_B_loop) + NT_in_loop == J_AB * alpha_AB;

% Link BC
LoopNeqn3 = -NF_B_loop + NF_C_loop + WBC == MassBC * a_s2_A_l;
LoopNeqn4 = cross(B_new-S2_new, -NF_B_loop) + cross(C_new-S2_new, NF_C_loop) == J_BC * alphaBC_vector_new;

% Link CDE
LoopNeqn5 = -NF_C_loop + NF_D_loop + NF_E_loop + WCDE == MassCDE * a_s3_D_l;
LoopNeqn6 = cross(C_new-S3_new, -NF_C_loop) + cross(D-S3_new, NF_D_loop) + cross(E_new-S3_new, NF_E_loop) == J_CDE * alphaCDE_vector_new;

% Link EF
LoopNeqn7 = -NF_E_loop + NF_F_loop + WEF + ArtifactWeight == MassEF * a_s4_D_l;
LoopNeqn8 = cross(E_new-S4_new, -NF_E_loop) + cross(F_new-S4_new, NF_F_loop) == J_EF * alphaFE_vector_new;

% Link FG
LoopNeqn9 = -NF_F_loop + NF_G_loop + WFG == MassGF * a_s5_G_l;
LoopNeqn10 = cross(F_new-S5_new, -NF_F_loop) + cross(G-S5_new, NF_G_loop) == J_FG * alphaGF_vector_new;

LoopDynamiceqnMatrix = [LoopNeqn1 LoopNeqn2 LoopNeqn3 LoopNeqn4 LoopNeqn5 LoopNeqn6 LoopNeqn7 LoopNeqn8 LoopNeqn9 LoopNeqn10];
LoopDynamicSolution = solve(LoopDynamiceqnMatrix, [NFAx_l NFAy_l NFBx_l NFBy_l NFCx_l NFCy_l NFDx_l NFDy_l NFEx_l NFEy_l NFFx_l NFFy_l NFGx_l NFGy_l NTin_l]);

% Store dynamic forces for later plotting
NForceAx_l_all(theta) = LoopDynamicSolution.NFAx_l;
NForceAy_l_all(theta) = LoopDynamicSolution.NFAy_l;
NForceBx_l_all(theta) = LoopDynamicSolution.NFBx_l;
NForceBy_l_all(theta) = LoopDynamicSolution.NFBy_l;
NForceCx_l_all(theta) = LoopDynamicSolution.NFCx_l;
NForceCy_l_all(theta) = LoopDynamicSolution.NFCy_l;
NForceDx_l_all(theta) = LoopDynamicSolution.NFDx_l;
NForceDy_l_all(theta) = LoopDynamicSolution.NFDy_l;
NForceEx_l_all(theta) = LoopDynamicSolution.NFEx_l;
NForceEy_l_all(theta) = LoopDynamicSolution.NFEy_l;
NForceFx_l_all(theta) = LoopDynamicSolution.NFFx_l;
NForceFy_l_all(theta) = LoopDynamicSolution.NFFy_l;
NForceGx_l_all(theta) = LoopDynamicSolution.NFGx_l;
NForceGy_l_all(theta) = LoopDynamicSolution.NFGy_l;
NInputTorque(theta) = LoopDynamicSolution.NTin_l;
end

% Joint trajectories graph
figure;            
hold on;
grid on;
axis equal;

plot(new_B_x, new_B_y, 'b-', 'LineWidth', 1.5);
plot(new_C_x, new_C_y, 'r-', 'LineWidth', 1.5);
plot(new_E_x, new_E_y, 'g-', 'LineWidth', 1.5);
plot(new_F_x, new_F_y, 'm-', 'LineWidth', 1.5);
xticks(0:30:360);

xlabel('X Position');
ylabel('Y Position');
title('Joint Trajectories');
legend('B', 'C', 'E', 'F');

hold off;

%Plot Forces (static)

% X Component
forceXData = {ForceAx_l_all, ForceBx_l_all, ForceCx_l_all, ...
              ForceDx_l_all, ForceEx_l_all, ForceFx_l_all, ...
              ForceGx_l_all};

forceXNames = {'ForceAx', 'ForceBx', 'ForceCx', ...
               'ForceDx', 'ForceEx', 'ForceFx', 'ForceGx'};

figure;
hold on;
grid on;

for i = 1:length(forceXData)
    plot(forceXData{i}, 'LineWidth', 1.5);
end

xlabel('Theta (degrees)');
ylabel('Force X (Newtons)');
title('Forces in X-Direction vs Theta');

legend(forceXNames, 'Location', 'best');

hold off;

% Y component

forceYData = {ForceAy_l_all, ForceBy_l_all, ForceCy_l_all, ...
    ForceDy_l_all, ForceEy_l_all, ForceFy_l_all, ...
    ForceGy_l_all};

forceYNames = {'ForceAy', 'ForceBy', 'ForceCy', ...
    'ForceDy', 'ForceEy', 'ForceFy', 'ForceGy'};

figure;
hold on;
grid on;

for i = 1:length(forceYData)
    plot(forceYData{i}, 'LineWidth', 1.5);
end

xlabel('Theta (degrees)');
ylabel('Force Y (Newtons)');
title('Forces in Y-Direction vs Theta');

legend(forceYNames, 'Location', 'best');

hold off;

% Angular Velocities graph
figure;            
hold on;
grid on;
axis equal;

position_of_crank = 1:1:360;

plot(position_of_crank, new_omega_BC, 'b-', 'LineWidth', 1.5);
plot(position_of_crank, new_omega_CDE, 'r-', 'LineWidth', 1.5);
plot(position_of_crank, new_omega_GF, 'g-', 'LineWidth', 1.5);
plot(position_of_crank, new_omega_FE, 'm-', 'LineWidth', 1.5);
xticks(0:30:360);

xlabel('Position of Crank');
ylabel('Angular Velocity of link');
title('Angular Velocity vs Crank Position');
legend('BC', 'CDE', 'GF', 'FE');

hold off;

% Linear Velocities graph
figure;            
hold on;
grid on;
axis equal;

position_of_crank = 1:1:360;

plot(position_of_crank, new_vB, 'b-', 'LineWidth', 1.5);
plot(position_of_crank, new_vC, 'r-', 'LineWidth', 1.5);
plot(position_of_crank, new_vE, 'g-', 'LineWidth', 1.5);
plot(position_of_crank, new_vF, 'm-', 'LineWidth', 1.5);
xticks(0:90:360);

xlabel('Position of crank');
ylabel('Velocity of Joint');
title('Joint Velocities vs Crank Position');
legend('B', 'C', 'E', 'F');

hold off;


% Angular Acceleration graph
figure;            
hold on;
grid on;
axis equal;

position_of_crank = 1:1:360;

plot(position_of_crank, new_alpha_BC, 'b-', 'LineWidth', 1.5);
plot(position_of_crank, new_alpha_CDE, 'r-', 'LineWidth', 1.5);
plot(position_of_crank, new_alpha_GF, 'g-', 'LineWidth', 1.5);
plot(position_of_crank, new_alpha_FE, 'm-', 'LineWidth', 1.5);
xticks(0:90:360);

xlabel('Position of crank');
ylabel('Angular Acceleration of Link');
title('Angular Acceleration of Link vs Crank Position');
legend('BC', 'CDE', 'GF', 'FE');

hold off;

% Linear Acceleration of Joints graph
figure;            
hold on;
grid on;
axis equal;

position_of_crank = 1:1:360;

plot(position_of_crank, new_aB, 'b-', 'LineWidth', 1.5);
plot(position_of_crank, new_aC, 'r-', 'LineWidth', 1.5);
plot(position_of_crank, new_aE, 'g-', 'LineWidth', 1.5);
plot(position_of_crank, new_aF, 'm-', 'LineWidth', 1.5);
xticks(0:90:360);

xlabel('Position of crank');
ylabel('Acceleration of Joint');
title('Acceleration of Joint vs Crank Position');
legend('B', 'C', 'E', 'F');

hold off;

%Plot Forces (Dynamic)

DynamicforceDatax = {NForceAx_l_all, NForceBx_l_all, NForceCx_l_all, NForceDx_l_all, NForceEx_l_all, NForceFx_l_all, NForceGx_l_all, };

DynamicforceNamesx = {'NForceAx', 'NForceBx', 'NForceCx', 'NForceDx', 'NForceEx', 'NForceFx', 'NForceGx'};

figure;
hold on;
grid on;

for i = 1:length(DynamicforceDatax)

    plot(DynamicforceDatax{i}, 'LineWidth', 1.5);

end

xlabel('Theta (degrees)');
ylabel('Force (Newtons)');
title('Dynamic Forces vs Theta');

legend(DynamicforceNamesx, 'Location', 'best');

hold off;

%Plot Forces (Dynamic)

DynamicforceDatay = {NForceAy_l_all,NForceBy_l_all, NForceCy_l_all, NForceDy_l_all, NForceEy_l_all, NForceFy_l_all, NForceGy_l_all};

DynamicforceNamesy = {'NForceAy', 'NForceBy', 'NForceCy', 'NForceDy','NForceEy', 'NForceFy', 'NForceGy'};

figure;
hold on;
grid on;

for i = 1:length(DynamicforceDatay)

    plot(DynamicforceDatay{i}, 'LineWidth', 1.5);

end

xlabel('Theta (degrees)');
ylabel('Force (Newtons)');
title('Dynamic Forces vs Theta');

legend(DynamicforceNamesy, 'Location', 'best');

hold off;
