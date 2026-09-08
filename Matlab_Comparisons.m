clc
clear

theta = linspace(1, 360, 360);

PMKSData = importdata('Kinematics Data.xlsx');

clc
clear
theta = linspace(1, 360, 360);
PMKSData = importdata('Kinematics Data.xlsx');
 
%% Joint
% Joint A
PositionA     = [PMKSData.data(:, 1),  PMKSData.data(:, 2)];
VelocityA     = [PMKSData.data(:, 3),  PMKSData.data(:, 4)];
AccelerationA = [PMKSData.data(:, 5),  PMKSData.data(:, 6)];
 
% Joint B
PositionB     = [PMKSData.data(:, 7),  PMKSData.data(:, 8)];
VelocityB     = [PMKSData.data(:, 9),  PMKSData.data(:, 10)];
AccelerationB = [PMKSData.data(:, 11), PMKSData.data(:, 12)];
 
% Joint D
PositionD     = [PMKSData.data(:, 13), PMKSData.data(:, 14)];
VelocityD     = [PMKSData.data(:, 15), PMKSData.data(:, 16)];
AccelerationD = [PMKSData.data(:, 17), PMKSData.data(:, 18)];
 
% Joint E
PositionE     = [PMKSData.data(:, 19), PMKSData.data(:, 20)];
VelocityE     = [PMKSData.data(:, 21), PMKSData.data(:, 22)];
AccelerationE = [PMKSData.data(:, 23), PMKSData.data(:, 24)];
 
% Joint C
PositionC     = [PMKSData.data(:, 25), PMKSData.data(:, 26)];
VelocityC     = [PMKSData.data(:, 27), PMKSData.data(:, 28)];
AccelerationC = [PMKSData.data(:, 29), PMKSData.data(:, 30)];
 
% Joint G
PositionG     = [PMKSData.data(:, 31), PMKSData.data(:, 32)];
VelocityG     = [PMKSData.data(:, 33), PMKSData.data(:, 34)];
AccelerationG = [PMKSData.data(:, 35), PMKSData.data(:, 36)];
 
% Joint P
PositionP     = [PMKSData.data(:, 37), PMKSData.data(:, 38)];
VelocityP     = [PMKSData.data(:, 39), PMKSData.data(:, 40)];
AccelerationP = [PMKSData.data(:, 41), PMKSData.data(:, 42)];
 
% Joint F
PositionF     = [PMKSData.data(:, 43), PMKSData.data(:, 44)];
VelocityF     = [PMKSData.data(:, 45), PMKSData.data(:, 46)];
AccelerationF = [PMKSData.data(:, 47), PMKSData.data(:, 48)];
 
%% Links (angular)
% Link AB
AngPos_AB = PMKSData.data(:, 49);   % deg
AngVel_AB = PMKSData.data(:, 50);   % rad/s
AngAcc_AB = PMKSData.data(:, 51);   % rad/s^2
 
% Link CB
AngPos_CB = PMKSData.data(:, 52);
AngVel_CB = PMKSData.data(:, 53);
AngAcc_CB = PMKSData.data(:, 54);
 
% Link DEC
AngPos_DEC = PMKSData.data(:, 55);
AngVel_DEC = PMKSData.data(:, 56);
AngAcc_DEC = PMKSData.data(:, 57);
 
% Link FE
AngPos_FE = PMKSData.data(:, 58);
AngVel_FE = PMKSData.data(:, 59);
AngAcc_FE = PMKSData.data(:, 60);

% Link GPF
AngPos_GPF = PMKSData.data(:, 61);
AngVel_GPF = PMKSData.data(:, 62);
AngAcc_GPF = PMKSData.data(:, 63);

MatlabData = load('MechanismResults.mat');


%% Comparisons

% Position differences [x, y]
N = min(size(PositionB,1), size(MatlabData.DatasetArrayPosition,1));

deltaPositionB = PositionB(1:N,:) - MatlabData.DatasetArrayPosition(1:N,1:2);
deltaPositionC = PositionC(1:N,:) - MatlabData.DatasetArrayPosition(1:N,3:4);
deltaPositionE = PositionE(1:N,:) - MatlabData.DatasetArrayPosition(1:N,5:6);
deltaPositionF = PositionF(1:N,:) - MatlabData.DatasetArrayPosition(1:N,7:8);
% Velocity
deltaVelocityB = vecnorm(VelocityB(1:N,:),2,2) - MatlabData.DatasetArrayLinVel(1:N,1);
deltaVelocityC = vecnorm(VelocityC(1:N,:),2,2) - MatlabData.DatasetArrayLinVel(1:N,2);
deltaVelocityE = vecnorm(VelocityE(1:N,:),2,2) - MatlabData.DatasetArrayLinVel(1:N,3);
deltaVelocityF = vecnorm(VelocityF(1:N,:),2,2) - MatlabData.DatasetArrayLinVel(1:N,4);

% Acceleration
deltaAccelerationB = vecnorm(AccelerationB(1:N,:),2,2) - MatlabData.DatasetArrayLinAcc(1:N,1);
deltaAccelerationC = vecnorm(AccelerationC(1:N,:),2,2) - MatlabData.DatasetArrayLinAcc(1:N,2);
deltaAccelerationE = vecnorm(AccelerationE(1:N,:),2,2) - MatlabData.DatasetArrayLinAcc(1:N,3);
deltaAccelerationF = vecnorm(AccelerationF(1:N,:),2,2) - MatlabData.DatasetArrayLinAcc(1:N,4);

% Angular Velocity
deltaAngVel_CB  = AngVel_CB(1:N)  - MatlabData.DatasetArrayAngVel(1:N,1);
deltaAngVel_DEC = AngVel_DEC(1:N) - MatlabData.DatasetArrayAngVel(1:N,2);
deltaAngVel_GPF = AngVel_GPF(1:N) - MatlabData.DatasetArrayAngVel(1:N,3);
deltaAngVel_FE  = AngVel_FE(1:N)  - MatlabData.DatasetArrayAngVel(1:N,4);

% Angular Acceleration
deltaAngAcc_CB  = AngAcc_CB(1:N)  - MatlabData.DatasetArrayAngAcc(1:N,1);
deltaAngAcc_DEC = AngAcc_DEC(1:N) - MatlabData.DatasetArrayAngAcc(1:N,2);
deltaAngAcc_GPF = AngAcc_GPF(1:N) - MatlabData.DatasetArrayAngAcc(1:N,3);
deltaAngAcc_FE  = AngAcc_FE(1:N)  - MatlabData.DatasetArrayAngAcc(1:N,4);

% Averages

avgDeltaVelocity = mean([deltaVelocityB; deltaVelocityC; deltaVelocityE; deltaVelocityF], 'omitnan')

avgDeltaAcceleration = mean([deltaAccelerationB; deltaAccelerationC; deltaAccelerationE; deltaAccelerationF], 'omitnan')

avgDeltaAngVel = mean([deltaAngVel_CB; deltaAngVel_DEC; deltaAngVel_GPF; deltaAngVel_FE], 'omitnan')

avgDeltaAngAcc = mean([deltaAngAcc_CB; deltaAngAcc_DEC; deltaAngAcc_GPF; deltaAngAcc_FE], 'omitnan')

