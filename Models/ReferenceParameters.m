%% Reference parameters
CellVolume = 93.78; % micron^3
RefLength = ((3*CellVolume)/(4*pi))^(1/3); % micron
RefShearRate = 100; % 1/s
RefVelocity = RefLength*RefShearRate; % micron/s
RefViscosity = 1.2*10^(-3); % Pa.s = N.s/m^2
RefPressure = RefViscosity*RefShearRate; % Pa = N/m^2