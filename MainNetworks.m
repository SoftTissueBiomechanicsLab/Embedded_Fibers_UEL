%% Main Script for Beam-to-solid coupling modeling
% Based on the work by Steinbrecher et al., 2020
% Please follow tutorial on GitHub for more info. Thanks!
% Sotiris Kakaletsis, 2022-2023
clear; close all; clc
addpath('./Functions')

%% ======================== INPUT =========================================
%------- CONFIGURATION OPTIONS --------------------------------------------
JobNum = 9891; % Abaqus Job Number [0-9999]
JobDescription = 'Network Example'; % Analysis Description comment
MatrixPaths = '/home/datastore/Sotiris/'; % Directory to save matrix *.txt

%------- GEOMETRY OPTIONS -------------------------------------------------
% Solid Domain Dimensions
Lm = 26; % Cube scale for solid domain
W = Lm;  % Global X-direction
D = Lm;  % Global Y-direction
H = Lm;  % Global Z-direction


% Discretization
hsolid = 1/25*Lm; % Solid element length
hbeam = 1.1*hsolid; % Target beam element length
dist_tol = 1e-8; % Distance tolerance for merging beam nodes

% Load Fiber Data
FiberData = load('./FiberDataDemo.mat');

% Undulations Options
Ampl_c = 0.1; % Undulation (fiber crimp) amplitude  
wx = 2; % Periodicity of Undulations
wy = 1;

% Loop individual fibers
exclude_branches = [89];
BeamAxes = {};
for i = 1:size(FiberData.Network.Segment,1)    
    
    r0 = FiberData.Network.Segment(i,1:3)*Lm;
    rf = FiberData.Network.Segment(i,4:6)*Lm;
    lfib(i) = norm(r0-rf); % Nominal fiber length
    
    curve_type='straight';
    if ~ismember(i,exclude_branches)
        curve_type='sin_und';
    end
            
    Ax = Ampl_c*lfib(i);
    Ay = Ampl_c*lfib(i);
        
    BeamAxes = [BeamAxes;{r0, rf, Ax, Ay, wx, wy, [], curve_type, hbeam}];    
end

%------- CONSTITUTIVE MODEL OPTIONS ---------------------------------------
PenaltyConst = 5e-1; % Coupling penalty parameter 
nGP = 6; % Number of Gauss points for beam domain integration

% Beam geometry
BeamOrder = 2; % Linear OR quadratic Timoshenko Beams OR Euler-Bernoulli 
BeamRadius = 0.075;

% Beam material parameters
mod_ratio = 5e3; % Beam-to-solid stiffness ratio
beam_nu = 0.495; % Poisson's Ratio for beams
Ebeam  = 6.5;

% Solid Material Parameters
Matl = 'NHK'; %  NHK (Neo Hookean)
C10 = Ebeam/(4*mod_ratio*(1+beam_nu));
kmurat = 1e3; % Bulk-to-shear modulus ratio (nearly incompressible)
 

%------- BVP AND SOLVER OPTIONS -------------------------------------------
% Displacement Mode
RPL = 'SSyx'; % Displacement mode (SSij:simple shear) OR (UEi: uniax.) for i,j = x,y,z
X_disp = Lm/2; % Prescribed displacement

% Abaqus Solver Options Steps
num_threads = 32; % Number of threads for Abaqus
nSteps = 1;
SolverOptions.InitStep = 1e-2;% Initial step
SolverOptions.minStep = 1e-4; % Min step
SolverOptions.maxStep = 1e-2; % Max step
SolverOptions.n_incr = 50; % Number of increments to ouput results
SolverOptions.STBL = true; % Stabilization boolean
SolverOptions.STBLfac = 2e-4; % If on, stabilization factor

%------- ADVANCED PIPELINE SETTINGS ---------------------------------------
ReDoDiscretization = true; % Boolean to re-calculate coupling matrices
LMOrder = BeamOrder; % Lagrange Multiplier Order, equal to beam order
s_resol = 6; % Segmentation resolution per beam element
COLim = 1e-7; % Cut off limit for  stiffness matrix entries
n_group = 3; % Group DOFs at super element
plot_segmentation  = false; % Boolean to inspect segmentation
nplotbeams = 7; % Beam domain interpolation data points for post-processing
LM_bool = true; % Plot LM field

%--------------------------------------------------------------------------
%---------------------- END OF INPUT --------------------------------------
%--------------------------------------------------------------------------

% EXPORT ALL DATA and call misc. functions
% Make an new working directory
WorkDir = './AbaqusWorkDir/'; % Default, local work directory for Abaqus
AbaqusNewDir = [WorkDir,'Job',num2str(JobNum,'%.4d'),'/'];
if ~exist(AbaqusNewDir, 'dir')
   mkdir(AbaqusNewDir)
end
save([AbaqusNewDir,RPL,'ConfigSettings'])
loadNrun(AbaqusNewDir,RPL)

