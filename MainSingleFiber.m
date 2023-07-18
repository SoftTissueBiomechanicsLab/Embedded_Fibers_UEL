%% Main Script for Beam-to-solid coupling modeling
% Based on the work by Steinbrecher et al., 2020
% Please follow tutorial on GitHub for more info. Thanks!
% Sotiris Kakaletsis, 2022-2023
clear; close all; clc
addpath('./Functions')

%% ======================== INPUT =========================================
%------- CONFIGURATION OPTIONS --------------------------------------------
JobNum = 9888; % Abaqus Job Number [0-9999]
JobDescription = 'Helix Coil Example'; % Analysis Description comment
MatrixPaths = '/home/datastore/Sotiris/'; % Directory to save matrix *.txt

%------- GEOMETRY OPTIONS -------------------------------------------------
% Solid Domain Dimensions
W = 1;  % Global X-direction
D = 2;  % Global Y-direction
H = 1;  % Global Z-direction

% Beam/Fiber Geometry
r0 = [W/2, 0.05*D, H/2]; % Start point
rf = [W/2, 0.95*D, H/2]; % End point

% Discretization
hsolid = 1/17; % Solid element length
hbeam = 1.5*hsolid; % Beam element length
dist_tol = 1e-8; % Distance tolerance for merging beam nodes

% Introduce Undulations
curve_type='helix'; % 'straight' OR 'helix' OR 'sin_und'
Ax = min([W,H])*0.4; Ay=[]; % Undulations Amplitude
wx=[]; wy=[];  % Unudations frequency (applicable only for 'sin_und'
nloops=3; % Complete loops (only for 'helix')
BeamAxes = {r0, rf, Ax, Ay, wx, wy, nloops, curve_type, hbeam}; % Gather

%------- CONSTITUTIVE MODEL OPTIONS ---------------------------------------
PenaltyConst = 1e-1; % Coupling penalty parameter 
nGP = 6; % Number of Gauss points for beam domain integration

% Solid Material Parameters
Matl = 'NHK'; %  NHK (Neo Hookean)
C10 = 1e-3;
kmurat = 1e3; % Bulk-to-shear modulus ratio (nearly incompressible)
 
% Beam geometry
BeamOrder = 2; % Linear OR quadratic Timoshenko Beams OR Euler-Bernoulli 
BeamRadius = 0.017;

% Beam material parameters
mod_ratio = 1.5e4; % Beam-to-solid stiffness ratio
beam_nu = 0.495; % Poisson's Ratio for beams
Ebeam  = 4*C10*mod_ratio*(1+beam_nu);



%------- BVP AND SOLVER OPTIONS -------------------------------------------
% Displacement Mode
RPL = 'UEy'; % Displacement mode (SSij:simple shear) OR (UEi: uniax.) for i,j = x,y,z
X_disp = 2; % Prescribed displacement

% Abaqus Solver Options Steps
num_threads = 32; % Number of threads for Abaqus
nSteps = 1;
SolverOptions.InitStep = 1e-2;% Initial step
SolverOptions.minStep = 1e-4; % Min step
SolverOptions.maxStep = 5e-2; % Max step
SolverOptions.n_incr = 50; % Number of increments to ouput results
SolverOptions.STBL = true; % Stabilization boolean
SolverOptions.STBLfac = 2e-4; % If on, stabilization factor

%------- ADVANCED PIPELINE SETTINGS ---------------------------------------
ReDoDiscretization = true; % Boolean to re-calculate coupling matrices
LMOrder = BeamOrder; % Lagrange Multiplier Order, equal to beam order
s_resol = 17; % Segmentation resolution per beam element
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

