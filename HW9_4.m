close all; clear all; clc;

% domain parameters
R = 0.01; % m
Lx = 0.04; % m
Ly = 0.03; % m

% discretization parameters
N_radial = 5; % number of radial lines
N_circ = 5; % number of circumferential lines
plot_on = true;

% create the mesh and count nodes and elements
[xnode, ynode, B] = mesh2D(N_radial, N_circ, R, Lx, Ly, plot_on);
N = size(xnode, 2);
E = size(B, 1);



