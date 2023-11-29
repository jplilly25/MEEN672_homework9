close all; clear all; clc;

% discretization parameters
N = 5; % number of radial lines
M = 5; % number of circumferential lines

% check that N is odd
if mod(N, 2) == 0
    error('N must be odd')
end

% discretization parameters used to split domain
N_half = 0.5 * (N-1);
N_mid = round(0.5 * N);

% count total nodes and total elements
n_nodes = N*M;
E = (N-1)*(M-1);

% domain parameters
R = 0.01; % m
Lx = 0.04; % m
Ly = 0.03; % m

total_arc_length = R * (pi / 2); % m
delta = total_arc_length / N; % local arc lengths

% initialize three groups to put edge nodes in
xc = zeros(1, N);
yc = zeros(1, N);
x_side = zeros(1, N_half);
y_side = zeros(1, N_half);
x_top = zeros(1, N_half);
y_top = zeros(1, N_half);

% assign the corner nodal coordinates
x_mid = Lx;
y_mid = Ly;

% organize edge nodes in three groups (center nodes, side nodes, top nodes)
for k = 1:N
    xc(k) = R * cos(circle_angle_loc(k, N, R));
    yc(k) = R * sin(circle_angle_loc(k, N, R));

    if k < N_mid
        x_side(k) = Lx;
        y_side(k) = (2*Ly)/(N-1) * (k-1);
    elseif k > N_mid && k <= N
        x_top(k-N_mid) = -Lx * ((2/(N-1)*(k-N)));
        y_top(k-N_mid) = Ly;
    end
end

% assemble the outer edge nodes in counter-clockwise order
x_CCW = [x_side, x_mid, x_top];
y_CCW = [y_side, y_mid, y_top];

% generate nodal connectivities matrix
B = zeros(E, 4);
for k = 1:1:N-1
    for l = 1:1:M-1
        e = (k-1) * (M-1) + l; % element counter
        B(e,1) = M * (k-1) + l; 
        B(e,2) = B(e,1) + M; 
        B(e,3) = B(e,1) + 1; 
        B(e,4) = B(e,2) + 1;
    end
end

% organize all nodes in order
xnode = zeros(1, n_nodes);
ynode = zeros(1, n_nodes);
for k = 1:N
    x_line = linspace(xc(k), x_CCW(k), M);
    y_line = linspace(yc(k), y_CCW(k), M);
    for l = 1:M
        ic = (k-1) * M + l;
        xnode(ic) = -x_line(l);
        ynode(ic) = y_line(l);
    end
end

% create the mesh check plot
for e = 1:1:E
    for i = 1:1:4
        % node 1
        if i == 1
            xp(1) = xnode(B(e,1));
            yp(1) = ynode(B(e,1));
            xp(2) = xnode(B(e,2));
            yp(2) = ynode(B(e,2));
        end
        % node 2
        if i == 2
            xp(1) = xnode(B(e,2));
            yp(1) = ynode(B(e,2));
            xp(2) = xnode(B(e,4));
            yp(2) = ynode(B(e,4));
        end
        % node 3
        if i == 3
            xp(1) = xnode(B(e,3));
            yp(1) = ynode(B(e,3));
            xp(2) = xnode(B(e,4));
            yp(2) = ynode(B(e,4));
        end
        % node 4
        if i == 4
            xp(1) = xnode(B(e,3));
            yp(1) = ynode(B(e,3));
            xp(2) = xnode(B(e,1));
            yp(2) = ynode(B(e,1));
        end
        
        % plot the nodes of the e^th element
        plot(xp, yp, 'k-')
        hold on
    end
end

% add plot title and axes labels
title('Mesh Check Plot')
xlabel('X axis (m)');
ylabel('Y axis (m)');
axis([-1.125 * Lx, 0.125 * Lx, -0.125 * Ly, 1.125 * Ly]);

% ax_low = [-1.01, -0.01];
% ax_high = [0.01, 1.01];
% axis([ax_low(q) * Lx, ax_high(q) * Lx, ...
%     ax_low(q) * Ly, ax_high(q) * Ly]);

% initial plotting of radial lines
% for k = 1:N
%     x_line = linspace(xc(k), x_CCW(k), M);
%     y_line = linspace(yc(k), y_CCW(k), M);
%     plot(x_line, y_line, '-o')
%     hold on
% end
