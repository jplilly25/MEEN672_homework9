close all; clear all; clc;

% INPUT NUMERICAL VALUES
Cline = 7 ;% no. of circumferential lines
Rline = 11 ;% no. of radial lines
edof = 4; % number of nodes per element
nG = 2 ; % Order of Guass Quadrature
Gpoint = [-1/sqrt(3) 1/sqrt(3)]; % Gauss Points
Gwgt = [1 1] ; % Gauss Weghts

R = 0.01;
Lx = 0.04;
Ly = 0.02;
delta = (pi/2) / (Cline-1);

% Generate nodes in tha problem domain
N = Cline * Rline; % No. of nodes in the model
xnode = zeros(1, N);
ynode = zeros(1, N);
rnode = zeros(1, N);
thetanode = zeros(1, N);


for k = 1:1:Cline
    for l = 1:1:0.5*(Rline-1)
        ic = (k-1) * Rline + l; % node counter
%         xnode(ic) = 3/(Rline-1) * (k-1);
%         ynode(ic) = 2/(Cline-1) * (l-1);
        rnode(ic) = R + ((Lx/cos((l-1)*delta))-R)/(Cline-1) * (k-1);
        thetanode(ic) = (l-1) * delta;
    end
end

xnode = rnode .* cos(thetanode);
ynode = rnode .* sin(thetanode);
plot(xnode, ynode)
figure
% Generate element nodal connectivities in tha problem domain
for k = 1:1:Cline-1
    for l = 1:1:Rline-1
        e = (k-1) * (Rline-1) + l; % element counter
        B(e,1) = Rline * (k-1) + l; 
        B(e,2) = B(e,1) + Rline; 
        B(e,3) = B(e,1) + 1; 
        B(e,4) = B(e,2) + 1;
    end
end

% Store Node Numbers and Values for Boundary Conditions
% for i = 1:1:Rline
%     dp(i) = i * Cline; % Nodes on top boundary
%     Uprescribed(i) = cos(pi*(i-1)*3/(Rline-1)/6); % Prescribed T's on top boundary
% end
% 
% for i = 1:1:Cline-1
%     dp(Rline+i) = (Rline-1) * Cline + i; % Nodes on right side boundary
%     Uprescribed(Rline+i) = 0 ; % Prescribed T's on right side boundary
% end
E = size(B,1); % number of elements in the model
% Npd = size(dp,2) ; % number of nodes with prescribed u's


for e = 1:1:E
    for i = 1:1:4
        if i == 1
            xp(1) = xnode(B(e,1));
            yp(1) = ynode(B(e,1));
            xp(2) = xnode(B(e,2));
            yp(2) = ynode(B(e,2));
        end
        if i == 2
            xp(1) = xnode(B(e,2));
            yp(1) = ynode(B(e,2));
            xp(2) = xnode(B(e,4));
            yp(2) = ynode(B(e,4));
        end
        if i == 3
            xp(1) = xnode(B(e,3));
            yp(1) = ynode(B(e,3));
            xp(2) = xnode(B(e,4));
            yp(2) = ynode(B(e,4));
        end
        if i == 4
            xp(1) = xnode(B(e,3));
            yp(1) = ynode(B(e,3));
            xp(2) = xnode(B(e,1));
            yp(2) = ynode(B(e,1));
        end

        if e == 1
            xp, yp
        end
        plot(xp, yp, 'k-')
        hold on
    end % end i loop
    
    if e == 1
        title('Mesh Check Plot')
        xlabel(' X axis - units=a ');
        ylabel(' Y axis - units=a ');
        axis([-0.01 * Lx, 1.01 * Lx, -0.01 * Ly, 1.01 * Ly]);
        hold on
    end
end % end e loop