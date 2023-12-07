close all; clear all; clc;

% domain parameters
R = 0.01; % m
Lx = 0.04; % m
Ly = 0.02; % m
kx = 10; % W/(m*C)
ky = 10; % W/(m*C)
betaT = 40; % W/(m^2*C)
T_inf = 20; % C

% discretization parameters
N_radial = 19; % number of radial lines
N_circ = 19; % number of circumferential lines
plot_on = false;

% create the mesh and count nodes and elements
[xnode, ynode, nodes, B, delta] = mesh2D(N_radial, N_circ, R, Lx, Ly, plot_on);
N = size(xnode, 2);
E = size(B, 1);
nodes_per_element = 4;

% store node numbers and values for boundary condition on left boundary
left_nodes = zeros(1,  0.5 * (N_radial-1));
left_temps = zeros(size(left_nodes));
left_temp = 300;
j = 1;
for i = 1:N
    if xnode(i) == -Lx
        left_nodes(j) = nodes(i);
        left_temps(j) = left_temp;
        j = j + 1;
    end
end
nodes_prescribed_bc = left_nodes(:);
prescribed_val = left_temps(:);
Npd = size(left_nodes, 2);

% store node numbers for convection boundary condition on inner circle
circ_nodes = zeros(1, N_radial);
circ_elements = zeros(1, N_radial-1);

j = 1;
for e = 1:E
    if sqrt(xnode(B(e, 1))^2 + ynode(B(e, 1))^2) == R
        circ_nodes(j) = nodes(B(e, 1));
        circ_elements(j) = e;
        j = j + 1;
    end
end

% define the j and i integer arrays 
% define the j array
N_nonprescribed = 0;
for i = 1:N
    flag=0;
    for k=1:Npd 
        if i == nodes_prescribed_bc(k)
            flag = 1;
        end
    end
    if flag == 0 
        N_nonprescribed = N_nonprescribed + 1;
        jarray(N_nonprescribed) = i;
    end
end
 
% Define l array
larray = zeros(1, N);
for i = 1:N_nonprescribed
    istore = jarray(i);
    larray(1, istore) = i;
end

% define a vector with the prescribed value of u at its array location
Pdisp = zeros(1,N);
Pdisp(1, nodes_prescribed_bc) = prescribed_val;
 
% evaluate the derivatives of the shape functions at the four Gauss points
nG = 2;
Gpoint = [-1/sqrt(3), 1/sqrt(3)];
Gwgt = [1 1];
for i = 1:nG
    psi_der_zeta(1,i) = 0.25*(Gpoint(i)-1);
    psi_der_zeta(2,i) = -0.25*(Gpoint(i)-1);
    psi_der_zeta(3,i) = -0.25*(Gpoint(i)+1);
    psi_der_zeta(4,i) = 0.25*(Gpoint(i)+1);
    psi_der_eta(1,i) = 0.25*(Gpoint(i)-1);
    psi_der_eta(2,i) = -0.25*(Gpoint(i)+1);
    psi_der_eta(3,i) = -0.25*(Gpoint(i)-1);
    psi_der_eta(4,i) = 0.25*(Gpoint(i)+1);
end

% Initialize Kc and fc to zero
Kc = zeros(N_nonprescribed, N_nonprescribed);
fc = zeros(N_nonprescribed, 1);

% loop through all elements and update Kc and fc
for e = 1:E
    % form the nodal coordinate array for this element 
    xn(1) = xnode(B(e,1));
    xn(2) = xnode(B(e,2));
    xn(3) = xnode(B(e,3));
    xn(4) = xnode(B(e,4)); 
    yn(1) = ynode(B(e,1));
    yn(2) = ynode(B(e,2));
    yn(3) = ynode(B(e,3));
    yn(4) = ynode(B(e,4));
    XYmat = [xn(1) yn(1); xn(2) yn(2); xn(3) yn(3); xn(4) yn(4)];

    % form element stiffness matrix for element e by ISO-P method
    Ke = zeros(4, 4);
    for r=1:4
        for s=1:4
            
            % perform Gauss numerical integration 
            for i = 1:nG
                for j = 1:nG
                    psider=[psi_der_zeta(1,j) psi_der_zeta(2,j) ...
                            psi_der_zeta(3,j) psi_der_zeta(4,j); 
                            psi_der_eta(1,i)  psi_der_eta(2,i)  ...
                            psi_der_eta(3,i)  psi_der_eta(4,i)];
                    jacobian_transform = psider * XYmat;        % Jacobian
                    detJmat = det(jacobian_transform);          % determinant of Jacobian
                    jacobian_inverse = inv(jacobian_transform); % inverse of Jacobian
                    Frs = detJmat*(... 
                            (jacobian_inverse(1,1)*psi_der_zeta(r,j) + ...
                            jacobian_inverse(1,2)*psi_der_eta(r,i))*kx* ...
                            (jacobian_inverse(1,1)*psi_der_zeta(s,j) + ...
                            jacobian_inverse(1,2)*psi_der_eta(s,i)) + ...
                            (jacobian_inverse(2,1)*psi_der_zeta(r,j) + ...
                            jacobian_inverse(2,2)*psi_der_eta(r,i))*ky* ...
                            (jacobian_inverse(2,1)*psi_der_zeta(s,j) + ...
                            jacobian_inverse(2,2)*psi_der_eta(s,i)));
                    Ke(r,s) = Ke(r,s) + Frs*Gwgt(i)*Gwgt(j);
                end
            end
        end
    end

    % form element H matrix and element P vector by ISO-P method for
    %   convective boundary condition on edge 1-2
    He = zeros(4, 4);
    Pe = zeros(4, 1);
    
    for i = 1:size(circ_elements, 2)
        if e == circ_elements(i)
            He(1,1) = 2;
            He(1,2) = 1;
            He(2,1) = 1;
            He(2,2) = 2;

            Pe(1,1) = 1;
            Pe(2,1) = 1;
        end
    end
    h12 = delta;
    He = h12 * betaT / 6 * He;
    Pe = h12 * betaT * T_inf / 2 * Pe;
    
    if e == 1
        print_output('K1', Ke)
        print_output('H1', He)
        print_output('P1', Pe)
    end

    % assemble Kc 
    for r = 1:nodes_per_element
        for s = 1:nodes_per_element
            lBer = larray(B(e,r));
            lBes = larray(B(e,s));
            if lBer ~= 0
                if lBes ~= 0
                    Kc(lBer, lBes) = Kc(lBer, lBes) + Ke(r,s) + He(r,s);
                end
            end
             
            % assemble load contributions due to prescribed u's
            if lBer ~= 0
                if lBes == 0
                    fc(lBer,1) = fc(lBer,1) - Pdisp(B(e,s))*Ke(r,s) + Pe(r,1); 
                end
            end
         
        end
    end
end

% solve for the condense system u's
Uc = Kc \ fc;

% assemble all u's
% add in prescribed values first
for i = 1:Npd
    Usystem(1, nodes_prescribed_bc(i)) = prescribed_val(i); 
end

% add in calculated values second
for i = 1:N_nonprescribed
    k = jarray(i);
    Usystem(1, k) = Uc(i);
end

% make a contour plot of the temperature field
figure

XS = zeros(N_radial, N_circ);
YS = zeros(N_radial, N_circ);

y_zero = zeros(1, N_circ);

for i = 1:N_circ
    for j = 1:N_radial
        ic = (i-1) * N_circ + j;
        XS(i, j) = xnode(ic);
        YS(i, j) = ynode(ic);
        ZS(i, j) = Usystem(N_circ*(i-1) + j);

        if YS(i,j) == 0
            ic;
            y_zero(j) = ZS(i,j);
        end
    end
end

print_output('y_zero', y_zero')

contourf(XS, YS, ZS)
c = colorbar;
c.Label.String = 'Temperature (C)';
colormap ((jet));
title('T Contour Plot ISO-P')
xlabel('x');
ylabel('y');
