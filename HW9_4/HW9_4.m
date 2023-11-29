close all; clear all; clc;

% domain parameters
R = 0.01; % m
Lx = 0.04; % m
Ly = 0.03; % m
kx = 10; % W/(m*C)
ky = 10; % W/(m*C)
betaT = 40; % W/(m^2*C)
T_inf = 20; % C

% discretization parameters
N_radial = 25; % number of radial lines
N_circ = 25; % number of circumferential lines
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
% j = 1;
% for n = 1:N
%     if sqrt(xnode(n)^2 + ynode(n)^2) == R
%         circ_nodes(j) = nodes(n);
%         j = j + 1;
%     end
% end

j = 1;
for e = 1:E
    if sqrt(xnode(B(e, 1))^2 + ynode(B(e, 1))^2) == R
        circ_nodes(j) = nodes(B(e, 1));
        circ_elements(j) = e;
        j = j + 1;
    end
end

% DEFINE INTEGER ARRAYS 
% Define j array from the nodes with prescribed u's
Nnpd=0 ; % counter for number of non-prescribed nodes
for i = 1:1:N
 flag=0;
 for k=1:1:Npd 
     if i == nodes_prescribed_bc(k)
        flag = 1;
     end
 end
     if flag == 0 
        Nnpd = Nnpd + 1;
     jarray(Nnpd) = i;
     end
end

 % Note at this point Nnpd is the total number of non-prescribed u's .
 % These u's must be solved for.
 
% Define L array
larray=zeros(1,N);
for i=1:1:Nnpd
    istore = jarray(i);
    larray(1,istore)=i;
end

% Define a Nx1 vector that has the prescribed value of u at its respective
% array location and zeros at all other array locations 
Pdisp = zeros(1,N);
Pdisp(1, nodes_prescribed_bc) = prescribed_val;
% for i = 1:1:Npd
%     Pdisp(1, nodes_prescribed_bc(i))=prescribed_val(i);
% end
 
 
% Evaluate the derivatives of the shape functions at the Gauss Points 
% These values do not change between elements, so do only once for all elements
nG = 2;
Gpoint = [-1/sqrt(3), 1/sqrt(3)];
Gwgt = [1 1];
for i=1:1:nG
    psi_der_zeta(1,i) = 0.25*(Gpoint(i)-1) ;
    psi_der_zeta(2,i) = -0.25*(Gpoint(i)-1) ;
    psi_der_zeta(3,i) = -0.25*(Gpoint(i)+1) ;
    psi_der_zeta(4,i) = 0.25*(Gpoint(i)+1) ;
    psi_der_eta(1,i) = 0.25*(Gpoint(i)-1) ;
    psi_der_eta(2,i) = -0.25*(Gpoint(i)+1) ;
    psi_der_eta(3,i) = -0.25*(Gpoint(i)-1) ;
    psi_der_eta(4,i) = 0.25*(Gpoint(i)+1) ;
end


% ASSEMBLE THE CONDENSED FORM OF THE SYSTEM (GLOBAL) STIFFNESS
% MATRIX AND FORCE VECTOR
% Initialize Kc and fc to zero
Kc = zeros(Nnpd, Nnpd);
fc = zeros(Nnpd, 1);
 
for e = 1:E % Element Loop 
    % Form the nodal coordinate array for this element 
    xn(1) = xnode(B(e,1));
    xn(2) = xnode(B(e,2));
    xn(3) = xnode(B(e,3));
    xn(4) = xnode(B(e,4)); 
    yn(1) = ynode(B(e,1));
    yn(2) = ynode(B(e,2));
    yn(3) = ynode(B(e,3));
    yn(4) = ynode(B(e,4));
    XYmat=[xn(1) yn(1); xn(2) yn(2); xn(3) yn(3); xn(4) yn(4)];

    % form element stiffness matrix for element e by ISO-P method
    Ke = zeros(4, 4);
    for r=1:1:4
        for s=1:1:4
            Ke(r, s) = 0;
            
            % Perform Gauss Numerical Integration 
            for i=1:1:nG % zeta Gauss Point Loop
                for j=1:1:nG % eta Gauss Point Loop
                    psider=[psi_der_zeta(1,j) psi_der_zeta(2,j) psi_der_zeta(3,j) psi_der_zeta(4,j); 
                            psi_der_eta(1,i)  psi_der_eta(2,i)  psi_der_eta(3,i)  psi_der_eta(4,i)];
                    jacobian_transform = psider * XYmat; % Jacobian
                    detJmat = det(jacobian_transform); % determinant of Jacobian
                    jacobian_inverse = inv(jacobian_transform) ; % inverse of Jacobian
                    Frs = detJmat*(... 
                            (jacobian_inverse(1,1)*psi_der_zeta(r,j) + jacobian_inverse(1,2)*psi_der_eta(r,i))*kx* ...
                            (jacobian_inverse(1,1)*psi_der_zeta(s,j) + jacobian_inverse(1,2)*psi_der_eta(s,i)) + ...
                            (jacobian_inverse(2,1)*psi_der_zeta(r,j) + jacobian_inverse(2,2)*psi_der_eta(r,i))*ky* ...
                            (jacobian_inverse(2,1)*psi_der_zeta(s,j) + jacobian_inverse(2,2)*psi_der_eta(s,i)));
                    Ke(r,s) = Ke(r,s) + Frs*Gwgt(i)*Gwgt(j);
                end % for j=1:1:nG
            end % for i=1:1:nG
         
        end % for s=1:1:4
    end % for i=1:1:nG

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
    
    % Assemble the condensed stiffness matrix Kc 
    for r = 1:nodes_per_element
        for s = 1:nodes_per_element
            lBer=larray( B(e,r) );
            lBes=larray( B(e,s) );
            if lBer ~= 0
                if lBes ~= 0
                    Kc(lBer,lBes) = Kc(lBer,lBes) + Ke(r,s) + He(r,s); % assemble Kc
                end
            end
             
             
            if lBer~=0
                if lBes ==0
                    % assemble load contributions due to prescribed u's
                    fc(lBer,1) = fc(lBer,1) - Pdisp(B(e,s))*Ke(r,s) + Pe(r,1); 
                end
            end
         
        end % s loop
    end % r loop
    Kestore(1:4,1:4,e) = Ke;
 end % e (element) loop

% SOLVE FOR THE CONDENSED SYSTEM u's (i.e. the non prescribed values of u)
Uc = Kc \ fc;
% ASSEMBLE AN Nx1 VECTOR OF ALL u's (prescribed and non-prescribed)
for i=1:1:Npd % first, enter the prescribed values 
 Usystem(1,nodes_prescribed_bc(i))= prescribed_val(i); 
end
for i=1:1:Nnpd % second, enter the calculated (non-prescribed) values 
 k=jarray(i);
 Usystem(1,k)=Uc(i);
end
Usystem ; % display the Uc and Usystem vectors



% MAKE A SURFACE PLOT OF T/To
figure
[XS, YS] = meshgrid(xnode, ynode);

XS = reshape(xnode, [N_radial N_circ]);
YS = reshape(ynode, [N_radial N_circ]);

for i = 1:N_circ
    for j = 1:N_radial
        ic = (i-1) * N_circ + j;
        ZS(i, j) = Usystem(N_circ*(i-1) + j);
    end
end



% subplot(2,1,2)
contourf(XS, YS, ZS)
xlim([-0.045, 0.005])
ylim([-0.005 0.035])

% surf(XS,YS,ZS)
title('T/To Surface Plot ISO-P')
% xlabel('x/a');
% ylabel('y/a');
% zlabel('T/To');
% view(-10,50);