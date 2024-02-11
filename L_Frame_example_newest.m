% Statically indeterminate __| frame
% Assumes that the frame is fixed at both ends
% Jónas Þór Snæbjörnsson, Feb. 2024
%------------------------------------------------------------------------
clear vars
close all

% Structural parameters
Ex = 68.9e9;            % N/m^2 
% Cross-sectional properties of hollow circular tubes (24mm OD, 2mm thick)
outer_diameter = 24e-3;  % Outer diameter (m)
thickness = 2e-3;        % Thickness (m)
inner_diameter = outer_diameter - 2*thickness;  % Inner diameter (m)
Ar = pi*(outer_diameter^2 - inner_diameter^2)/4;  % Cross-sectional area (m^2)
Iz = pi*(outer_diameter^4 - inner_diameter^4)/64;   % Second moment of area (m^4)

% Iz = 280*1e6 * 1e-12;  % m^4
% Ar =  12*1e3 * 1e-6;   % m^2

E = [Ex, Ex, Ex];
A = [Ar, Ar, Ar];
I = [Iz, Iz, Iz];

% Nodal point no., x-coord., y-coord. for all nodes
% N = [1, 0, 0; 2, 3, 0; 3, 6, 0; 4, 6, 4 ; 5.....;6];
N = [1 4, 4 ;2 -1.5, 3.5;3 4.5, 3;4 0,0;5 6, 0.5;6 -4.5,0.5];

% Element no., 1.node no., 2. node no., for all elements
% C = [1 1 2; 2 2 3; 3 3 4; 4 ... ;5 ...; 6....;7 ...   ];
C = [1 1 2; 2 1 3; 3 3 4; 4 2 4; 5 3 5; 6 2 6; 7 4 6];

% Count number of elements, nodal points and degrees of freedom
[nel,~]=size(C);  % nel = number of elements
E = ones(nel,1)*Ex;
A = ones(nel,1)*Ar;
A(6) = 2*Ar;
A(7) = 2*Ar;
I = ones(nel,1)*Iz;
I(6) = 2*Iz;
I(7) = 2*Iz;



[nnp,nc]=size(N);  % nnp = number of nodal points
ndof = 3*nnp;      % 3 dof at each node
% Evaluate element length and directional cosines
L = zeros(1,nel); lx = zeros(1,nel); ly =zeros(1,nel);
for i=1:C(end,1)
    L(i) =sqrt((N(C(i,3),2)-N(C(i,2),2))^2+(N(C(i,3),3)-N(C(i,2),3))^2);
    lx(i)=(N(C(i,3),2)-N(C(i,2),2))/L(i);
    ly(i)=(N(C(i,3),3)-N(C(i,2),3))/L(i);
end

lambda=[lx' ly'];

% Definition of frame (note that element no. 3 is vertical):
% Node no:       1o-------------o2o-----------o3o-----------o4
% Element no:            1               2             3
% DOF no:        1,2,3         4,5,6         7,8,9         10,11,12
% DOF 1, 4, 7 and 10 are horizontal displacements (x-disp.)
% DOF 2, 5, 8 and 11 are vertical displacements (y-disp.)
% DOF 3, 6, 9 and 12 are rotations about the z-axis

% Boundary conditions (bc's) 
% Translations and rotations are fixed at nodes 1 and 4 
% Free DOF are: 4, 5, 6, 7, 8, 9 
% Restrained DOF are: 1,2,3,10,11,12
% Vector defining the bc's: 0 = known dof , 1 = unknown dof')
ir = [1, 1, 1, ... % Node 1 is free
      1, 1, 1, ... % Node 2 is free
      1, 1, 1, ... % Node 3 is free
      1, 1, 1, ... % Node 4 is free
      0, 0, 0, ... % Node 5 is fixed
      0, 0, 0];    % Node 6 is fixed

% Element stiffness matrix
ke = zeros(6,6,nel);      % preallocating matrix size
Te = zeros(6,6,nel);
Ke = zeros(6,6,nel);
for i=1:length(L)
    
    k1 = A(i)*E(i)/L(i);
    k2 = 12*E(i)*I(i)/(L(i)^3);
    k3 = 6*E(i)*I(i)/(L(i)^2);
    k4 = 4*E(i)*I(i)/L(i);
    k5 = 2*E(i)*I(i)/L(i);

    ke(:,:,i) = [ k1  0   0    -k1   0    0
                  0   k2  k3    0   -k2   k3
                  0   k3  k4    0   -k3   k5
                 -k1  0   0     k1   0    0
                  0  -k2 -k3    0    k2  -k3
                  0   k3  k5    0   -k3   k4];
    
    Te(:,:,i) = [ lambda(i,1) lambda(i,2)  0    0            0         0
                 -lambda(i,2) lambda(i,1)  0    0            0         0
                    0            0         1    0            0         0
                    0            0         0  lambda(i,1)  lambda(i,2) 0
                    0            0         0 -lambda(i,2)  lambda(i,1) 0
                    0            0         0    0            0         1 ];
                
    Ke(:,:,i) = Te(:,:,i)'*ke(:,:,i)*Te(:,:,i);                
end

% Expand element stiffness matrices to accunt for all DOF

K = zeros(ndof,ndof,nel);  % Zero expanded element matrices in Global coord.

for ne=1:nel
    i=C(ne,2)*3-2;
    j=C(ne,3)*3-2;
    K(i:i+2,i:i+2,ne)=K(i:i+2,i:i+2,ne)+Ke(1:3,1:3,ne);
    K(j:j+2,j:j+2,ne)=K(j:j+2,j:j+2,ne)+Ke(4:6,4:6,ne);
    K(i:i+2,j:j+2,ne)=K(i:i+2,j:j+2,ne)+Ke(1:3,4:6,ne);
    K(j:j+2,i:i+2,ne)=K(j:j+2,i:i+2,ne)+Ke(4:6,1:3,ne);
end

% Global stiffness matrix evaluated by combining element stifness matrices
KG = sum(K,3);

% Boundary conditions accounted for through the ir-vektor
% Stiffness matrix after implementing boundary values by 
% zero-ing out lines and columns corresponding to known dof

KR=KG;
for i = 1:ndof
    if ir(i) == 0
        KR(i,:) = 0;  % Zero out the i-th row
        KR(:,i) = 0;  % Zero out the i-th column
        KR(i,i) = 1;  % Set the diagonal term to 1
    end
end

% External loading (N & Nm), 
% note that reaction forces at the supports are unknown
% dof: 1 2 3 4  5      6 7 8 9 10 11 12
F = [-80, 0, 0,... 
    -540,0,0,...
    0,0,0,...
    -180,0,0,...
    0,0,0,...
    0,0,0]';
% Displacement vektor
D = KR\F;                % solving for the uknown displacements
displ = D*1000               % mm    

% External Loading and Reaction Forces
Q=KG*D;                      % solving for the unkown reaction forces
Forces = Q/1000              % kN   

%% Plot the frame, undeformed and deformed shape
x=N(:,2); y=N(:,3);
x1=x+D(1:3:end-1)*500; y1=y+D(2:3:end)*1;  % Displacements magnified 1 times

Nn=C(:,2); Nf=C(:,3);

figure; hold
% axis([-1 22 -1 10])
for i = 1:C(end,1)
    plot([x(Nn(i)) x(Nf(i))],[y(Nn(i)) y(Nf(i))],'-o');
    plot([x1(Nn(i)) x1(Nf(i))],[y1(Nn(i)) y1(Nf(i))],'--o');
end

% Calculate Internal Forces: shearforce and moment
Nx = zeros(nel,2); Vy = zeros(nel,2); Mz = zeros(nel,2);
for el=1:nel
    da = D(C(el,2)*3-2:C(el,2)*3);
    db = D(C(el,3)*3-2:C(el,3)*3);
    q=ke(:,:,el)*Te(:,:,el)*[da; db];    % solving for internal moment and shear

    Nx(el,:) = q(1:3:4);        % =[q1(1:3:4); q2(1:3:4); q3(1:3:4)];
    Vy(el,:) = q(2:3:5);        %  [q1(2:3:5); q2(2:3:5); q3(2:3:5)];
    Mz(el,:) = q(3:3:6);        % =[q1(3:3:6); q2(3:3:6); q3(3:3:6)];

    Axial_Shear_Moment(el*2-1:el*2,:) = [Nx(el,1) Vy(el,1) Mz(el,1); Nx(el,1) Vy(el,1) Mz(el,1)]/1000;  % kN
    Element_number(el*2-1:el*2,1) = [C(el,1); C(el,1)];
    Node_number(el*2-1:el*2,1) = [C(el,2); C(el,3)];
end
  
%  Table reporting
    table(Element_number,Node_number,Axial_Shear_Moment)

%  Manually collecting N(x), V(x) and M(x) for plotting internal forces
x1=[0 L(1)];
N1=[-Nx(1,1)  Nx(1,2)]/1000; % kN
V1=[-Vy(1,1)  Vy(1,2)]/1000; % kN
M1=[Mz(1,1)  -Mz(1,2)]/1000; % kNm

x2=[L(1) L(1)+L(2)];
N2=[-Nx(2,1)  Nx(2,2)]/1000;
V2=[-Vy(2,1)  Vy(2,2)]/1000;
M2=[Mz(2,1)  -Mz(2,2)]/1000;

y2=[0 L(3)];
N3=[-Nx(3,1) Nx(3,2)]/1000;
V3=[-Vy(3,1) Vy(3,2)]/1000;
M3=[Mz(3,1) -Mz(3,2)]/1000;

figure;

subplot(3,2,1); plot(x1,N1,x2,N2,[0 6],[0 0]); 
axis([0 6  -10  50]);   
ylabel('Axial (kN)'); xlabel('Length (m)'); pbaspect([2,1,1])
subplot(3,2,3); plot(x1,V1,x2,V2,[x1(2) x2(1)],[V1(2) V2(1)],[0 6],[0 0]); 
axis([0 6 -100 100]);
ylabel('Shear (kN)'); xlabel('Length (m)'); pbaspect([2,1,1])
subplot(3,2,5); plot(x1,M1,x2,M2,[0 6],[0 0]);
axis([0 6 -100 100]);
ylabel('Moment (kNm)'); xlabel('Length (m)'); pbaspect([2,1,1])

subplot(3,2,2); plot(N3,y2,[0 0],[0 4]); axis([-10 50 0 4]) 
xlabel('Axial (kN)'); ylabel('Length (m)'); pbaspect([1,1,1])
subplot(3,2,4); plot(V3,y2,[0 0],[0 4]); axis([-50 10 0 4])
xlabel('Shear (kN)'); ylabel('Length (m)'); pbaspect([1,1,1])
subplot(3,2,6); plot(M3,y2,[0 0],[0 4]); axis([-50  50 0 4]); 
xlabel('Moment (kNm)'); ylabel('Length (m)'); pbaspect([1,1,1])

%%
% -----------------------------------------------------------------------
% Dynamics
%------------------------------------------------------------------------
% Element mass matrix in local coordinates & 
% Transformation matrix (same as for the stiffness matrix) ->
% Element mass matrix in global coordinates

% Density
hro = 2700;   % kg/m^3  

% Define the size of the element matrices
mc = zeros(6,6,nel); ml = zeros(6,6,nel);
Mc = zeros(6,6,nel); Ml = zeros(6,6,nel);

% Build the element matrices
for i=1:length(L)

% Lumped mass matrix
    ml1 = hro*A(i)*L(i)/2;
    ml2 = (hro*A(i)*L(i)/2)*L(i)^2/39;

    ml(:,:,i) = [ ml1   0    0       0       0     0
                  0     ml1  0       0       0     0
                  0        0   ml2    0      0     0
                  0        0    0      ml1   0     0
                  0        0    0      0      ml1  0
                  0        0    0      0        0   ml2];  
% Consistent mass matrix
    m1 = hro*A(i)*L(i)*2/6;
    m2 = hro*A(i)*L(i)*156/420;
    m3 = hro*A(i)*L(i)*22*L(i)/420;
    m4 = hro*A(i)*L(i)*54/420;
    m5 = hro*A(i)*L(i)*13*L(i)/420;
    m6 = hro*A(i)*L(i)*4*L(i)^2/420;
    m7 = hro*A(i)*L(i)*3*L(i)^2/420;
    
    mc(:,:,i) = [ m1   0   0   m1/2    0    0
                  0   m2  m3    0     m4  -m5
                  0   m3  m6    0     m5  -m7
                 m1/2  0   0   m1      0    0
                  0   m4  m5    0     m2  -m3
                  0  -m5 -m7    0    -m3   m6];
    
%    T(:,:,i) = [ lx(i) ly(i)  0   0      0     0
%                -ly(i) lx(i)  0   0      0     0
%                  0     0     1   0      0     0
%                  0     0     0  lx(i)  ly(i)  0
%                  0     0     0 -ly(i)  lx(i)  0
%                  0     0     0   0      0     1 ];

% Use the same transformation as for the stiffness matrix
    Ml(:,:,i) = Te(:,:,i)'*ml(:,:,i)*Te(:,:,i);                
    Mc(:,:,i) = Te(:,:,i)'*mc(:,:,i)*Te(:,:,i);                

end

% Combined Global Mass Matrix in Global coordinates
% Using either Lumped (Ml) or Consistent (M) element mass matrices

% Predefine the size of the global mass matrix
Mcg  = zeros(ndof,ndof,nel);
Mlg = zeros(ndof,ndof,nel);  

for ne=1:nel
    i=C(ne,2)*3-2;
    j=C(ne,3)*3-2;
    
    Mlg(i:i+2,i:i+2,ne)=Mlg(i:i+2,i:i+2,ne)+Ml(1:3,1:3,ne);
    Mlg(j:j+2,j:j+2,ne)=Mlg(j:j+2,j:j+2,ne)+Ml(4:6,4:6,ne);
    Mlg(i:i+2,j:j+2,ne)=Mlg(i:i+2,j:j+2,ne)+Ml(1:3,4:6,ne);
    Mlg(j:j+2,i:i+2,ne)=Mlg(j:j+2,i:i+2,ne)+Ml(4:6,1:3,ne);
    
    Mcg(i:i+2,i:i+2,ne)=Mcg(i:i+2,i:i+2,ne)+Mc(1:3,1:3,ne);
    Mcg(j:j+2,j:j+2,ne)=Mcg(j:j+2,j:j+2,ne)+Mc(4:6,4:6,ne);
    Mcg(i:i+2,j:j+2,ne)=Mcg(i:i+2,j:j+2,ne)+Mc(1:3,4:6,ne);
    Mcg(j:j+2,i:i+2,ne)=Mcg(j:j+2,i:i+2,ne)+Mc(4:6,1:3,ne);
end

% Global stiffness matrix evaluated by combining element stifness matrices
MlG = sum(Mlg,3);
McG = sum(Mcg,3);

% Boundary conditions accounted for through the ir-vektor
% Use reduced matrices to elimate false eigenvalues, 
% i.e. select only active DOF & eliminate rows and columns for fixed DOF

index = find(ir>0);

Mlr = MlG(index,index); 
Mcr = McG(index,index);
Kr  = KG(index,index);

%Solve the eigenvalue problem using lumped and consistent mass matrices
[Msl,Frl] = eig(Kr,Mlr); 
[Msc,Frc] = eig(Kr,Mcr);

[wl,Index] = sort(sqrt(diag(Frl)));   % rad/s
fl = wl/(2*pi);         % 1/s = Hz  
Tl = fl.^-1;            % s
MSl = zeros(ndof,length(index)); MSl(index,:) = Msl(:,Index);

[wc,Index] = sort(sqrt(diag(Frc))); % rad/s
fc = wc/(2*pi);         % 1/s = Hz  
Tc = fc.^-1;            % s  
MSc = zeros(ndof,length(index)); MSc(index,:) = Msc(:,Index);

% Compare eigenfrequencies from Lumped and Consistent mass matrices
mode_no = (1:1:length(wl))';
table(mode_no, wl, wc, fl, fc, Tl, Tc)

% plot mode shapes for the frame, undeformed and deformed shape (ignoring rotation)
% using the mode shapes from the consistent equation of motion
Ms1 = zeros(ndof,length(index));  
Ms1(index,:) = Msc;
Nn=C(:,2); Nf=C(:,3);      % node number near o---- and ----o far

figure
for j = 1:6
    x2 = x+Ms1(1:3:end-2,j)*50; y2 = y+Ms1(2:3:end-1,j)*50;   % Displacements magnified 300 times
    subplot(3,2,j); hold; axis([-1 10 -2 6])
    for i = 1:nel
    plot([x(Nn(i)) x(Nf(i))],[y(Nn(i)) y(Nf(i))],'-o');
    plot([x2(Nn(i)) x2(Nf(i))],[y2(Nn(i)) y2(Nf(i))],'--o');
    legend(['f =',num2str(fc(j))],'Location','Northwest');
    end
end

