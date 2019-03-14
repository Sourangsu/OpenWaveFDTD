%% ========================================================================
% Written by: Sourangsu Banerji, University of Utah
% Verified by: Manjunath Machnoor, University of Southern California
%% ========================================================================

%% code for 1D FDTD (sinusoidal pulse hitting a dielectric medium-pml boundary condition)
%% workspace definition
close all;
clear all;
clc;

%%parameter definition (material - source - structure definition - boundary condition)
MaX = 200;                                                                 %number of cells to be used

%material definition
eps = 4;                                                                   %relative permittivity of material

%source definition
kStart = 100;                                                              %start of the structure
to = 40;                                                                   %center of the incident pulse
spread = 12;                                                               %width of the incident pulse
ddx = 0.01;                                                                %spatial sampling
dt = ddx/(2*3e8);                                                          %temporal interval (could be derived from courant stability factor)
freq_in = 700*1e6;                                                         %frequency of the excitation pulse
pulse_start_grid_point = 5;                                                %grid point of the excitation pulse

%structure definition
def_structure = [zeros(1,(MaX)/2),ones(1,(MaX)/2)];
grid_cells = linspace(1,MaX,MaX);

%boundary condition
Ex_low_m1 = 0;                                                             %leftmost boundary condition
Ex_low_m2 = 0;                                                             %leftmost boundary condition
Ex_high_m1 = 0;                                                            %rightmost boundary condition
Ex_high_m2 = 0;                                                            %rightmost boundary condition

T = 0;
Nsteps = 1;

%%field definition
Ex = zeros(1,MaX);                                                         %electric field
Hy = zeros(1,MaX);                                                         %magnetic field
cB = zeros(1,MaX);                                                         

for k = 1:MaX
    Ex(k) = 0;
    Hy(k) = 0;
end

for k = 1:MaX
    cB(k) = 0.5;
end

for k = kStart:MaX
    cB(k) = 0.5/eps;
end


%% Warning!! Don't change code from here!!
while (Nsteps > 0)
    n = 0;
    
    for n = 1:Nsteps                                                       %Nsteps is the number of times the main loop has executed
        T =T+1;                                                            %T keeps track of the timesteps
        %main fdtd loop
        
        %calculate the Ex field
        for k = 2:MaX
            Ex(k) = Ex(k) + cB(k)*(Hy(k-1)-Hy(k));
        end
            
        %put gaussian pulse in the middle
        pulse =  sin(2*pi*freq_in*dt*T);
        Ex(pulse_start_grid_point) = Ex(pulse_start_grid_point)+pulse;
        
        
        %PML boundary condition
        Ex(1) = Ex_low_m2;
        Ex_low_m2 = Ex_low_m1;
        Ex_low_m1 = Ex(2);
        
        Ex(MaX-1) = Ex_high_m2;
        Ex_high_m2 = Ex_high_m1;
        Ex_high_m1 = Ex(MaX-2);
            
        %calculate the Hy field
        for k = 1:MaX-1
            Hy(k) = Hy(k) + 0.5*(Ex(k)-Ex(k+1));
        end
    end
    
    subplot(2,1,1);
    plot(grid_cells,Ex,grid_cells,def_structure,'--');
    xlabel('FDTD cells');
    ylabel('Ex');
    subplot(2,1,2);
    plot(grid_cells,Hy,grid_cells,def_structure,'--');
    xlabel('FDTD cells');
    ylabel('Hy');
    pause(0.2);
    fprintf('Timestep = %f \n',T);
end