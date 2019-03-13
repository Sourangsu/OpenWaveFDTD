%% code for 1D FDTD (pulse hitting a dielectric medium-pml boundary condition)
%% workspace definition
close all;
clear all;
clc;

MaX = 200;                                                                 %number of cells to be used

Ex_low_m1 = 0;
Ex_low_m2 = 0;
Ex_high_m1 = 0;
Ex_high_m2 = 0;
eps = 4;
kStart = 100;
kc = MaX/2;                                                                %center of the problem space
to = 40;                                                                   %center of the incident pulse
spread = 12;                                                               %width of the incident pulse
T = 0;
Nsteps = 1;
def_structure = [zeros(1,100),ones(1,100)];
grid_cells = linspace(1,200,200);


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
        pulse =  exp(-0.5*((to-T)/spread)^2);
        Ex(5) = Ex(5)+pulse;
        %fprintf('%f %f \n',(to-T),Ex(kc));
        
        %%PML boundary condition
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