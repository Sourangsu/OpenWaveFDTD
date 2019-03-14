%% ========================================================================
% Written by: Sourangsu Banerji, University of Utah, March 2019 
% Verified by: Manjunath Machnoor, University of Southern California
%% ========================================================================

%% code for 1D FDTD (in free space-PML boundary condition)
%% workspace definition
close all;
clear all;
clc;

%%parameter definition (material - source - structure definition - boundary condition)
MaX = 200;                                                                 %number of cells to be used

%source definition
kc = MaX/2;                                                                %center of the problem space
to = 40;                                                                   %center of the incident pulse
spread = 12;                                                               %width of the incident pulse

%boundary condition
Ex_low_m1 = 0;
Ex_low_m2 = 0;
Ex_high_m1 = 0;
Ex_high_m2 = 0;

T = 0;
Nsteps = 1;

%%field definition
Ex = zeros(1,MaX);                                                         %electric field
Hy = zeros(1,MaX);                                                         %magnetic field

for k = 1:MaX
    Ex(k) = 0;
    Hy(k) = 0;
end


%% Warning!! Don't change code from here!!
while (Nsteps > 0)
    n = 0;
    
    for n = 1:Nsteps                                                       %Nsteps is the number of times the main loop has executed
        T =T+1;                                                            %T keeps track of the timesteps
        %main fdtd loop
        
        %calculate the Ex field
        for k = 2:MaX
            Ex(k) = Ex(k) + 0.5*(Hy(k-1)-Hy(k));
        end
            
        %put gaussian pulse in the middle
        pulse =  exp(-0.5*((to-T)/spread)^2);
        Ex(kc) = Ex(kc)+pulse;
        
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
    plot(Ex);
    xlabel('FDTD cells');
    ylabel('Ex');
    subplot(2,1,2);
    plot(Hy);
    xlabel('FDTD cells');
    ylabel('Hy');
    pause(0.2);
    fprintf('Timestep = %f \n',T);
end