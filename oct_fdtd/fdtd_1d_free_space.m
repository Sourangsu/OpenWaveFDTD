%% code for 1D FDTD (in free space)
%% workspace definition
close all;
clear all;
clc;

MaX = 200;

%%field definition
Ex = zeros(1,MaX);                                                         %electric field
Hy = zeros(1,MaX);                                                         %magnetic field


for k = 1:MaX
    Ex(k) = 0;
    Hy(k) = 0;
end

kc = MaX/2;                                                                %center of the problem space
to = 40;                                                                   %center of the incident pulse
spread = 12;                                                               %width of the incident pulse
T = 0;
Nsteps = 1;

while (Nsteps > 0)
    %fprintf('Nsteps -----> %d \n',Nsteps);
    n = 0;
    
    for n = 1:Nsteps
        T =T+1;
        
        %main fdtd loop
        
        %calculate the Ex field
        for k = 2:MaX
            Ex(k) = Ex(k) + 0.5*(Hy(k-1)-Hy(k));
        end
            
        %put gaussian pulse in the middle
        pulse =  exp(-0.5*((to-T)/spread)^2);
        Ex(kc) = pulse;
        %fprintf('%f %f \n',(to-T),Ex(kc));
            
        %calculate the Hy field
        for k = 1:MaX-1
            Hy(k) = Hy(k) + 0.5*(Ex(k)-Ex(k+1));
        end
    end
    
    plot(Ex);
    pause(0.5);
    fprintf('Timestep = %f \n',T);
end

