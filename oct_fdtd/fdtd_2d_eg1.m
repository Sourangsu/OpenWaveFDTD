%% ========================================================================
% Written by: Sourangsu Banerji, University of Utah 
% Verified by: Manjunath Machnoor, University of Southern California
%% ========================================================================

%% code for 2D FDTD (TM mode wave propagation)
%% workspace definition
close all;
clear all;
clc;

%%parameter definition (material - source - structure definition - boundary condition)
IE = 60;                                                                   %number of cells to be used
JE = 60;                                                                   %number of cells to be used

%material definition
epsz = 8.85419e-12;                                                        

%source definition
ic = IE/2;
jc = JE/2;
pi = 3.14159;
to = 20;                                                                   %center of the incident pulse
spread = 6;                                                                %width of the incident pulse
ddx = 0.01;                                                                %spatial sampling
dt = ddx/(2*3e8);                                                          %temporal interval (could be derived from courant stability factor)

T = 0;
Nsteps = 1;

for j = 1:IE
    for i = 1:JE
        dz(i,j) = 0;
        ez(i,j) = 0;
        hx(i,j) = 0;
        hy(i,j) = 0;
        ga(i,j) = 1;
    end
    
end

%% Warning!! Don't change code from here!!
while (Nsteps > 0)
    n = 0;
    
    for n = 1:Nsteps                                                       %Nsteps is the number of times the main loop has executed
        T =T+1;                                                            %T keeps track of the timesteps
        %main fdtd loop
        
        %calculate the Dz field
        for j = 2:IE
            for i = 2:IE
                dz(i,j) = dz(i,j) + 0.5*(hy(i,j)-hy(i-1,j)-hx(i,j)+hx(i,j-1));
            end
        end
            
        %put pulse in the specified grid position
        pulse =  exp(-0.5*((to-T)/spread)^2);
        dz(ic,jc) = pulse;
        
        %calculate Ez from Dz
        for j = 2:JE
            for i = 2:IE
                ez(i,j) = ga(i,j) * dz(i,j);
            end
        end
        
        %calculate the Hx field
        for j = 1:JE-1
            for i = 1:IE-1
                hx(i,j) = hx(i,j) + 0.5*(ez(i,j)-ez(i,j+1));
            end
        end
       
        %calculate the Hy field
        for j = 1:JE-1
            for i = 1:IE-1
                hy(i,j) = hy(i,j) + 0.5*(ez(i+1,j)-ez(i,j));
            end
        end
    end 

pause(0.2);
fprintf('Timestep = %f \n',T);
surf(ez); 

end


