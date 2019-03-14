%% ========================================================================
% Written by: Sourangsu Banerji, University of Utah 
% Verified by: Manjunath Machnoor, University of Southern California
%% ========================================================================

%% code for 1D FDTD (gaussian pulse hitting a lossy dielectric medium-pml boundary condition-using electric flux density)
%% Fourier transform to calculate amplitude and phase
%% workspace definition
close all;
clear all;
clc;

%%parameter definition (material - source - structure definition - boundary condition)
MaX = 200;                                                                 %number of cells to be used

%material definition
eps = 4;                                                                   %relative permittivity of material
epsz = 8.85419e-12;                                                        
sig = 0.04;                                                                %conductivity of the material

%source definition
pi = 3.14159;
kStart = 100;                                                              %start of the structure
to = 50;                                                                   %center of the incident pulse
spread = 10;                                                               %width of the incident pulse
ddx = 0.01;                                                                %spatial sampling
dt = ddx/(2*3e8);                                                          %temporal interval (could be derived from courant stability factor)
freq_in = [100,200,500]*1e6;                                               %frequency of the excitation pulse
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
cA = zeros(1,MaX);                                                         

for k = 1:MaX
    Ga(k) = 1;
    Gb(k) = 0;
    Ex(k) = 0;
    Dx(k) = 0;
    Hy(k) = 0;
    Ix(k) = 0;
    mag(k) = 0;
    
    for m = 1:3
        real_pt(m,k) = 0;
        imag_pt(m,k) = 0;
        ampn(m,k) = 0;
        phasen(m,k) = 0;
    end
    
end

%fourier transform of the input pulse
for m = 1:3
    real_in(m) = 0;
    imag_in(m) = 0;
end

for m = 1:3
    arg(m) = 2*pi*dt*freq_in(m);
end

for k = kStart:MaX
    Ga(k) = 1/(eps+((sig*dt)/epsz));
    Gb(k) = (sig*dt)/epsz;
end


%% Warning!! Don't change code from here!!
while (Nsteps > 0)
    n = 0;
    
    for n = 1:Nsteps                                                       %Nsteps is the number of times the main loop has executed
        T =T+1;                                                            %T keeps track of the timesteps
        %main fdtd loop
        
        %calculate the Ex field
        for k = 2:MaX
            Dx(k) = Dx(k) + 0.5*(Hy(k-1)-Hy(k));
        end
            
        %put pulse in the specified grid position
        pulse =  exp(-0.5*((to-T)/spread)^2);
        Dx(pulse_start_grid_point) = Dx(pulse_start_grid_point)+pulse;
        
        %calculate Ex from Dx
        for k = 1:MaX-1
            Ex(k) = Ga(k)*(Dx(k)-Ix(k));
            Ix(k) = Ix(k)+(Gb(k)*Ex(k));
        end
        
        %calculate the fourier transform of Ex
        for k = 1:MaX
            for m = 1:3
                real_pt(m,k) = real_pt(m,k) + cos(arg(m)*T)*Ex(k);
                imag_pt(m,k) = imag_pt(m,k) - sin(arg(m)*T)*Ex(k);
            end
        end
        
        if T < 100
           for m = 1:3
                real_in(m) = real_in(m) + cos(arg(m)*T)*Ex(10);
                imag_in(m) = imag_in(m) - sin(arg(m)*T)*Ex(10);
            end
        end  
        
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
 
for m = 1:3
        amp_in(m) = sqrt(imag_in(m).^2 + real_in(m).^2);
        phase_in(m) = atan2(imag_in(m),real_in(m)); 

        for k = 1:MaX
                ampn(m,k) = (1/amp_in(m))*sqrt(real_pt(m,k).^2 + imag_pt(m,k).^2);
                phasen(m,k) = atan2(imag_pt(m,k),real_pt(m,k))-phase_in(m); 
        end
end

end


