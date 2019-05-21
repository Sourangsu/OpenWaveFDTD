%% ========================================================================
% Written by: Sourangsu Banerji, University of Utah
% Verified by: Manjunath Machnoor, University of Southern California
%% ========================================================================

%% code for 2D FDTD (TM mode plane wave propagation-pml boundary condition)
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
ia = 7;                                                                    %total field/scattered field
ib = IE-ia-1;                                                              %total field/scattered field
ja = 7;                                                                    %total field/scattered field
jb = JE-ja-1;                                                              %total field/scattered field
pi = 3.14159;
to = 20;                                                                   %center of the incident pulse
spread = 8;                                                                %width of the incident pulse
ddx = 0.01;                                                                %spatial sampling
dt = ddx/(2*3e8);                                                          %temporal interval (could be derived from courant stability factor)

T = 0;
Nsteps = 1;

%boundary condition
ez_inc_low_m1 = 0;
ez_inc_low_m2 = 0;
ez_inc_high_m1 = 0;
ez_inc_high_m2 = 0;

ez_inc = zeros(JE,1);
hx_inc = zeros(JE,1);


for j = 1:JE
    for i = 1:IE
        dz(i,j) = 0;
        hx(i,j) = 0;
        hy(i,j) = 0;
        ihx(i,j) = 0;
        ihy(i,j) = 0;
        ga(i,j) = 1;
    end
end

for i = 1:IE
    gi2(i) = 1;
    gi3(i) = 1;
    fi1(i) = 0;
    fi2(i) = 1;
    fi3(i) = 1;
end

for j = 1:IE
    gj2(j) = 1;
    gj3(j) = 1;
    fj1(j) = 0;
    fj2(j) = 1;
    fj3(j) = 1;
end

prompt = 'What is the number of PML cells you wish to use? ';
npml = input(prompt);

for i = 1:npml
    xnum = npml - i;
    xd = npml;
    xxn = xnum/xd;
    xn = (0.33)*power(xxn,3);
    gi2(i) = 1/(1+xn);
    gi2(IE-1-i) = 1/(1+xn);
    gi3(i) = (1-xn)/(1+xn);
    gi3(IE-i-1) = (1-xn)/(1+xn);
    xxn = (xnum-0.5)/xd;
    xn = (0.25)*power(xxn,3);
    fi1(i) = xn;
    fi1(IE-2-i) = xn;
    fi2(i) = 1/(1+xn);
    fi2(IE-2-i) = 1/(1+xn);
    fi3(i) = (1-xn)/(1+xn);
    fi3(IE-2-i) = (1-xn)/(1+xn);
end

for j = 1:npml
    xnum = npml - j;
    xd = npml;
    xxn = xnum/xd;
    xn = (0.33)*power(xxn,3);
    gj2(j) = 1/(1+xn);
    gj2(JE-1-j) = 1/(1+xn);
    gj3(j) = (1-xn)/(1+xn);
    gj3(JE-j-1) = (1-xn)/(1+xn);
    xxn = (xnum-0.5)/xd;
    xn = (0.25)*power(xxn,3);
    fj1(j) = xn;
    fj1(JE-2-j) = xn;
    fj2(j) = 1/(1+xn);
    fj2(JE-2-j) = 1/(1+xn);
    fj3(j) = (1-xn)/(1+xn);
    fj3(JE-2-j) = (1-xn)/(1+xn);
end


%% Warning!! Don't change code from here!!
while (Nsteps > 0)
    n = 0;
    
    for n = 1:Nsteps                                                       %Nsteps is the number of times the main loop has executed
        T =T+1;                                                            %T keeps track of the timesteps
        %main fdtd loop
        
        %calculate the Dz field
        for j = 2:JE
                ez_inc(j) = ez_inc(j) + 0.5*(hx_inc(j-1)-hx_inc(j));
        end
        
        %pml boundary for the incident buffer
        ez_inc(1) = ez_inc_low_m2;
        ez_inc_low_m2 = ez_inc_low_m1;
        ez_inc_low_m1 = ez_inc(2);
        
        ez_inc(JE-1) = ez_inc_high_m2;
        ez_inc_high_m2 = ez_inc_high_m1;
        ez_inc_high_m1 = ez_inc(JE-2);
        
        %put pulse in the specified grid position
        pulse =  exp(-0.5*((to-T)/spread)^2);
        ez_inc(3) = pulse;
        
        %calculate incident Dz values
        for i = ia:ib
            dz(i,ja) = dz(i,ja) + (0.5*hx_inc(ja-1));
            dz(i,jb) = dz(i,jb) + (0.5*hx_inc(jb));
        end
        
            
        %calculate Ez field
        for j = 2:JE-1
            for i = 2:IE-1
                ez(i,j) = ga(i,j) * (dz(i,j) - iz(i,j));
                iz(i,j) = iz(i,j) + gb(i,j)*ez(i,j);
            end
        end
        
       
        
        for j = 1:JE
            hx_inc(j) = hx_inc(j) + 0.5*(ez_inc(j)-ez_inc(j+1));
        end
        
        %calculate the Hx field
        for j = 1:JE-1
            for i = 1:IE-1
                curl_e = ez(i,j) - ez(i,j+1);
                ihx(i,j) = ihx(i,j) + fi1(i)*curl_e;
                hx(i,j) = fj3(j)*hx(i,j) + fj2(j)*0.5*(curl_e + ihx(i,j));
            end
        end
        
        %incident Hx field values
        for i = ia:ib
            hx(i,ja-1) = hx(i,ja-1) + 0.5*ez_inc(ja);
            hx(i,jb) = hx(i,jb) - 0.5*ez_inc(jb);
        end
        
        %calculate the Hy field
        for j = 1:JE-1
            for i = 1:IE-1
                curl_e = ez(i+1,j) - ez(i,j);
                ihy(i,j) = ihy(i,j) + fj1(j)*curl_e;
                hx(i,j) = fj3(i)*hy(i,j) + fi2(i)*0.5*(curl_e + ihy(i,j));
            end
        end
        
        %incident Hy field values
        for i = ja:jb
            hy(ia-1,j) = hy(ia-1,j) - 0.5*ez_inc(j);
            hy(ib,j) = hy(ib,j) + 0.5*ez_inc(j);
        end
        
    end
    
    
    pause(0.2);
    fprintf('Timestep = %f \n',T);
    surf(ez);

end

