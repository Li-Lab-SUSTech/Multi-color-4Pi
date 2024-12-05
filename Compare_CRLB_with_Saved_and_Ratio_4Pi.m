% Copyright (c) 2024 Li Lab, Southern University of Science and Technology, Shenzhen
% author: Jianwei Chen 
% email: 12149038@mail.sustech.edu.cn
% date: 2024.11.28
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
clear
close all
clc
%% hyper parameters for PSF model used for fit
times = 2000;
paraSim.NA = 1.35;                                                % numerical aperture of obj             
paraSim.refmed = 1.40;                                            % refractive index of sample medium
paraSim.refcov = 1.518;                                           % refractive index of converslip
paraSim.refimm = 1.40;                                           % refractive index of immersion oil
paraSim.lambda = 680;                                             % wavelength of emission  568 & 647
paraSim.objStage0_upper = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.objStage0_lower = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.zemit0_upper = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_upper);  % reference emitter z position, nm, distance of molecule to coverslip
paraSim.zemit0_lower = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_lower);  % reference emitter z position, nm, distance of molecule to coverslip


paraSim. pixelSizeX = 120;                                        % nm, pixel size of the image
paraSim. pixelSizeY = 120;                                        % nm, pixel size of the image
paraSim.Npupil = 64;                                             % sampling at the pupil plane

paraSim.aberrations(:,:,1) = [2,-2,0.0; 2,2,-0.06; 3,-1,0.0; 3,1,0.0; 4,0,0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,:,2) = [2,-2,0.0; 2,2,0.06; 3,-1,0.0; 3,1,0.0; 4,0,0.0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];

paraSim.aberrations(:,3,:) =  paraSim.aberrations(:,3,:)*paraSim.lambda;

paraSim.offset = [0 0];
paraSim.phaseshift = [0 ,pi/2, pi, 3*pi/2];

%% parameters for molecules for simulation
Nmol = 201;
Npixels = 31;
% Nphotons = 5000 +10000*rand(1,Nmol);
% bg = 10 + 10*rand(1,Nmol);
paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;

paraSim.xemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.yemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.zemit = linspace(-1000,1000,Nmol)*1;                                      %nm
paraSim.objStage = linspace(-1000,1000,Nmol)*0;                                  %nm

[PSFs,PSFsUpper,PSFsLower,WaberrationUpper, WaberrationLower] = vectorPSF_4Pi(paraSim);

% imageslicer(PSFs);

%% generate IAB model


%phaseshift, lambdanm and zcand (z range in nm) are 
%saved in the file with the zstack
ipalm_im  = PSFs;

phaseshift = paraSim.phaseshift;
% k = 2 * pi*paraSim.refimm / paraSim.lambda; %lambda in nm
k = 2 * pi / paraSim.lambda; %lambda in nm
zcand = paraSim.zemit;% if move obj stage use paraSim.objStage
zstep = zcand(2) - zcand(1);

imsz = paraSim.sizeX;

%PSF.ipalm_im = PSF.ipalm_im * 4 / sum(reshape(PSF.ipalm_im(:, :, round(length(PSF.zcand) / 2), :), 1, [])); %normalize PSF
%imageslicer(PSF.ipalm_im);

%make I, A, B and their splines
I = squeeze((ipalm_im(:, :, :, 1) + ipalm_im(:, :, :, 3)) / 2);


kz2 = 2 * k * zcand';
kz2 = permute(repmat(kz2, 1, imsz, imsz), [2, 3, 1]);

F_phi1 = squeeze(ipalm_im(:, :, :, 1)) - I;
F_phi2 = squeeze(ipalm_im(:, :, :, 2)) - I;
phi1 = phaseshift(1);
phi2 = phaseshift(2);

A = (F_phi1 .* sin(kz2 + phi2) - F_phi2 .* sin(kz2 + phi1)) / sin(phi2 - phi1);
B = (-F_phi1 .* cos(kz2 + phi2) + F_phi2 .* cos(kz2 + phi1)) / sin(phi2 - phi1);


%check: restore PSF in the 4th quadrant
phi4 = phaseshift(4);
check_PSF = I + A .* cos(kz2 + phi4) + B .* sin(kz2 + phi4);
% imageslicer(squeeze(PSF.ipalm_im(:, :, :, 4)));
% imageslicer(check_PSF);
imageslicer([check_PSF;PSFs(:,:,:,4);(check_PSF-PSFs(:,:,:,4))]);
Ispline = Spline3D_interp(I);
Aspline = Spline3D_interp(A);
Bspline = Spline3D_interp(B);

%% test fitter
lambdanm = paraSim.lambda;
dz = zcand(2) - zcand(1);
z0 = round(length(zcand) / 2);
% k = 2 * pi * paraSim.refimm / lambdanm; %lambda in nm
NV = 6;
PSF.Ispline = Ispline;
PSF.Aspline = Aspline;
PSF.Bspline = Bspline;
Npixels = 13;
bg = [20 20 20 20];
zrange = paraSim.zemit;
Nfits = length(zrange);
ratio1=1;
ratio2=1;
NP=1000;
phasediff= [0];


ratio1_all = [0.716,1,1];
ratio2_all = [0.716,0.538,1];

for ii = 1:size(phasediff,2)
    for jj = 1:size(ratio1_all,2)
        phi0=[0+phasediff(ii)/2 ,pi/2-phasediff(ii)/2, pi+phasediff(ii)/2, 3*pi/2-phasediff(ii)/2];
        Nphotons = [NP*ratio1_all(jj) NP*ratio2_all(jj) NP*ratio1_all(jj) NP*ratio2_all(jj)];

        ground_truth.x(:,1) = Npixels/2 -1 +zrange*0;
        ground_truth.y(:,1) = Npixels/2 -1 +zrange*0;
        ground_truth.N = repmat(Nphotons, [Nfits 1]);
        ground_truth.bg =repmat(bg, [Nfits 1]);
        ground_truth.znm = zrange';
        ground_truth.zspline = ground_truth.znm / dz + z0;
        ground_truth.phi =  wrapTo2Pi(2 * k * ground_truth.znm);

        scale =0;
        ground_truth.x(:,2)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
        ground_truth.y(:,2)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;

        ground_truth.x(:,3)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
        ground_truth.y(:,3)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;

        ground_truth.x(:,4)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
        ground_truth.y(:,4)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;

        % coordinates for simulation
        coordinates = zeros(Nfits,4,length(phi0));
        for kk=1:1:length(phi0)
            coordinates(:,:,kk) = [ground_truth.x(:,kk) ground_truth.y(:,kk) ground_truth.zspline ground_truth.phi];
        end

        true_theta = zeros(Nfits,length(phi0),NV);
        true_theta(:,:,1) = ground_truth.x;
        true_theta(:,:,2) = ground_truth.y;
        true_theta(:,:,3) = ground_truth.N;
        true_theta(:,:,4) = ground_truth.bg;
        true_theta(:,:,5) = repmat(ground_truth.zspline,1,4);
        true_theta(:,:,6) = repmat(ground_truth.phi,1,4);
        true_theta =permute(true_theta,[3 2 1]);

        %

        shared=[1,1,1,1,1,1];

        CRLB_multicolor = calculate_CRLB_YL_shared(length(zrange), PSF, Npixels, phi0, true_theta,shared);
        CRLBx(:,ii,jj)=sqrt(CRLB_multicolor(:,1))*paraSim. pixelSizeX;
        CRLBy(:,ii,jj)=sqrt(CRLB_multicolor(:,2))*paraSim. pixelSizeY;
        CRLBz(:,ii,jj)=sqrt(CRLB_multicolor(:,5))*dz;
        CRLBphi(:,ii,jj)=sqrt(CRLB_multicolor(:,6))/(2*k);
    end

end

 %% plot

figure,
h1 = plot(zrange, CRLBx(:,1,1));
for ii =2:size(ratio1_all,2)
    hold on
    plot(zrange, CRLBx(:,1,ii));
end
x_label = xlabel('Z (nm)');
y_label = ylabel('CRLB X (nm)');
title('DY634', 'FontSize', 23,'LineWidth',1,'FontWeight','bold');
lgd = legend('Location', 'north');
axis ( [-600 600 0 12] );
lgd=legend("Salvaged", "Ratiometric", "Phase-based");
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');
set(x_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',18,'FontWeight','bold');
legend('boxoff')

figure,
h1 = plot(zrange, CRLBy(:,1,1));
for ii =2:size(ratio1_all,2)
    hold on
    plot(zrange, CRLBy(:,1,ii));
end
x_label = xlabel('Z (nm)');
y_label = ylabel('CRLB Y (nm)');
title('DY634', 'FontSize', 23,'LineWidth',1,'FontWeight','bold');
lgd = legend('Location', 'north');
axis ( [-600 600 0 12] );
lgd=legend("Salvaged", "Ratiometric", "Phase-based");
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');
set(x_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',18,'FontWeight','bold');
legend('boxoff')

figure,
h1 = plot(zrange, CRLBphi(:,1,1));
for ii =2:size(ratio1_all,2)
    hold on
    plot(zrange, CRLBphi(:,1,ii));
end
x_label = xlabel('Z (nm)');
y_label = ylabel('CRLB Z (nm)');
title('DY634', 'FontSize', 23,'LineWidth',1,'FontWeight','bold');
lgd = legend('Location', 'north');
axis ( [-600 600 0 8] );
lgd=legend("Salvaged", "Ratiometric", "Phase-based");
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');
set(x_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',18,'FontWeight','bold');
legend('boxoff')
