% Copyright (c) 2024 Li Lab, Southern University of Science and Technology, Shenzhen
% author: Jianwei Chen
% email: 12149038@mail.sustech.edu.cn
% date: 2024.11.28
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
%%
clear
close all
clc
%% hyper parameters for PSF model used for fit
lambda_all = [583, 658, 698];                                      % CF568, DY644, and CF680
paraSim.NA = 1.35;                                                % numerical aperture of obj             
paraSim.refmed = 1.40;                                            % refractive index of sample medium
paraSim.refcov = 1.518;                                           % refractive index of converslip
paraSim.refimm = 1.40;                                           % refractive index of immersion oil
paraSim.lambda = lambda_all(1);                                             % wavelength of emission
paraSim.objStage0_upper = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.objStage0_lower = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.zemit0_upper = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_upper);  % reference emitter z position, nm, distance of molecule to coverslip
paraSim.zemit0_lower = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_lower);  % reference emitter z position, nm, distance of molecule to coverslip


paraSim. pixelSizeX = 120;                                        % nm, pixel size of the image
paraSim. pixelSizeY = 120;                                        % nm, pixel size of the image
paraSim.Npupil = 64;                                             % sampling at the pupil plane

paraSim.aberrations(:,:,1) = [2,-2,0.0; 2,2,-0.1; 3,-1,0.0; 3,1,0.0; 4,0,0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,:,2) = [2,-2,0.0; 2,2,0.1; 3,-1,0.0; 3,1,0.0; 4,0,0.0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];

paraSim.aberrations(:,3,:) =  paraSim.aberrations(:,3,:)*paraSim.lambda;

paraSim.offset = [0 0];

%% parameters for molecules for simulation
Nmol = 101;
Npixels = 31;
paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;

paraSim.xemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.yemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.zemit = linspace(-600,600,Nmol)*1;                                      %nm
paraSim.objStage = linspace(-1000,1000,Nmol)*0;                                  %nm

phasediff = [-0.15,0,0.15];


for ii = 1: length(phasediff)
    paraSim.phaseshift = [0 ,pi/2, pi, 3*pi/2];
    paraSim.lambda =lambda_all(ii);
    [PSFs(:,:,:,:,ii),PSFsUpper,PSFsLower,WaberrationUpper, WaberrationLower] = vectorPSF_4Pi(paraSim);


    % imageslicer(PSFs(:,:,:,:,ii));


    ipalm_im  = PSFs(:,:,:,:,ii);

    phaseshift = paraSim.phaseshift;
    % k = 2 * pi*paraSim.refimm / paraSim.lambda; %lambda in nm
    k = 2 * pi / paraSim.lambda; %lambda in nm
    zcand = paraSim.zemit;% if move obj stage use paraSim.objStage
    zstep = zcand(2) - zcand(1);

    imsz = paraSim.sizeX;


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
     
    IAB(:,:,:,:,3,ii) =Spline3D_interp(I);
    IAB(:,:,:,:,1,ii) =Spline3D_interp(A);
    IAB(:,:,:,:,2,ii) =Spline3D_interp(B);
end
%% test fitter

dz = zcand(2) - zcand(1);
z0 = round(length(zcand) / 2);
% k = 2 * pi * paraSim.refimm / lambdanm; %lambda in nm
NV = 6;
Npixels = 13;
bg = [20 20 20 20];
Nfits = 10001;
zrange = linspace(-600,600,Nmol)*1;

load('Photons_CF568_DY634_AF647_CF680.mat');

[CF568_no,CF568_xo]=hist(CF568,1000);
distribution1=zeros(1000,2);
distribution1(:,1) = CF568_xo';
distribution1(:,2) = CF568_no';
N(:,1) = rand_arb_cjw(Nfits, distribution1, 0);

[DY634_no, DY634_xo]=hist(DY634,1000);
distribution1=zeros(1000,2);
distribution1(:,1) = DY634_xo';
distribution1(:,2) = DY634_no';
N(:,2) = rand_arb_cjw(Nfits, distribution1, 0);

[CF680_no, CF680_xo]=hist(CF680,1000);
distribution2=zeros(1000,2);
distribution2(:,1) = CF680_xo';
distribution2(:,2) = CF680_no';
N(:,3) = rand_arb_cjw(Nfits, distribution2, 0);



P = [];
CRLB = [];
LogL =[];
Ph = [];
CRLBh = [];
LogLh =[];


for jj = 1:size(N,2)

    phi0 = [0+phasediff(jj)/2 ,pi/2-phasediff(jj)/2, pi+phasediff(jj)/2, 3*pi/2-phasediff(jj)/2];


    ground_truth.x(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
    ground_truth.y(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
    ground_truth.bg =repmat(bg, [Nfits 1]);
    ground_truth.znm = linspace(-600,600,Nfits)';
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






    ground_truth.N  = [N(:,jj) N(:,jj) N(:,jj) N(:,jj)];
    true_theta = zeros(Nfits,length(phi0),NV);
    true_theta(:,:,1) = ground_truth.x;
    true_theta(:,:,2) = ground_truth.y;
    true_theta(:,:,3) = ground_truth.N;
    true_theta(:,:,4) = ground_truth.bg;
    true_theta(:,:,5) = repmat(ground_truth.zspline,1,4);
    true_theta(:,:,6) = repmat(ground_truth.phi,1,4);
    true_theta =permute(true_theta,[3 2 1]);



    PSF.Ispline = IAB(:,:,:,:,3,jj);
    PSF.Aspline = IAB(:,:,:,:,1,jj);
    PSF.Bspline = IAB(:,:,:,:,2,jj);
    

    shared=[1,1,1,1,1,1];

    CRLB_multicolor = calculate_CRLB_YL_shared(length(zrange), PSF, Npixels, phi0, true_theta,shared);
    CRLBx(:,ii)=sqrt(CRLB_multicolor(:,1))*paraSim. pixelSizeX;
    CRLBy(:,ii)=sqrt(CRLB_multicolor(:,2))*paraSim. pixelSizeY;
    CRLBz(:,ii)=sqrt(CRLB_multicolor(:,5))*dz;
    CRLBphi(:,ii)=sqrt(CRLB_multicolor(:,6))/(2*k);






    imstack(:,:,:,:,jj) = simSplinePSF(Npixels, PSF, ground_truth.N, ground_truth.bg, coordinates, phi0);%simulate images



    zstart = linspace(0,2*z0-1,11);
    numberchannel = 4;
    dTAll=zeros(Nfits,NV,numberchannel);
    for m =1:Nfits
        for i=1:2
            for j=2:numberchannel

                dTAll(m,i,j)= coordinates(m,i,j)-coordinates(m,i,1);
            end
        end
    end




    dTAll=permute(dTAll,[2 3 1]);
    initZA=repmat(zstart(1),[1 Nfits]);
    initPhase = linspace(0,2*pi,11);
    initPhaseA=repmat(initPhase(1),[1 Nfits]);


    shared=[1 1 1 1 1 1]; % shared = [x, y, NP, bg, z, phi]. when shared(i)=1, link all channels. when shared(i)=0, individual fit. when shared(i)=2, link channel 1 and 3.
    if shared(3)==0
        offset= 4*(2-(shared(3))-shared(4));
    else
        offset= 4*(2-(1/shared(3))-shared(4));
    end
    sharedA = repmat(shared',[1 Nfits]);


    for ii = 1: length(phasediff)
        phi0_fit = [0+phasediff(ii)/2 ,pi/2-phasediff(ii)/2, pi+phasediff(ii)/2, 3*pi/2-phasediff(ii)/2];

        phi0A=repmat( phi0_fit',[1 Nfits]);

        initZA=repmat(zstart(1),[1 Nfits]);
        initPhaseA=repmat(initPhase(1),[1 Nfits]);


        tic % GPU fitter

        [P(:,:,ii,jj),CRLB(:,:,ii,jj), LogL(:,ii,jj)] = GPUmleFit_LM_4Pi(single(imstack(:,:,:,:,jj)),uint32(sharedA),100,single(IAB(:,:,:,:,:,ii)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
        tGPU=toc;
        disp(['GPU fitter speed' num2str(Nfits/tGPU)])
        if length(zstart)>1
            for i=2:length(zstart)
                initZA=repmat(zstart(i),[1 Nfits]);
                initPhaseA=repmat(initPhase(i),[1 Nfits]);
                [Ph(:,:,ii,jj),CRLBh(:,:,ii,jj), LogLh(:,ii,jj)] = GPUmleFit_LM_4Pi(single(imstack(:,:,:,:,jj)),uint32(sharedA),100,single(IAB(:,:,:,:,:,ii)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
                indbetter = [];
                indbetter=squeeze(LogLh(:,ii,jj))-squeeze(LogL(:,ii,jj))>1e-4; %copy only everything if LogLh increases by more than rounding error.
                P(indbetter(:),:,ii,jj)=Ph(indbetter(:),:,ii,jj);
                CRLB(indbetter(:),:,ii,jj)=CRLBh(indbetter(:),:,ii,jj);
                LogL(indbetter(:),ii,jj)=LogLh(indbetter(:),ii,jj);
            end
        end
    end
end



%% find dye


[maxValues1, colIndices1]= max(LogL(:,:,1),[],2);


[maxValues2, colIndices2]= max(LogL(:,:,2),[],2);

[maxValues3, colIndices3]= max(LogL(:,:,3),[],2);

%% Figure


% Calculate the proportion of correct and incorrect classifications for each category
correct1 = sum(colIndices1 == 1) / length(colIndices1) * 100;
wrong1_2 = sum(colIndices1 == 2) / length(colIndices1) * 100;
wrong1_3 = sum(colIndices1 == 3) / length(colIndices1) * 100;

correct2 = sum(colIndices2 == 2) / length(colIndices2) * 100;
wrong2_1 = sum(colIndices2 == 1) / length(colIndices2) * 100;
wrong2_3 = sum(colIndices2 == 3) / length(colIndices2) * 100;

correct3 = sum(colIndices3 == 3) / length(colIndices3) * 100;
wrong3_1 = sum(colIndices3 == 1) / length(colIndices3) * 100;
wrong3_2 = sum(colIndices3 == 2) / length(colIndices3) * 100;



% composeing data
data = [correct1,   wrong1_2+wrong1_3;...
        correct2,   wrong2_1+wrong2_3;...
        correct3,   wrong3_1+wrong3_2;...
       ];




figure;
h = bar(data, 'stacked');


set(gca, 'XTickLabel', {'CF568','DY634',  'CD680'}, 'FontSize', 18,'FontWeight','bold');
set(gca, 'YTick', 0:25:100, 'YTickLabel', 0:25:100, 'FontSize', 18,'FontWeight','bold');


x_label = xlabel('True Category', 'FontSize', 12);
y_label = ylabel('Accuracy (%)', 'FontSize', 12);
ylim([0,120])



% add text
lgd=legend('Correct', 'Wrong', 'Location', 'NorthWest');

% show data
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        text(i, sum(data(i, 1:j)) - data(i, j)/2, sprintf('%0.1f%%', data(i, j)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w');
    end
end

set(lgd,'FontName','time','FontSize',12,'LineWidth',1,'FontWeight','bold');
set(x_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
