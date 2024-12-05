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
lambda_all = [583, 665];                                          % CF568 and AF647
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
Nmol = 201;
Npixels = 31;
paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;

paraSim.xemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.yemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.zemit = linspace(-600,600,Nmol)*1;                                      %nm
paraSim.objStage = linspace(-1000,1000,Nmol)*0;                                  %nm

phasediff = [-0.15,0.15];
NP_all = [50,100,150,200,250,300,350,400,500,600,800,1000,1500,2000,3000,4000,5000];


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
Nfits = 2001;
zrange = linspace(-600,600,Nfits)*1;
P = [];
CRLB = [];
LogL =[];
Ph = [];
CRLBh = [];
LogLh =[];

for nn = 1:size(NP_all,2)

    N(:,1) = NP_all(nn)*ones(Nfits,1);
    N(:,2) = NP_all(nn)*ones(Nfits,1);

    for jj = 1:size(N,2)

        phi0 = [0+phasediff(jj)/2 ,pi/2-phasediff(jj)/2, pi+phasediff(jj)/2, 3*pi/2-phasediff(jj)/2];


        ground_truth.x(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
        ground_truth.y(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
        ground_truth.bg =repmat(bg, [Nfits 1]);
        %     ground_truth.znm = -500+1000*rand([Nfits 1]);
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


        zstart = linspace(0,2*z0-1,4);
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
        initPhase = linspace(0,2*pi,4);
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

            [P(:,:,ii,jj,nn),CRLB(:,:,ii,jj,nn), LogL(:,ii,jj,nn)] = GPUmleFit_LM_4Pi(single(imstack(:,:,:,:,jj)),uint32(sharedA),100,single(IAB(:,:,:,:,:,ii)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
            tGPU=toc;
            disp(['GPU fitter speed' num2str(Nfits/tGPU)])
            if length(zstart)>1
                for i=2:length(zstart)
                    initZA=repmat(zstart(i),[1 Nfits]);
                    initPhaseA=repmat(initPhase(i),[1 Nfits]);
                    [Ph(:,:,ii,jj,nn),CRLBh(:,:,ii,jj,nn), LogLh(:,ii,jj,nn)] = GPUmleFit_LM_4Pi(single(imstack(:,:,:,:,jj)),uint32(sharedA),100,single(IAB(:,:,:,:,:,ii)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
                    indbetter = [];
                    indbetter=squeeze(LogLh(:,ii,jj,nn))-squeeze(LogL(:,ii,jj,nn))>1e-4; %copy only everything if LogLh increases by more than rounding error.
                    P(indbetter(:),:,ii,jj,nn)=Ph(indbetter(:),:,ii,jj,nn);
                    CRLB(indbetter(:),:,ii,jj,nn)=CRLBh(indbetter(:),:,ii,jj,nn);
                    LogL(indbetter(:),ii,jj,nn)=LogLh(indbetter(:),ii,jj,nn);
                end
            end
        end
    end
end

%% find dye

for kk = 1:size(NP_all,2)
    [maxValues1(:,kk), colIndices1(:,kk)]= max(LogL(:,:,1,kk),[],2);
    correct1(kk) = sum(colIndices1(:,kk) == 1) / length(colIndices1(:,kk)) * 100;
    wrong1(kk) = sum(colIndices1(:,kk) ~= 1) / length(colIndices1(:,kk)) * 100;


    [maxValues2(:,kk), colIndices2(:,kk)]= max(LogL(:,:,2,kk),[],2);
    correct2(kk) = sum(colIndices2(:,kk) == 2) / length(colIndices2(:,kk)) * 100;
    wrong2(kk) = sum(colIndices2(:,kk) ~= 2) / length(colIndices2(:,kk)) * 100;


end

%% Figure


% composeing data
data = [correct1,   wrong1;...
        correct2,   wrong2;...
       ];

figure;
plot(NP_all,smooth(correct1,0.5,'loess'));
hold on
plot(NP_all,smooth(correct2,0.5,'loess'));
xlim([0,5000])
ylim([0,110])
yticks([0,25,50,75,100])
x_label = xlabel('Photo count in 20 bg');
y_label = ylabel('Accuracy (%)');
lgd = legend('Location', 'north');
lgd=legend('CF568', 'AF647');
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');
set(x_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',23,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',18,'FontWeight','bold');
legend('boxoff')
