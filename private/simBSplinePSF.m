function [out data] = simBSplinePSF(Npixels,b4pi,I,bg,cor,phi0)
t=tic;
if (nargin <5)
    error('Minimal usage: simSplinePSF(Npixels,coeff,I,bg,cor)');
end

Nfits = size(cor,1);
spline_xsize = b4pi.dataSize(1)  ;
spline_ysize = b4pi.dataSize(2);
spline_zsize = b4pi.dataSize(3);
spline_wsize = b4pi.dataSize(4);
off = floor(((spline_xsize+1)-Npixels)/2);
data = zeros(Npixels,Npixels,Nfits,'single');

for kk = 1:Nfits
    xcenter = cor(kk,1,:);
    ycenter = cor(kk,2,:);
    zcenter = cor(kk,3,:);
    
    xc = -1*(xcenter - Npixels/2+0.5);
    yc = -1*(ycenter - Npixels/2+0.5);
    zc = zcenter - floor(zcenter);
    
    xstart = floor(xc);
    xc = xc - xstart;
    
    ystart = floor(yc);
    yc = yc - ystart;
    
    
    zstart = floor(zcenter);
    
    
    phi=cor(kk,4,1);
    %     [delta_f1,delta_dxf1,delta_dyf1,delta_dzf1] = computeDelta3D_bSpline_v2_single(single(xc1),single(yc1),single(zc));
    %
    %     [delta_f2,delta_dxf2,delta_dyf2,delta_dzf2] = computeDelta3D_bSpline_v2_single(single(xc2),single(yc2),single(zc));
    for ss=1:length(phi0)
        Phi = wrapTo2Pi(phi+phi0(ss));
        wc=Phi./b4pi.dphidN+1;
        wstart=floor(wc);
        wc=wc-wstart;
        [delta_f,delta_dxf,delta_dyf,delta_dzf,delta_dwf] = Fan_computeDelta4D_bSpline(single(xc(ss)),single(yc(ss)),single(zc(ss)),single(wc));% new
        
        newTheta1 = [xc(ss),yc(ss),I(kk,ss),bg(kk,ss),zc(ss),wc];
        for ii = 0:Npixels-1
            for jj = 0:Npixels-1
                
%                 [~, model] = kernel_Derivative_bSpline4D_4pi(xstart+jj+off,ystart+ii+off,zstart,wT(1),spline_xsize,spline_ysize,spline_zsize,spline_wsize,b4pi.coeffs,delta_f,delta_dxf,delta_dyf,delta_dzf,delta_dwf,newTheta1);% nwe
                [~,model] = kernel_Derivative_bSpline4D_v2(xstart(ss)+jj+off,ystart(ss)+ii+off,zstart(ss),wstart,spline_xsize,spline_ysize,spline_zsize,spline_wsize,...
                    b4pi.coeffs,delta_f,delta_dxf,delta_dyf,delta_dzf,delta_dwf,newTheta1);

                data(ii+1, jj+1, kk, ss) = model;
                
            end
        end
        if toc(t)>1
            disp(kk/Nfits)
            t=tic;
        end
        
    end
    scale = [1 1 1 1];
    
    
for p = 1:1:length(phi0)
    out(:, :, :, p) = (poissrnd(data(:, :, :, p)*scale(p),Npixels,Npixels,Nfits)); 
   %out = data;
end
end


