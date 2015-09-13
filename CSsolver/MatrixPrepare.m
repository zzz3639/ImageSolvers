function [A,b,c]=MatrixPrepare(Img, div,margin)
% Currently this function only work for square image.
% Modified from CSSTORM of paper "faster storm using compressed sensing"
    boxsize = size(Img,1);
    mag = 2 / boxsize;

    % generate grid for local estimation
    x_pos_est = linspace(-1,1, div*boxsize+1);
    y_pos_est = linspace(-1,1, div*boxsize+1);
    x_pos_est = x_pos_est(1:div*boxsize);
    y_pos_est = y_pos_est(1:div*boxsize);

    x_inx = x_pos_est((div/2) : div : end);
    y_inx = y_pos_est((div/2) : div : end);

    % provide pixel padding for local CS estimate
    % This should correspond to the box size around peak intensity
    % -- Begin padding the estimated local patch
    dx = x_pos_est(2)-x_pos_est(1);
    dy = y_pos_est(2)-y_pos_est(1);

    for mm = 1:margin
        x_pos_est = [x_pos_est(1)-dx x_pos_est];
        x_pos_est = [x_pos_est x_pos_est(end)+dx];
        y_pos_est = [y_pos_est(1)-dy y_pos_est];
        y_pos_est = [y_pos_est y_pos_est(end)+dy];
    end

    [x_inx y_inx] = meshgrid(x_inx,y_inx);
    [x_pos_est y_pos_est] = meshgrid(x_pos_est,y_pos_est);

    % generate measurement matrix
    len = length(x_pos_est(:));
    %A = zeros(x_inx, len + 1);
    for ii = 1:len
        img_kernel = MolKernel(x_inx, y_inx, x_pos_est(ii), y_pos_est(ii), mag);
        A(:,ii) = img_kernel(:);
    end
    
    c = sum(A);
    PSF_integ = max(c); % integration of the PSF over space, used for normalization
    c = c./ PSF_integ; % normalize to 1
    A = A./ PSF_integ;
    
    % add the extra optimization variable for the estimation of the background
    len = len + 1;
    c(len) = 0;
    A(:,len) = 1;
    b=Img(:);
end
