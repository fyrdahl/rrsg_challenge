function [last_img,all_imgs,inital_img,delta] = do_sense_recon(rawdata,k,csm,nIter)

    % Contribution to the reproducible research study group initiative to reproduce [1]
    %
    % [1] Pruessmann, K. P.; Weiger, M.; Boernert, P. and Boesiger, P.
    % Advances in sensitivity encoding with arbitrary k-space trajectories.
    % Magn Reson Med 46: 638-651 (2001)
    %
    % By Alexander Fyrdahl (alexander.fyrdahl@ki.se)
    % Code made available for the ISMRM Reproducible Research Study Group
    %
    
    dims = size(csm);
    imSize = dims(1:2);
    delta = zeros(1,nIter);

    FT = NUFFT(k,abs(k),imSize,dims(end)); % NUFFT Operator
    inital_img = FT'*rawdata;

    E = SENSE(k,abs(k(:)),csm); % SENSE Operator

    % Conjugate gradient
    % ----------------------

    a = E'*rawdata;
    b = zeros(prod(imSize),1);
    p = a;
    r = a;

    for i = 1:nIter
        q = E'*(E*p);
        delta(i) = abs(sqrt(sum(E'*(E*b(:,i))-a)^2)/sqrt(sum(a(:))^2));
        b(:,i+1) = b(:,i) + (r(:,i)'*r(:,i))/(p'*q)*p;
        r(:,i+1) = r(:,i) - (r(:,i)'*r(:,i))/(p'*q)*q;
        p = r(:,i+1) + (r(:,i+1)'*r(:,i+1))/(r(:,i)'*r(:,i))*p;
    end
    
    % ----------------------
    
    all_imgs = reshape(b(:,2:end),[imSize(1) imSize(2) nIter]);
    last_img = all_imgs(:,:,end);
    
end
