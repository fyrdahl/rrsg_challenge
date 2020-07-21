classdef SENSE
    
    % Contribution to the reproducible research study group initiative to reproduce [1]
    %
    % [1] Pruessmann, K. P.; Weiger, M.; Boernert, P. and Boesiger, P.
    % Advances in sensitivity encoding with arbitrary k-space trajectories.
    % Magn Reson Med 46: 638-651 (2001)
    %
    % By Alexander Fyrdahl (alexander.fyrdahl@ki.se)
    % Code made available for the ISMRM Reproducible Research Study Group
    %
    
    properties
        w
        st
        csm
        conj_csm
        nCh
        imSize
        adjoint
        dataSize
    end
    
    methods
        
        function obj = SENSE(k,w,csm)
            
            k = double(k);
            csm = double(csm);
            dims = size(csm);
            Nd = dims(1:end-1);
            Jd = repmat(6,[1 length(dims)-1]);
            Kd = floor(Nd*1.5);
            om = [real(k(:)), imag(k(:))]*2*pi;
            
            obj.st = nufft_init(om, Nd, Jd, Kd, Nd/2,'kaiser');
            obj.w = w;
            obj.csm = csm;
            obj.conj_csm = conj(csm);
            obj.nCh = size(obj.csm, 3);
            obj.imSize = Nd;
            obj.adjoint = 0;
            obj.dataSize = size(k);
        end
        
        function out = mtimes(obj,in)
            
            % Reproduce algorithm from [1]

            if obj.adjoint
                
                kspace = in; % Input is k-space
                
                % Density correction (D)
                kspace = reshape(kspace, [], obj.nCh).* repmat(obj.w,[1 obj.nCh]);
                
                % Apply E^H
                % ----------------
                
                % Fourier Transform to image domain (FT1)
                img = nufft_adj( kspace, obj.st )/sqrt(prod(obj.imSize));
                
                % Multiply by complex conjugate of sensitivity (S*)
                img = img .* obj.conj_csm;
                
                % Complex sum of input images (SUM)
                img = sum(img, 3);
                
                % ----------------

                out = img(:);
                
            else
                
                img = in;
                
                % Apply E
                % ----------------

                % Multiply by sensitivity (S)
                img = repmat(reshape(img,[obj.imSize(1),obj.imSize(2)]),[1 1 obj.nCh]) .* obj.csm;
                
                % Fourier Transform to k-space (FT2)
                kspace = nufft(img, obj.st)/sqrt(prod(obj.imSize));
                
                out = kspace(:);
                
                % ----------------

                
            end
        end

        function obj = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
        end
        
        function obj = transpose(obj)
            obj.adjoint = xor(obj.adjoint,1);
        end
        
    end
end

