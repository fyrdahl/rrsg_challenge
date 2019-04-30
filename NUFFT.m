classdef NUFFT
    
    % Matlab wrapper for the NUFFT toolbox from Jeff Fessler.
    % Contributors: Miki Lustig, Li Feng, Ricardo Otazo,
    % Claudio Santelli and Alexander Fyrdahl
    %
    % Inputs:
    %		k      - normalized kspace coordinates (-0.5 to 0.5)
    %`		w      - density compensation
    %		imSize - the image size
    %		nc     - number of coils for multicoil gridding [optional]
    %
    %	Outputs:
    %		obj    - NUFFT operator
    %
    % Classdef syntax by Alexander Fyrdahl (alexander.fyrdahl@ki.se)
    % Code made available for the ISMRM Reproducible Research Study Group
    %
    
    properties
        k
        w
        st
        nCoils
        imSize
        adjoint
        dataSize
    end
    
    methods
        function obj = NUFFT(k,w,imSize,nc)
            Nd = imSize;
            Jd = [6,6];
            Kd = floor(Nd*1.5);
            om = [real(k(:)), imag(k(:))]*2*pi;
            if nargin < 4; nc = 1; end 
            obj.k = k;
            obj.w = w(:);
            obj.st = nufft_init(om, Nd, Jd, Kd, Nd/2,'kaiser');
            obj.nCoils = nc;
            obj.imSize = imSize;
            obj.adjoint = 0;
            obj.dataSize = size(k);
        end
        
        % Operator overloading
        function res = mtimes(obj,b)
            if obj.adjoint
                res = reshape(nufft_adj(reshape(b,[],obj.nCoils).*repmat(obj.w(:),[1 obj.nCoils]),obj.st)/sqrt(prod(obj.imSize)),obj.imSize(1),obj.imSize(2),obj.nCoils);
            else
                res = reshape(nufft(b,obj.st)/sqrt(prod(obj.imSize)),obj.dataSize(1),obj.dataSize(2));
            end
        end

        function obj = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
        end
        
        function obj = transpose(obj)
            obj.adjoint = xor(obj.adjoint,1);
        end
        
        % Optional methods go here
        function w = calculate_dcf(obj)
            w = abs(1 ./ (obj.st.p * (obj.st.p' * ones(size(obj.k(:))))));
        end
    end
end

