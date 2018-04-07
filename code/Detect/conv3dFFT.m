function C = conv3dFFT(A,B)
    sizeA = size(A(:,:,:,1));
    sizeB = size(B(:,:,:,1));
    sizeConv = sizeA + sizeB - 1;
    C = zeros(sizeA - sizeB + 1);
    for dim = 1:size(A,4),
        % pad, multiply and transform back
        CPad = ifftn(fftn(A(:,:,:,dim),sizeConv).* fftn(B(:,:,:,dim),sizeConv));

        % frequency-domain convolution result
        C = C + CPad(sizeB(1):end-sizeB(1)+1,sizeB(2):end-sizeB(2)+1,sizeB(3):end-sizeB(3)+1); 
    end
end





