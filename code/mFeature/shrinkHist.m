function hist2 = shrinkHist(hist1,stepSize)
    if numel(stepSize) == 1, stepSize = repmat(stepSize,[1 ndims(hist1)]); end
    nSteps = size(hist1) ./ stepSize;
    for i = 1:ndims(hist1),
        if stepSize ~= round(stepSize), error('step size must be an integer!'); end
        if nSteps(i) ~= round(nSteps(i)), error('histogram size must be multiples of step size!'); end
    end
    if ndims(hist1) ~= 3, error('histogram must be 3D!'); end
    
    [X,Y,Z] = ndgrid(1:size(hist1,1),1:size(hist1,2),1:size(hist1,3));
    idx = sub2ind(size(hist1),X(:),Y(:),Z(:));
    X = ceil(X(:)/stepSize(1));
    Y = ceil(Y(:)/stepSize(2));
    Z = ceil(Z(:)/stepSize(3));
    hist2 = accumarray([X,Y,Z],hist1(idx),nSteps);
end