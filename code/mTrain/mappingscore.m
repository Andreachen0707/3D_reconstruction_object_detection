function Y2 = mappingscore(calibratefun,confs)
    X =calibratefun(1,:); % scores
    Y =calibratefun(2,:); % precisions 

    X2 = confs(:);
    if numel(X) ~= numel(Y), error('size of X and Y must be the same.'); end

    % sort lookup table
    [~,idx] = sort(X,'ascend');
    X = X(idx);
    Y = Y(idx);

    % look up
    diff = bsxfun(@minus,X,X2) >= 0;
    diffdiff = imfilter(diff,[-1 1 0],'same','replicate');
    [J,I] = find(diffdiff);

    % linear interpolation
    Y2 = nan(numel(X2),1);
    if ~isempty(J),
        wLeft = X(I)' - X2(J);
        wRight = X2(J) - X(I-1)';
        Y2(J) = (Y(I-1)' .* wLeft + Y(I)' .* wRight) ./ (wLeft + wRight);
    end

    % add solution when X2 = X(1)
Y2(X2 <= X(1)) = Y(1);
Y2(X2 >= X(end)) = Y(end);
    
    % reshape Y2 to be consistent with X2
    Y2 = reshape(Y2,size(X2));
end