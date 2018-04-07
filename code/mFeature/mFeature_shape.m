function feature = mFeature_shape(bincellDim,cellIdx,xyzShifted,Space,Qnum,KNN,featurepara,raw)

load('./featurecenters/noNomalizedcenter/featurecenter_shape.mat','shape_center');

supergridSize =2;
if exist('raw','var')&&raw
    feature = zeros([bincellDim,84]);
else 
    raw =0;
    feature = zeros([bincellDim,size(shape_center,1)]);
end
[ind1,ind2,ind3]=ndgrid(1:size(feature,1),1:size(feature,2),1:size(feature,3));
for i=1:length(ind1(:))
    ptsIncellidx = cellIdx(:,2)==ind1(i)&cellIdx(:,1)== ind2(i)&cellIdx(:,3)==ind3(i);
    ptsIncell = xyzShifted(ptsIncellidx,:);
    shapefea_raw = getShapefea(ptsIncell,ind2(i),ind1(i),ind3(i),Space.s,supergridSize*Space.s/Qnum,Qnum/supergridSize);
    if featurepara.normalize||featurepara.highdim
        shapefea_norm = (shapefea_raw(:)-norm_mean(:))./norm_scale(:);
    else
        shapefea_norm =shapefea_raw;
        stdf =4e-04;
    end
    if raw
        feature(ind1(i),ind2(i),ind3(i),1:end)  =shapefea_raw(:);
    else
        shapefea_code = sparse_coding(shapefea_norm',shape_center,stdf,KNN);
        feature(ind1(i),ind2(i),ind3(i),1:end)  =shapefea_code(:);
    end
end
end

function fea =getShapefea(pts,x,y,z,cellsize,stepSize,Gridnum)
        fea = get_shape_feature(pts');
        xyzShifted = bsxfun(@minus,pts,([x,y,z]-1)*cellsize);
        gridIdx = floor(bsxfun(@rdivide,xyzShifted,[stepSize,stepSize,stepSize])) + 1; % matlab index starts from 1
        [ind1,ind2,ind3]=ndgrid(1:Gridnum,1:Gridnum,1:Gridnum);
        for i=1:length(ind1(:))
            ptsIncellidx = gridIdx(:,2)==ind1(i)&gridIdx(:,1)== ind2(i)&gridIdx(:,3)==ind3(i);
            ptsIncell = pts(ptsIncellidx,:);
            shapefea= get_shape_feature(ptsIncell');
            fea =[fea;shapefea];
        end
end

function shapefea= get_shape_feature(data)
%data 3*N
    if size(data,2)<5
        shapefea =zeros(3,1);
    else
        [~,N] = size(data); 
        mn = mean(data,2); 
        data = data - repmat(mn,1,N); 
        covariance = 1 / (N-1) * data * data'; 
        [PC, V] = eig(covariance); 
        V = real(diag(V)); 
        % sort the variances in decreasing order 
        [~, rindices] = sort(-1*V); 
        V = V(rindices);
        C =PC(:,rindices);
        r1 =V(1);
        r2 =(V(1)-V(2));
        r3 =(V(2)-V(3));
        shapefea =[r1;r2;r3];

    end
end