function out = getFeatureMap(im_patch, feature_type, cf_response_size, hog_cell_size,w2c)

% code from DSST

% allocate space
switch feature_type
    case 'fhog'
        im=im_patch;
        temp = fhog(single(im_patch), hog_cell_size);
        h = cf_response_size(1);
        w = cf_response_size(2);
        out = zeros(h, w, 28, 'single');
        out(:,:,2:28) = temp(:,:,1:27);
        if hog_cell_size > 1
            im_patch = mexResize(im_patch, [h, w] ,'auto');
        end
        % if color image
        if size(im_patch, 3) > 1
            im_patch = rgb2gray(im_patch);
        end
        out(:,:,1) = single(im_patch)/255 - 0.5;
        im_patch = imresize(im, [h w]);
        out_pca = get_feature_map(im_patch, 'cn', w2c);
        %size(out_pca)
        %out=out_pca(:,:,1:2);
    
        
        out=cat(3,out,out_pca);%这里不考虑灰度特征，融合10通道的CN特征
       
    case 'hogcolor'
        im=im_patch;
        temp = fhog(single(im_patch), hog_cell_size);
        h = cf_response_size(1);
        w = cf_response_size(2);
        out = zeros(h, w, 28, 'single');
        out(:,:,2:28) = temp(:,:,1:27);
        if hog_cell_size > 1
            im_patch = mexResize(im_patch, [h, w] ,'auto');
        end
        % if color image
        if size(im_patch, 3) > 1
            im_patch = rgb2gray(im_patch);
        end
        out(:,:,1) = single(im_patch)/255 - 0.5;
		%sz = size(out);
      
		%im_patch = imresize(im, [sz(1) sz(2)]);
		%out_npca = get_feature_map(im_patch, 'gray', w2c);%这里只是将图像resize到和hog特征图一样的大小，基本上是h和w缩小四倍，然后再简单提取hog特征
		%im_patch = imresize(im, [h w]);
        %out_pca = get_feature_map(im_patch, 'cn', w2c);
        

    
		%out=cat(3,out,out_pca);%这里不考虑灰度特征，融合10通道的CN特征
        
        
    case 'gray'
        if hog_cell_size > 1, im_patch = mexResize(im_patch,cf_response_size,'auto');   end
        if size(im_patch, 3) == 1
            out = single(im_patch)/255 - 0.5;
        else
            out = single(rgb2gray(im_patch))/255 - 0.5;
        end        
end
        
end

