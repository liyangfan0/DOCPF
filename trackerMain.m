function [results] = trackerMain(p, im, bg_area, fg_area, area_resize_factor)
%TRACKERMAIN contains the main loop of the tracker, P contains all the parameters set in runTracker
    %% INITIALIZATION
    temp = load('w2crs');
    w2c = temp.w2crs;
    num_frames = numel(p.img_files);
    % used for OTB-13 benchmark
    OTB_rect_positions = zeros(num_frames, 4);
	pos = p.init_pos;
    target_sz = p.target_sz;
	num_frames = numel(p.img_files);
    % patch of the target + padding
    patch_padded = getSubwindow(im, pos, p.norm_bg_area, bg_area);
    % initialize hist model
    new_pwp_model = true;
    [bg_hist, fg_hist] = updateHistModel(new_pwp_model, patch_padded, bg_area, fg_area, target_sz, p.norm_bg_area, p.n_bins, p.grayscale_sequence);
    new_pwp_model = false;
    % Hann (cosine) window
    if isToolboxAvailable('Signal Processing Toolbox')
        hann_window = single(hann(p.cf_response_size(1)) * hann(p.cf_response_size(2))');
    else
        hann_window = single(myHann(p.cf_response_size(1)) * myHann(p.cf_response_size(2))');
    end
    % gaussian-shaped desired response, centred in (1,1)
    % bandwidth proportional to target size
    output_sigma = sqrt(prod(p.norm_target_sz)) * p.output_sigma_factor / p.hog_cell_size;
    y = gaussianResponse(p.cf_response_size, output_sigma);
    yf = fft2(y);
   %定义kalman的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A=[1,0,1,0;0,1,0,1;0,0,1,0;0,0,0,1];
    C=[1,0,0,0;0,1,0,0];
    mu=[pos,0,0]';
    Sigma=1000*eye(4);
    mu2=[pos,0,0]'; %定义第二个运动模型的参数
    Sigma2=1000*eye(4);
    %cQ=1e-5;
    cQ=1e-7;
    cR=1e-2;
    Q=cQ*[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    R=cR*[1,0;0,1];
    %cQ2=0;%只相信预测
    %cR2=1e+5;
    %Q2=cQ2*[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    %R2=cR2*[1,0;0,1];
    %设置APCE的阈值
    %APCE_eta=0.57;
    APCE_eta=0.1;
    APCE_ref=0;
    F=0;
    S1=[0,0];
    S2=[0,0];
    S=cell(1,2);
    flag=1;%定义目标的运动阶段，1代表正常无遮挡，2代表进入遮挡状态，3代表刚出遮挡
    %p.scale_adaptation=1;
    %% SCALE ADAPTATION INITIALIZATION
    if p.scale_adaptation
        % Code from DSST
        scale_factor = 1;
        base_target_sz = target_sz;
        scale_sigma = sqrt(p.num_scales) * p.scale_sigma_factor;
        ss = (1:p.num_scales) - ceil(p.num_scales/2);
        ys = exp(-0.5 * (ss.^2) / scale_sigma^2);
        ysf = single(fft(ys));
        if mod(p.num_scales,2) == 0
            scale_window = single(hann(p.num_scales+1));
            scale_window = scale_window(2:end);
        else
            scale_window = single(hann(p.num_scales));
        end;

        ss = 1:p.num_scales;
        scale_factors = p.scale_step.^(ceil(p.num_scales/2) - ss);

        if p.scale_model_factor^2 * prod(p.norm_target_sz) > p.scale_model_max_area
            p.scale_model_factor = sqrt(p.scale_model_max_area/prod(p.norm_target_sz));
        end

        scale_model_sz = floor(p.norm_target_sz * p.scale_model_factor);
        % find maximum and minimum scales
        min_scale_factor = p.scale_step ^ ceil(log(max(5 ./ bg_area)) / log(p.scale_step));
        max_scale_factor = p.scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ target_sz)) / log(p.scale_step));
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    num_particle=5;
    t_imread = 0;
    particle=cell(num_particle+1,1);%每一个cell{pos,particle_S1,particle_S2,mu,Sigma,flag,APCE_value}
    particle_S1=cell(num_particle+1,1);
    particle_S2=cell(num_particle+1,1);
    particle_S1{1}=[0,0];
    particle_S2{1}=[0,0];
   mu_particle=pos;%pos是第一帧的真实位置
   sigma_particle=[0.2,0;0,0.2];
   pos_particle=mvnrnd(mu_particle,sigma_particle,num_particle);
   %把pos也加上
   pos_particle=[pos;pos_particle];
   for i=1:num_particle
       particle{i}{1}=pos_particle(i,:);
       particle{i}{2}=[0,0];
       particle{i}{3}=[0,0];
       particle{i}{4}=mu;
       particle{i}{5}=Sigma;
       particle{i}{9}=mu2;
       particle{i}{10}=Sigma2;
       
   end
  
   
    %% MAIN LOOP
    tic;
    for frame = 1:num_frames
        if frame>1
            tic_imread = tic;
            im = imread([p.img_path p.img_files{frame}]);
            t_imread = t_imread + toc(tic_imread);
            
	    %% TESTING step
          
            for i=1:num_particle

                [particle{i}{1},particle{i}{2},particle{i}{3},particle{i}{4},particle{i}{5},particle{i}{6},particle{i}{7},particle{i}{8},particle{i}{9},particle{i}{10}]=...
                 my_response(particle{i}{4},particle{i}{5},Q,R,A,C,particle{i}{2},particle{i}{3},particle{i}{1},p,im,bg_area,area_resize_factor,w2c,hann_window,hf_num_hog,hf_num_cn,hf_den_hog,hf_den_cn,bg_hist,fg_hist,target_sz,particle{i}{9},particle{i}{10});
                %i
                %particle{i}{6}
                %particle{i}{7}
            end
            %[pos,particle_S1{1},particle_S2{1},mu,Sigma,flag,F]=my_response(mu,Sigma,Q,R,A,C,particle_S1{1},particle_S2{1},pos,p,im,bg_area,area_resize_factor,w2c,hann_window,hf_num_hog,hf_num_cn,hf_den_hog,hf_den_cn,bg_hist,fg_hist)  
            
             
            F_all=[];
            pos_all=[];
           
            for i=1:num_particle
                
                
                pos_all=[pos_all;particle{i}{1}];
                F_all=[F_all,particle{i}{7}];
            end
            
            %取最大值的情况
            %index=find(F_all==max(F_all(:)),1);
            %pos=pos_all(index,:);
            %取前三个粒子的值加权
            index=find(F_all==max(F_all(:)),2);
            pos_max3=pos_all(index,:);%[[180,48],[180,48],[180,48]]
            weight0=F_all(index);
            weight=weight0/sum(weight0);
            
            pos=weight*pos_max3;
            
            
           
            
            %由当前粒子的位置值predict粒子的下一个位置
            %进行重采样resample
            F_weight=F_all/sum(F_all);
            pos2_all=pos_all';
            S_hat=[pos2_all;F_weight];
            S = resample(S_hat);
            for i=1:num_particle
                S1=S(1:2,:);
                S2=S1';
                particle{i}{1}=S2(i,:);
            end 
            %%%%%%%%%%%%%
            %predict
            
            for i=1:num_particle
                particle{i}{1}=mvnrnd(particle{i}{1},sigma_particle,1);
            end    
            %for i=1:num_particle
            %    particle{i}{1}=particle{i}{8};
            %end    
                   
    
           
            
            rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
         
            
            
            %% SCALE SPACE SEARCH
            if p.scale_adaptation
                im_patch_scale = getScaleSubwindow(im, pos, base_target_sz, scale_factor * scale_factors, scale_window, scale_model_sz, p.hog_scale_cell_size);
                xsf = fft(im_patch_scale,[],2);
                scale_response = real(ifft(sum(sf_num .* xsf, 1) ./ (sf_den + p.lambda) ));
                recovered_scale = ind2sub(size(scale_response),find(scale_response == max(scale_response(:)), 1));
                %set the scale
                scale_factor = scale_factor * scale_factors(recovered_scale);

                if scale_factor < min_scale_factor
                    scale_factor = min_scale_factor;
                elseif scale_factor > max_scale_factor
                    scale_factor = max_scale_factor;
                end
                % use new scale to update bboxes for target, filter, bg and fg models
                target_sz = round(base_target_sz * scale_factor);
                avg_dim = sum(target_sz)/2;
                bg_area = round(target_sz + avg_dim);
                if(bg_area(2)>size(im,2)),  bg_area(2)=size(im,2)-1;    end
                if(bg_area(1)>size(im,1)),  bg_area(1)=size(im,1)-1;    end

                bg_area = bg_area - mod(bg_area - target_sz, 2);
                fg_area = round(target_sz - avg_dim * p.inner_padding);
                fg_area = fg_area + mod(bg_area - fg_area, 2);
                % Compute the rectangle with (or close to) params.fixed_area and
                % same aspect ratio as the target bboxgetScaleSubwindow
                area_resize_factor = sqrt(p.fixed_area/prod(bg_area));
            end

            if p.visualization_dbg==1
                mySubplot(2,1,5,1,im_patch_cf,'FG+BG','gray');
                mySubplot(2,1,5,2,likelihood_map,'obj.likelihood','parula');
                mySubplot(2,1,5,3,response_cf,'CF response','parula');
                mySubplot(2,1,5,4,response_pwp,'center likelihood','parula');
                mySubplot(2,1,5,5,response,'merged response','parula');
                drawnow
            end
        end

        %% TRAINING
        % extract patch of size bg_area and resize to norm_bg_area
        im_patch_bg = getSubwindow(im, pos, p.norm_bg_area, bg_area);
      
        % compute feature map, of cf_response_size
        xt = getFeatureMap(im_patch_bg, p.feature_type, p.cf_response_size, p.hog_cell_size,w2c);
        % apply Hann window
        xt = bsxfun(@times, hann_window, xt);
        % compute FFT
        xtf = fft2(xt);
        xtf_hog=xtf(:,:,1:28);
        xtf_cn=xtf(:,:,29:38);
        %% FILTER UPDATE
        % Compute expectations over circular shifts,
        % therefore divide by number of pixels.
		new_hf_num_hog = bsxfun(@times, conj(yf), xtf_hog) / prod(p.cf_response_size);
		new_hf_den_hog = (conj(xtf_hog) .* xtf_hog) / prod(p.cf_response_size);
        new_hf_num_cn = bsxfun(@times, conj(yf), xtf_cn) / prod(p.cf_response_size);
		new_hf_den_cn = (conj(xtf_cn) .* xtf_cn) / prod(p.cf_response_size);

        if frame == 1
            % first frame, train with a single image
		    hf_den_hog = new_hf_den_hog;
		    hf_num_hog = new_hf_num_hog;
            hf_den_cn = new_hf_den_cn;
		    hf_num_cn = new_hf_num_cn;
		else
		    % subsequent frames, update the model by linear interpolation
        	hf_den_hog = (1 - p.learning_rate_cf) * hf_den_hog + p.learning_rate_cf * new_hf_den_hog;
	   	 	hf_num_hog = (1 - p.learning_rate_cf) * hf_num_hog + p.learning_rate_cf * new_hf_num_hog;
            hf_den_cn = (1 - p.learning_rate_cf) * hf_den_cn + p.learning_rate_cf * new_hf_den_cn;
	   	 	hf_num_cn = (1 - p.learning_rate_cf) * hf_num_cn + p.learning_rate_cf * new_hf_num_cn;

            %% BG/FG MODEL UPDATE
            % patch of the target + padding
            [bg_hist, fg_hist] = updateHistModel(new_pwp_model, im_patch_bg, bg_area, fg_area, target_sz, p.norm_bg_area, p.n_bins, p.grayscale_sequence, bg_hist, fg_hist, p.learning_rate_pwp);
        end

        %% SCALE UPDATE
        if p.scale_adaptation
            im_patch_scale = getScaleSubwindow(im, pos, base_target_sz, scale_factor*scale_factors, scale_window, scale_model_sz, p.hog_scale_cell_size);
            xsf = fft(im_patch_scale,[],2);
            new_sf_num = bsxfun(@times, ysf, conj(xsf));
            new_sf_den = sum(xsf .* conj(xsf), 1);
            if frame == 1,
                sf_den = new_sf_den;
                sf_num = new_sf_num;
            else
                sf_den = (1 - p.learning_rate_scale) * sf_den + p.learning_rate_scale * new_sf_den;
                sf_num = (1 - p.learning_rate_scale) * sf_num + p.learning_rate_scale * new_sf_num;
            end
        end

        % update bbox position
        if frame==1, rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])]; end

        rect_position_padded = [pos([2,1]) - bg_area([2,1])/2, bg_area([2,1])];

        OTB_rect_positions(frame,:) = rect_position;

        if p.fout > 0,  fprintf(p.fout,'%.2f,%.2f,%.2f,%.2f\n', rect_position(1),rect_position(2),rect_position(3),rect_position(4));   end

        %% VISUALIZATION
        if p.visualization == 1
            if isToolboxAvailable('Computer Vision System Toolbox')
                im = insertShape(im, 'Rectangle', rect_position, 'LineWidth', 4, 'Color', 'black');
                im = insertShape(im, 'Rectangle', rect_position_padded, 'LineWidth', 4, 'Color', 'yellow');
                % Display the annotated video frame using the video player object.
                %step(p.videoPlayer, im);
            else
                figure(1)
                imshow(im)
                rectangle('Position',rect_position, 'LineWidth',2, 'EdgeColor','g');
                rectangle('Position',rect_position_padded, 'LineWidth',2, 'LineStyle','--', 'EdgeColor','b');
                %drawnow
            end
        end
        
    end
    
    elapsed_time = toc;
    % save data for OTB-13 benchmark
    results.type = 'rect';
    results.res = OTB_rect_positions;
    results.fps = num_frames/(elapsed_time - t_imread);
end

% Reimplementation of Hann window (in case signal processing toolbox is missing)
function H = myHann(X)
    H = .5*(1 - cos(2*pi*(0:X-1)'/(X-1)));
end

% We want odd regions so that the central pixel can be exact
function y = floor_odd(x)
    y = 2*floor((x-1) / 2) + 1;
end

function y = ensure_real(x)
    assert(norm(imag(x(:))) <= 1e-5 * norm(real(x(:))));
    y = real(x);
end
