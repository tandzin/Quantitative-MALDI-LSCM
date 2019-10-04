% ----------------------------------------------------------------------------------------------------------
% Author: 	Jan Michálek,Dr.sc.techn.ETH
% 			Centre for Biomedical Image Analysis
% 			Faculty of Informatics, Masaryk University
% 			Botanicka 68a
% 			Brno, 602 00
% 			Czech Republic
% 			jan.michalek@fi.muni.cz
% ----------------------------------------------------------------------------------------------------------
% The algorithm is based on the theory presented in the paper:
%
% Michálek, Jan, Karel Št?pka, Michal Kozubek, Jarmila Navrátilová, Barbora Pavlatovská, Markéta Machálková, Jan Preisler, and Adam Pruška. n.d.
% 'Quantitative Assessment of Anti-Cancer Drug Efficacy From Coregistered Mass Spectrometry and Fluorescence Microscopy Images of Multicellular Tumor Spheroids.'
% Microscopy and Microanalysis. Cambridge University Press, 1–12. doi:10.1017/S1431927619014983.
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
% The program flow in the 'MALDI_IHC_peels_github.m' MATLAB script matches the visual flowchart in Fig.6 of the paper 
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
% To run the demos in the script
%		-> download the directory 'github' to a location of your choice
%		-> launch MATLAB (the script was developed with MATLAB R2016b)
%		-> change directory to 'github'
%		-> type 'MALDI_IHC_peels_github' on the command line
%		-> select one of the 3 demo datasets
%		-> results will be found in the subdirectory '4_Results' which will open in 
%          MATLAB's 'Current folder' and its subdirectories after the script terminates
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
% Example of the MATLAB dialog:
% 
%         >> cd 'I:\Topics\M&M 2019\revision 1\github'
%         >> MALDI_IHC_peels_github
%         image_nr [1...3]:2
%         LSCM_R_name=Spheroid_red.tif
%         LSCM_B_name=Spheroid_blue.tif
%         LSCM_T_name=MAX_Spheroid_transm.tif
%         ElapsedTime=66.9235
%         >>
% 
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
% If you wish to run the script with your own MALDI/LSCM data, you may place it in additional subdirectories (e.g. s4,s5..) 
% and extend the data selection dialog in the script 'Select_3_MALDI_IHC_github.m'
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
addpath(genpath(cd))
close all; clear variables;clear all;
system_dependent('precision', '24')
warning off images:initSize:adjustingMag
warning off MATLAB:MKDIR:DirectoryExists
warning off MATLAB:print:FigureTooLargeForPage
fullscreen = get(0,'ScreenSize');
screen_left=fullscreen(1);
screen_bottom=fullscreen(2);
screen_width=fullscreen(3);
screen_height=fullscreen(4);
Select_3_MALDI_IHC_github
fullscreen = get(0,'ScreenSize');
cd('3_TIFF+PNG_export')

% MALDI perifosine image 
MALDI_perifosine_float=imread(MALDI_perifosine_name);

% INPUT MALDI FIDUCIAL
MALDI_fiducial=single(rgb2gray(imread(MALDI_fiducial_name,'png')));% convert RGB to gray
% BINARIZE MALDI FIDUCIAL
MALDI_fiducial_binary=single(MALDI_fiducial>max(max(MALDI_fiducial))*MALDI_fiducial_threshold);% thresholding
figure('Name','MALDI_fiducial_binary','NumberTitle','off');imshow(MALDI_fiducial_binary,[]);

maldi_se=strel('disk',1)';
MALDI_fiducial_binary=imclose(MALDI_fiducial_binary,maldi_se);% morphological closing of dark spots in the fiducials
figure('Name','MALDI_fiducial_closed','NumberTitle','off');imshow(MALDI_fiducial_binary,[]);% binarized MALDI fiducial


% INPUT THE LSCM SPHEROID BOUNDARY MASK
LSCM_mask=single(imrotate(imread(LSCM_mask_name,'tif'),90)>0);

[H W]=size(LSCM_mask);
moving_mask=LSCM_mask;% prefix moving_ denotes all LSCM images that will be aligned with the fixed MALDI images

N_graphs=1; % perifosine means in the peels will be always plotted
if ~isempty (LSCM_R_name);% if red LSCM image is provided
    N_graphs=N_graphs+1;% then means of red signal will be plotted
    InfoImage=imfinfo(LSCM_R_name);
    LSCM_Img_R=zeros(InfoImage(1).Width,InfoImage(1).Height,length(InfoImage),'single');
    SUM_LSCM_Img_R=zeros(InfoImage(1).Width,InfoImage(1).Height,'single');
    for i=1:length(InfoImage)
        % INPUT CONFOCAL SLICES OF THE LSCM RED SIGNAL
        LSCM_Img_R(:,:,i)=single(imrotate(imread(LSCM_R_name,'Index',i),90));
        SUM_LSCM_Img_R=SUM_LSCM_Img_R+single(LSCM_Img_R(:,:,i));
    end
    disp(['LSCM_R_name=' LSCM_R_name]);
    LSCM_Img_R_spheroid=SUM_LSCM_Img_R/length(InfoImage);
    moving_red=LSCM_Img_R_spheroid;% will be aligned with the fixed MALDI images
    z_stack_length=length(InfoImage);
end;
if ~isempty (LSCM_G_name);% if green LSCM image is provided
    N_graphs=N_graphs+1;% then means of green signal will be plotted
    InfoImage=imfinfo(LSCM_G_name);
    LSCM_Img_G=zeros(InfoImage(1).Width,InfoImage(1).Height,length(InfoImage),'single');
    SUM_LSCM_Img_G=zeros(InfoImage(1).Width,InfoImage(1).Height,'single');
    for i=1:length(InfoImage)
        % INPUT CONFOCAL SLICES OF THE LSCM GREEN SIGNAL
        LSCM_Img_G(:,:,i)=single(imrotate(imread(LSCM_G_name,'Index',i),90));
        SUM_LSCM_Img_G=SUM_LSCM_Img_G+single(LSCM_Img_G(:,:,i));
    end
    disp(['LSCM_G_name=' LSCM_G_name]);
    LSCM_Img_G_spheroid=SUM_LSCM_Img_G/length(InfoImage);
    moving_green=LSCM_Img_G_spheroid;% will be aligned with the fixed MALDI images
    z_stack_length=length(InfoImage);
end;

if ~isempty (LSCM_B_name);% if blue LSCM image is provided
    N_graphs=N_graphs+1;% then means of blue signal will be plotted
    InfoImage=imfinfo(LSCM_B_name);
    LSCM_Img_B=zeros(InfoImage(1).Width,InfoImage(1).Height,length(InfoImage),'single');
    SUM_LSCM_Img_B=zeros(InfoImage(1).Width,InfoImage(1).Height,'single');
    for i=1:length(InfoImage)
        % INPUT CONFOCAL SLICES OF THE LSCM BLUE SIGNAL
        LSCM_Img_B(:,:,i)=single(imrotate(imread(LSCM_B_name,'Index',i),90));
        SUM_LSCM_Img_B=SUM_LSCM_Img_B+single(LSCM_Img_B(:,:,i));
    end
    disp(['LSCM_B_name=' LSCM_B_name]);
    LSCM_Img_B_spheroid=SUM_LSCM_Img_B/length(InfoImage);
    moving_blue=LSCM_Img_B_spheroid;% will be aligned with the fixed MALDI images
    z_stack_length=length(InfoImage);
end;
% if ~isempty (LSCM_T_name); we always have transmission
InfoImage=imfinfo(LSCM_T_name);
LSCM_Img_T=zeros(InfoImage(1).Width,InfoImage(1).Height,length(InfoImage),'single');
SUM_LSCM_Img_T=zeros(InfoImage(1).Width,InfoImage(1).Height,'single');
for i=1:length(InfoImage)
    % INPUT PATCHED CONFOCAL SLICES OF THE LSCM FIDUCIALS
    LSCM_Img_T(:,:,i)=single(imrotate(imread(LSCM_T_name,'Index',i),90));
    SUM_LSCM_Img_T=SUM_LSCM_Img_T+single(LSCM_Img_T(:,:,i));
end
LSCM_Img_T_spheroid=SUM_LSCM_Img_T/length(InfoImage);
disp(['LSCM_T_name=' LSCM_T_name]);

% INVERT AND THRESHOLD THE LSCM FIDUCIALS IN A SINGLE OPERATION
LSCM_fiducial_threshold=LSCM_fiducial_threshold/z_stack_length;
LSCM_fiducial=single((LSCM_Img_T_spheroid<LSCM_fiducial_threshold)&(LSCM_Img_T_spheroid>0));
figure('Name','LSCM_fiducial','NumberTitle','off');imshow(LSCM_fiducial,[]);
LSCM_fiducial_closed = imclose(LSCM_fiducial,se);% morphological closing of dark spots in the fiducials
LSCM_fiducial_open = imopen(LSCM_fiducial_closed,se);% morphological opening of white spots in the background
LSCM_fiducial_binary=single(LSCM_fiducial_open);
figure('Name','LSCM_fiducial_binary','NumberTitle','off');imshow(LSCM_fiducial_binary,[]);

moving_fiducial=LSCM_fiducial_binary;% LSCM fiducial will be aligned with the fixed MALDI images
% MATCH  LSCM  AND  MALDI  FIDUCIAL  CENTROIDS
% find the centroid of the moving fiducial image
[r_cntrd_mov, c_cntrd_mov]=find_image_centroid(moving_fiducial);
% find in the LSCM fiducial image the x-range and the y-range beyond
% ... which there are no white pixels
x_proj_moving=single(sum(moving_fiducial,1)~=0);
x_proj_moving=single(max(moving_fiducial,[],1));
x_proj_der_moving=abs(x_proj_moving(1:length(x_proj_moving)-1)-x_proj_moving(2:length(x_proj_moving)));
x_limits_moving=find(x_proj_der_moving);
x_range_moving=x_limits_moving(length(x_limits_moving))-x_limits_moving(1);
y_proj_moving=single(sum(moving_fiducial,2)~=0);
y_proj_moving=single(max(moving_fiducial,[],2));
y_proj_der_moving=abs(y_proj_moving(1:length(y_proj_moving)-1)-y_proj_moving(2:length(y_proj_moving)));
y_limits_moving=find(y_proj_der_moving);
y_range_moving=y_limits_moving(length(y_limits_moving))-y_limits_moving(1);

% PRE-SCALE THE FIXED AND THE MOVING FIDUCIAL IMAGES BASED ON THE LSCM AND MALDI IMAGE SIZE 
[H_fix,W_fix]=size(MALDI_fiducial_binary);
[H_mov,W_mov]=size(LSCM_fiducial_binary);
pre_resize_scale=min(H_mov/H_fix,W_mov/W_fix);
fixed_fiducial=MALDI_fiducial_binary;
fixed_fiducial_pre_resized=imresize(fixed_fiducial,pre_resize_scale,'nearest');

% FINE-SCALE THE FIXED AND THE MOVING FIDUCIAL IMAGES BASED ON THE LSCM AND MALDI IMAGE SIZE 
% put a frame of black pixels around the fixed fiducial image to guarantee
% ... that the white fiducial marks do not touch the border
[H_fix_pre,W_fix_pre]=size(fixed_fiducial_pre_resized);
fixed_fiducial_pre_resized(1,:)=0;
fixed_fiducial_pre_resized(H_fix_pre,:)=0;
fixed_fiducial_pre_resized(:,1)=0;
fixed_fiducial_pre_resized(:,W_fix_pre)=0;
% 
% find in the MALDI fiducial image the x-range and the y-range beyond
% ... which there are no white pixels
x_proj_fixed=sum(fixed_fiducial_pre_resized,1)~=0;
x_proj_der_fixed=abs(x_proj_fixed(1:length(x_proj_fixed)-1)-x_proj_fixed(2:length(x_proj_fixed)));
x_limits_fixed=find(x_proj_der_fixed);
x_range_fixed=x_limits_fixed(length(x_limits_fixed))-x_limits_fixed(1);
y_proj_fixed=sum(fixed_fiducial_pre_resized,2)~=0;
y_proj_der_fixed=abs(y_proj_fixed(1:length(y_proj_fixed)-1)-y_proj_fixed(2:length(y_proj_fixed)));
y_limits_fixed=find(y_proj_der_fixed);
y_range_fixed=y_limits_fixed(length(y_limits_fixed))-y_limits_fixed(1);
% the resize scale is the ratio of the moving and the fixed range
resize_scale=min(y_range_moving/y_range_fixed,x_range_moving/x_range_fixed);
fixed_fiducial_resized=imresize(fixed_fiducial_pre_resized,resize_scale,'nearest');

% find the centroid of the fixed fiducial image
[r_cntrd_fix, c_cntrd_fix]=find_image_centroid(fixed_fiducial_resized);
[H_fix,W_fix]=size(fixed_fiducial_resized);

% PAD BLACK PIXELS AROUND THE LSCM FIDUCIAL AND THE RESIZED MALDI FIDUCIAL
% ...TO MAKE THEIR PIXEL SIZES MATCH
% calculate the width and the height of both the MALDI and the LSCM padded
% ...image
W_pad=2*max([c_cntrd_mov W_mov-c_cntrd_mov c_cntrd_fix W_fix-c_cntrd_fix]);
H_pad=2*max([r_cntrd_mov H_mov-r_cntrd_mov r_cntrd_fix H_fix-r_cntrd_fix]);

% pad the MALDI fiducial
left_pad_fix=W_pad/2-c_cntrd_fix;% pixel columns to be padded on the left
right_pad_fix=W_pad-W_fix-left_pad_fix;% pixel columns to be padded on the right
upper_pad_fix=H_pad/2-r_cntrd_fix;% pixel rows to be padded above the top
lower_pad_fix=H_pad-H_fix-upper_pad_fix;% pixel rows to be padded below the bottom
[fixed_fiducial_padded,r_cntrd_pad_fix,c_cntrd_pad_fix]...
    =pad_around(fixed_fiducial_resized,left_pad_fix,right_pad_fix,upper_pad_fix,lower_pad_fix,r_cntrd_fix,c_cntrd_fix);

% pad the LSCM fiducial
left_pad_mov=W_pad/2-c_cntrd_mov;% pixel columns to be padded on the left
right_pad_mov=W_pad-W_mov-left_pad_mov;% pixel columns to be padded on the right
upper_pad_mov=H_pad/2-r_cntrd_mov;% pixel rows to be padded above the top
lower_pad_mov=H_pad-H_mov-upper_pad_mov;% pixel rows to be padded below the bottom
[moving_fiducial_padded,r_cntrd_pad_mov,c_cntrd_pad_mov]...
    =pad_around(moving_fiducial,left_pad_mov,right_pad_mov,upper_pad_mov,lower_pad_mov,r_cntrd_mov,c_cntrd_mov);
[moving_mask_padded,r_cntrd_pad_mov,c_cntrd_pad_mov]...
    =pad_around(moving_mask,left_pad_mov,right_pad_mov,upper_pad_mov,lower_pad_mov,r_cntrd_mov,c_cntrd_mov);

figure('Name','fixed_moving_fiducial_padded','NumberTitle','off');
% show a pseudo-color overlay of the preregistered MALDI and LSCM fiducial markers 
imshowpair(fixed_fiducial_padded, moving_fiducial_padded,'falsecolor');
hold on
% put a cross at the common centroid of the 2 preregistered images
plot(c_cntrd_pad_fix,r_cntrd_pad_fix,'-gx')
hold off
close all hidden

% REGISTER LSCM FIDUCIALS ONTO MALDI FIDUCIALS
transformType='similarity';
[optimizer,metric]  = imregconfig('multimodal'); 
optimizer.InitialRadius=optimizer.InitialRadius/10;% optimizer parameters may need some tweaking
tic
% BRING THE FIDUCIALS IN REGISTER 
tform = imregtform(moving_fiducial_padded,fixed_fiducial_padded,transformType,optimizer,metric);
% tform will be used to align also any other LSCM images of the same
% ... spheroid
ElapsedTime=toc;
disp(['ElapsedTime=' num2str(ElapsedTime)])
cd ..%root
[mkdir_status,mkdir_message] = mkdir('4_Results');
cd('4_Results')

% IF THE RED LSCM CHANNEL IS PRESENT, REGISTER IT
if ~isempty (LSCM_R_name);
    [moving_red_padded,r_cntrd_pad_mov,c_cntrd_pad_mov]...
        =pad_around(moving_red,left_pad_mov,right_pad_mov,upper_pad_mov,lower_pad_mov,r_cntrd_mov,c_cntrd_mov);
    moving_red_registered = imwarp(moving_red_padded,tform,'OutputView',imref2d(size(fixed_fiducial_padded)));% uses tform found by fiducial registration
    imwrite(uint16(moving_red_registered),'moving_red_registered.tif','Compression','lzw')
    end
% IF THE GREEN LSCM CHANNEL IS PRESENT, REGISTER IT
if ~isempty (LSCM_G_name);
    [moving_green_padded,r_cntrd_pad_mov,c_cntrd_pad_mov]...
        =pad_around(moving_green,left_pad_mov,right_pad_mov,upper_pad_mov,lower_pad_mov,r_cntrd_mov,c_cntrd_mov);
    moving_green_registered = imwarp(moving_green_padded,tform,'OutputView',imref2d(size(fixed_fiducial_padded)));% uses tform found by fiducial registration
    imwrite(uint16(moving_green_registered),'moving_green_registered.tif','Compression','lzw')
end
% IF THE BLUE LSCM CHANNEL IS PRESENT, REGISTER IT
if ~isempty (LSCM_B_name);
    [moving_blue_padded,r_cntrd_pad_mov,c_cntrd_pad_mov]...
        =pad_around(moving_blue,left_pad_mov,right_pad_mov,upper_pad_mov,lower_pad_mov,r_cntrd_mov,c_cntrd_mov);
    moving_blue_registered = imwarp(moving_blue_padded,tform,'OutputView',imref2d(size(fixed_fiducial_padded)));% uses tform found by fiducial registration
%     imwrite_tiff_single_LZW(single(moving_blue_registered),'moving_blue_registered.tif')
    imwrite(uint16(moving_blue_registered),'moving_blue_registered.tif','Compression','lzw')
end
imwrite(uint16(fixed_fiducial_padded),'fixed_fiducial_padded.tif','Compression','lzw')
imwrite(uint16(moving_fiducial_padded),'moving_fiducial_padded.tif','Compression','lzw')

% REGISTER THE BINARY LSCM FIDUCIAL
moving_fiducial_registered = imwarp(moving_fiducial_padded,tform,'OutputView',imref2d(size(fixed_fiducial_padded)));
imwrite(uint16(moving_fiducial_registered),'moving_fiducial_registered.tif','Compression','lzw')

% REGISTER THE BINARY SPHEROID MASK
moving_mask_registered = imwarp(moving_mask_padded,tform,'OutputView',imref2d(size(fixed_fiducial_padded)));% uses tform found by fiducial registration
imwrite(uint16(moving_mask_registered),'moving_mask_registered.tif','Compression','lzw')

% save a pseudo-color overlay of the registered MALDI and LSCM fiducial markers 
fixed_moving=imfuse(fixed_fiducial_padded, moving_fiducial_registered,'falsecolor');
imwrite(fixed_moving,'fixed_moving.tif','Compression','lzw');
close all
% show a pseudo-color overlay of the registered MALDI and LSCM fiducial markers 
figure('Name','fixed_moving_fiducial_registered','NumberTitle','off');
imshowpair(fixed_fiducial_padded, moving_fiducial_registered,'falsecolor');

MALDI_image_gray=MALDI_perifosine_float;
% RESIZE MALDI PERIFOSINE IMAGE
MALDI_image_gray=imresize(MALDI_image_gray,pre_resize_scale,'nearest');
MALDI_image_gray=imresize(MALDI_image_gray,resize_scale,'nearest');
% PAD MALDI PERIFOSINE IMAGE JUST LIKE WE PADDED THE FIDUCIALS
[MALDI_image_padded,r_cntrd_pad_fix,c_cntrd_pad_fix]...
    =pad_around(MALDI_image_gray,left_pad_fix,right_pad_fix,upper_pad_fix,lower_pad_fix,r_cntrd_fix,c_cntrd_fix);
% erase possible fiducial shadows from the perifosine image
MALDI_image=MALDI_image_padded.*~fixed_fiducial_padded;
[mkdir_status,mkdir_message] = mkdir('weighted_means');
cd 'weighted_means'
[mkdir_status,mkdir_message] = mkdir(MALDI_perifosine_name);% MALDI perifosine window may slightly vary and we respect it in the directory name
cd(MALDI_perifosine_name);
imwrite(uint16(MALDI_image),[MALDI_perifosine_name '_padded.tif'],'Compression','lzw')

% GENERATE THE PEELS
% compute the distance map
method='euclidean';
I_binary=single(~moving_mask_registered);
[DistMap_plus]=bwdist(I_binary,method);% for the inward direction starting from the spheroid boundary
I_binary=single(moving_mask_registered);
[DistMap_minus]=bwdist(I_binary,method);% for the outward direction starting from the spheroid boundary
DistMap=DistMap_plus-DistMap_minus;
min_dist=min(min(DistMap));
max_dist=max(max(DistMap));

% FOR EACH PEEL CALCULATE THE TO-PRO-WEIGHTED MEAN VALUE, RELATIVE STANDARD
% DIFFERENCE, AND RED, GREEN OR BLUE MINIMA OR MAXIMA
step=1;
dist_vector=-10:max_dist;% use also peels just slightly outside the boundary
% allocate memory for peel-specific quantities
MALDI_mean=zeros(size(dist_vector),'single');
red_weighted_mean=zeros(size(dist_vector),'single');
green_weighted_mean=zeros(size(dist_vector),'single');
blue_mean=zeros(size(dist_vector),'single');
red_RSD=zeros(size(dist_vector),'single');
green_RSD=zeros(size(dist_vector),'single');
blue_RSD=zeros(size(dist_vector),'single');
MALDI_and_blue_min=zeros(size(dist_vector),'single');
MALDI_and_blue_max=zeros(size(dist_vector),'single');
MALDI_and_red_min=zeros(size(dist_vector),'single');
MALDI_and_red_max=zeros(size(dist_vector),'single');
MALDI_and_green_min=zeros(size(dist_vector),'single');
MALDI_and_green_max=zeros(size(dist_vector),'single');
% evaluate peel-dependent quantities at peels defined by the distance vector 
for index=1:length(dist_vector)
    dist=dist_vector(index);% distance of the current peel from the boundary
    peel=single((DistMap>=dist)&(DistMap<(dist+step)));
    lin_ind_peel=find(peel);%finds linear indices of pixels that belong to the peel
    weight=ones(size(lin_ind_peel));% default, will be replaced by TO-PRO weights if they are present
    sum_of_weights=max(sum(weight),realmin('single'));% default, will be replaced by TO-PRO weights if they are present
    MALDI_peel=MALDI_image(lin_ind_peel);% perifosine values on the peel's pixels
    MALDI_mean(index)=mean(MALDI_peel,'omitnan');% mean perifosine value on the peel
    MALDI_RSD(index)= 100*std(MALDI_peel,'omitnan')/(mean(MALDI_peel,'omitnan'));% perifosine RSD on the peel
   
    if ~isempty (LSCM_B_name)% TO-PRO channel available
        blue_peel=moving_blue_registered(lin_ind_peel);% TO-PRO values on the peel's pixels
        blue_mean(index)=mean(blue_peel,'omitnan');% mean TO-PRO value on the peel
        blue_RSD(index)= 100*std(blue_peel,'omitnan')/(mean(blue_peel,'omitnan'));% TO-PRO RSD on the peel
%         [r,c,MALDI_and_blue]=find(MALDI_image.*double(moving_blue_registered>0).*peel);
        ind=find(blue_peel>0);% find perifosine ranges on pixels with non-zero TO-PRO
        MALDI_and_blue=MALDI_peel(ind);
        if ~isempty(MALDI_and_blue)
            MALDI_and_blue_min(index)=min(MALDI_and_blue);
            MALDI_and_blue_max(index)=max(MALDI_and_blue);
        end
        weight=blue_peel;% for weighted averaging using TO-PRO
        sum_of_weights=max(sum(weight),realmin('single'));
        MALDI_weighted_mean=sum(sum(weight.*MALDI_peel))/sum_of_weights;
    end
    if ~isempty (LSCM_R_name)
        red_peel=moving_red_registered(lin_ind_peel);
        red_weighted_mean(index)=sum(weight.*red_peel)/sum_of_weights;
        red_RSD(index)= 100*std(red_peel,'omitnan')/(mean(red_peel,'omitnan'));
%         [r,c,MALDI_and_red]=find(MALDI_image.*double(moving_red_registered>0).*peel);
        ind=find(red_peel>0);% find perifosine ranges on pixels with non-zero Ki67
        MALDI_and_red=MALDI_peel(ind);
        if ~isempty(MALDI_and_red)
            MALDI_and_red_min(index)=min(MALDI_and_red);
            MALDI_and_red_max(index)=max(MALDI_and_red);
        end
    end
    if ~isempty (LSCM_G_name)
        green_peel=moving_green_registered(lin_ind_peel);
        green_RSD(index)= 100*std(green_peel,'omitnan')/(mean(green_peel,'omitnan'));
        green_weighted_mean(index)=sum(sum(weight.*green_peel))/sum_of_weights;
%         [r,c,MALDI_and_green]=find(MALDI_image.*double(moving_green_registered>0).*peel);
        ind=find(green_peel>0);% find perifosine ranges on pixels with non-zero green
        MALDI_and_green=MALDI_peel(ind);
        if ~isempty(MALDI_and_green)
            MALDI_and_green_min(index)=min(MALDI_and_green);
            MALDI_and_green_max(index)=max(MALDI_and_green);
        end
    end
end

% PLOT PROFILES OF VARIOUS STATISTICAL QUATITIES OVER THE PEELS
figure('Name','Mean values along peels','NumberTitle','off',...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
row=1;
subplot(N_graphs,1,1); plot(dist_vector,MALDI_mean,'m.');title(MALDI_perifosine_name, 'Interpreter', 'none');
if ~isempty (LSCM_R_name)
    row=row+1;
    subplot(N_graphs,1,row); plot(dist_vector,red_weighted_mean,'r.');title('LSCM\_weighted\_red');
    csvwrite('LSCM_weighted_red_peels.csv',red_weighted_mean);
end
if ~isempty (LSCM_G_name)
    row=row+1;
    subplot(N_graphs,1,row); plot(dist_vector,green_weighted_mean,'g.');title('LSCM\_weighted\_green');
    csvwrite('LSCM_weighted_green_peels.csv',green_weighted_mean);
end
if ~isempty (LSCM_B_name)
    row=row+1;
    subplot(N_graphs,1,row); plot(dist_vector,blue_mean,'b.');title('LSCM\_blue');
    csvwrite('LSCM_blue_peels.csv',blue_mean);
end
% plotted profiles are saved to pdf as well as to comma separated values
saveas(gcf,'Weighted_mean_values_along_peels.pdf')
csvwrite( 'dist_vector_peels.csv',dist_vector);
csvwrite('MALDI_slice_peels.csv',MALDI_mean);

N_RSD_graphs=N_graphs;
figure('Name','RSD charts','NumberTitle','off',...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
row=1;
subplot(N_RSD_graphs,1,row); plot(dist_vector,MALDI_RSD,'m.');
title([MALDI_perifosine_name '_RSD[%]'], 'Interpreter', 'none');
csvwrite([MALDI_perifosine_name '_RSD.csv'],MALDI_RSD);
if ~isempty (LSCM_R_name)
    row=row+1;
    subplot(N_RSD_graphs,1,row); plot(dist_vector,red_RSD,'r.');
    title('red\_RSD[%]');
    csvwrite('red_RSD.csv',red_RSD);
end
if ~isempty (LSCM_G_name)
    row=row+1;
    subplot(N_RSD_graphs,1,row); plot(dist_vector,green_RSD,'g.');
    title('green\_RSD[%]');
    csvwrite('green_RSD.csv',green_RSD);
end
if ~isempty (LSCM_B_name)
    row=row+1;
    subplot(N_RSD_graphs,1,row); plot(dist_vector,blue_RSD,'b.');
    title('blue\_RSD[%]');
    csvwrite('blue_RSD.csv',blue_RSD);
end
saveas(gcf,'RSD_along_peels.pdf')

figure('Name','MALDI min max on non-zero LSCM','NumberTitle','off',...
    'Position',[5 105 .5*screen_width-5 screen_height-190])%[left,bottom,width,height];
row=1;

subplot(N_graphs,1,row); plot(dist_vector,MALDI_mean,'m.');title(['MALDI average ' MALDI_perifosine_name], 'Interpreter', 'none');
if ~isempty (LSCM_R_name)
    row=row+1;
    subplot(N_graphs,1,row); plot(dist_vector,MALDI_and_red_min,'r.',...
        dist_vector,MALDI_and_red_max,'r.')
    title('Bounds\_of\_MALDI\_and\_red');
    csvwrite('MALDI_and_red_min.csv',MALDI_and_red_min);
    csvwrite('MALDI_and_red_max.csv',MALDI_and_red_max);
end
if ~isempty (LSCM_G_name)
    row=row+1;
    subplot(N_graphs,1,row); plot(dist_vector,MALDI_and_green_min,'g.',...
        dist_vector,MALDI_and_green_max,'g.')
    title('Bounds\_of\_MALDI\_and\_green');
    csvwrite('MALDI_and_green_min.csv',MALDI_and_green_min);
    csvwrite('MALDI_and_green_max.csv',MALDI_and_green_max);
end
if ~isempty (LSCM_B_name)
    row=row+1;
    subplot(N_graphs,1,row); plot(dist_vector,MALDI_and_blue_min,'b.',...
        dist_vector,MALDI_and_blue_max,'b.')
    title('Bounds\_of\_MALDI\_and\_blue');
    csvwrite('MALDI_and_blue_min.csv',MALDI_and_blue_min);
    csvwrite('MALDI_and_blue_max.csv',MALDI_and_blue_max);
end
box on

saveas(gcf,'MALDI_minmax_for_nonzero_LSCM.pdf')

if ~isempty (LSCM_R_name)
    red_MALDI=imfuse(moving_red_registered, MALDI_image,'falsecolor');
    figure('Name','red_MALDI','NumberTitle','off');imshow(red_MALDI,[])
    imwrite(red_MALDI,'red_MALDI.tif','Compression','lzw')
end
if ~isempty (LSCM_G_name)
    green_MALDI=imfuse(moving_green_registered, MALDI_image,'falsecolor');
    figure('Name','green_MALDI','NumberTitle','off');imshow(green_MALDI,[])
    imwrite(green_MALDI,'green_MALDI.tif','Compression','lzw')
end
if ~isempty (LSCM_B_name)
    blue_MALDI=imfuse(moving_blue_registered, MALDI_image,'falsecolor');
    figure('Name','blue_MALDI','NumberTitle','off');imshow(blue_MALDI,[])
    imwrite(blue_MALDI,'blue_MALDI.tif','Compression','lzw')
end

cd ..%'weighted_means'
cd ..%4_Results
addpath(genpath(cd))
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------
