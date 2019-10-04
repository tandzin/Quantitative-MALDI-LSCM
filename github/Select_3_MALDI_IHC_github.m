LSCM_R_name=[];
LSCM_B_name=[];
LSCM_G_name=[];
LSCM_T_name=[];
LSCM_mask_name=[];
MALDI_perifosine_name=[];
MALDI_fiducial_name=[];


image_nr=str2num(input('image_nr [1...3]:','s'));

if image_nr==1 % XXXXX
    cd('s1')
    %  dataset 28
    LSCM_R_name='Spheroid_red.tif';
    LSCM_B_name='Spheroid_blue.tif';
    LSCM_T_name='MAX_Spheroid_transm.tif';
    LSCM_mask_name='MAX_Spheroid_blue_mask.tif';
    se = strel('disk',15);
    LSCM_fiducial_threshold=50000;
    MALDI_perifosine_name='462.3-462.7@5524.tif';
    MALDI_fiducial_name='553.2-553.6_fiducial.png';
    MALDI_fiducial_threshold=0.09;
elseif image_nr==2 % XXXXX
    cd('s2')
    % dataset 31
    LSCM_R_name='Spheroid_red.tif';
    LSCM_B_name='Spheroid_blue.tif';
    LSCM_T_name='MAX_Spheroid_transm.tif';
    LSCM_mask_name='MAX_Spheroid_blue_mask.tif';
    se = strel('disk',15);
    LSCM_fiducial_threshold=50000;
    MALDI_perifosine_name='462.4-462.8@3111.tif';
    MALDI_fiducial_name='553.2-553.6_fiducial.png';
    MALDI_fiducial_threshold=0.09;
    
elseif image_nr==3 % XXXXX
    cd('s3')
    %  dataset 9
    LSCM_R_name='Spheroid_red.tif';
    LSCM_B_name='Spheroid_blue.tif';
    LSCM_T_name='patchedAVG_Spheroid_transm.tif';
    LSCM_mask_name='MAX_Spheroid_blue_mask.tif';
    se = strel('disk',15);
    LSCM_fiducial_threshold=50000;
    MALDI_perifosine_name='462.2-462.6@3246.tif';
    MALDI_fiducial_name='553.2-553.6_fiducial.png';
    MALDI_fiducial_threshold=0.09;
    
else
    return
end

