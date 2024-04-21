%% 狗鼻影像切割
%% owner: yang shu chun
clear all;

dir='base';
path=['..\data\' dir '\'];       % 資料位置
result_p=['..\data\' dir '\'];   % 基準樣本/測試樣本輸出位置

img_path=[path 'image\resize\'];    % 原始圖片位置
recall_path=[path 'recall\'];       % YOLO框選左右鼻孔位置
txt_path=[result_p 'txt\'];         % 基準線輸出位置
mask_path=[result_p 'mask\'];       % 狗鼻遮罩輸出影像位置
output_path=[result_p 'output\'];   % 狗鼻紋輸出位置

%% 狗鼻紋切割參數
EH_r=0.9;               % 手動調整
LG_Gamma=4;             % 學姊論文p.58 [r]
LG_Wsize=60;            % 學姊論文p.58 [wg]
Min_Wsize=1;            % 3*3
RL_Wsize=10;            % 學姊論文p.58 [wl]
MS_Wsize=40;            % 學姊論文p.58 [wm]
MS_step=2;              % 學姊論文p.58 [step]
connect_TH=11;          % 學姊論文p.58 [TH]
connect_Maxlen=11;      % 學姊論文p.58 [max_len]
% lower=-0.4;           % 學姊論文p.54 [k_low]
lower=-0.6;             % 手動調整參數
upper=1.45;             % 學姊論文p.54 [k_high]

%%

files = dir([img_path '*.bmp']);
for i=1:length(files)
% parfor i=1:length(files)
%     tic 

    name_part_num = files(i).name(1:(end-4));
    
%   讀取狗鼻原始影像 ==========================================
    img_name=[img_path files(i).name];
    
    if (exist(img_name,'file'))
        img = imread(img_name);
    else
        continue;
    end
    
%   讀取YOLO左右鼻孔定位位置 ==================================
    left_txt_name = [recall_path 'left\left_cutimg_' name_part_num '.txt'];
    right_txt_name = [recall_path 'right\right_cutimg_' name_part_num '.txt'];

    if (exist(left_txt_name,'file') && exist(right_txt_name,'file'))
        txt1 = fopen(left_txt_name, 'r');
        txt2 = fopen(right_txt_name, 'r');
    else
        continue;
    end
    
    p1=fscanf(txt1, '%d');
    p2=fscanf(txt2, '%d');      

    if (~isempty(p1) && ~isempty(p2))
        [y, x, z]=size(img);
        [left_position]=image_position (0, 0, p1);          %取得Yolo鼻孔定位之座標
        [right_position]=image_position ((x/2), 0, p2);
    else
        continue;
    end

%   繪畫狗鼻遮罩、旋轉前之原始基準線、旋轉後之水平基準線 ======    
    [mask, L_n, new_L_n]=image_label (txt_path, mask_path, name_part_num, img, left_position, right_position, p1(3:4), p2(3:4));

%   切割狗鼻紋影像 ==============================================    
    gray=Gray_img(img);
    enhance=Enhance_img(gray, EH_r);
    histogram=Histogram_img(enhance);
    localgamma=Local_gamma (histogram, LG_Gamma, LG_Wsize);
    median=Median_filter (localgamma);
    minimum=Minimum (median, Min_Wsize);
    runlength=Run_Length (median, RL_Wsize);
    avg=Avg (minimum, runlength);
    multiscale=Multiscale_line (avg, mask, MS_Wsize, MS_step);
        
%   旋轉圖像 ==============================================    
    [rotation_multiscale]=image_rotate_inverse (multiscale, L_n, 1, 'uint8');  % 背景為為黑色(0) 白色(1 or 255)

    thredshold=Thredshold (rotation_multiscale);    
    fillhole=Fillhole (thredshold);
    closing=Closing(fillhole);
    thin=Thining (closing);
    connect=Connect (thin, connect_TH, connect_Maxlen);
    spur=Spur(connect);
    [roiarea, whitearea, result]=Adjust_seg(spur, new_L_n, lower, upper);

%   輸出鼻紋影像    
    if ~exist(output_path, 'dir')
        % Folder does not exist so create it.
        mkdir(output_path);
    end
    
    img_code=[name_part_num  '_EHr=' num2str(EH_r)  '_LG=' num2str(LG_Gamma) '_LW=' num2str(LG_Wsize) '_MW=' num2str(Min_Wsize) '_RW=' num2str(RL_Wsize)  '_MsW=' num2str(MS_Wsize) '_Mstep=' num2str(MS_step)  '_TH=' num2str(connect_TH) '_Maxlen=' num2str(connect_Maxlen)  '_lower=' num2str(lower) '_upper=' num2str(upper) '.bmp'];    

%     output_img(mask, '0. mask', output_path, img_code);
%     output_img(gray, '1. gray', output_path, img_code);
%     output_img(uint8(enhance), '2. enhance', output_path, img_code);
%     output_img(histogram, '3. histogram', output_path, img_code);
%     output_img(localgamma, '4. localgamma', output_path, img_code);
%     output_img(median, '5. median', output_path, img_code);
%     output_img(minimum, '6.1 minimum', output_path, img_code);
%     output_img(runlength, '6.2 runlength', output_path, img_code);
%     output_img(avg, '7. avg', output_path, img_code);
%     output_img(multiscale, '8.1 multiscale', output_path, img_code);
%     output_img(rotation_multiscale, '8.2 rotation_multiscale', output_path, img_code);
%     output_img(thredshold, '9. thredshold', output_path, img_code);
%     output_img(fillhole, '10. fillhole', output_path, img_code);
%     output_img(closing, '11. closing', output_path, img_code);
%     output_img(thin, '12. thin', output_path, img_code);
%     output_img(connect, '13. connect', output_path, img_code);
%     output_img(spur, '14. spur', output_path, img_code);
    output_img(result, '15. result', output_path, img_code);
    
    fclose('all');
    
%     toc
end

% 輸出指定影像
function output_img(img, folder_name, path, img_code)
    output_path=[path folder_name '\'];
    if ~exist(output_path, 'dir')
        % Folder does not exist so create it.
        mkdir(output_path);
    end
    
    img_name= [output_path folder_name '_' img_code];
    imwrite(img, img_name);
end

% yolo 框取之輪廓
function [position]=image_position (start_x, start_y, p)
        x_0=start_x+p(1);
        y_0=start_y+p(2);
        
        x_1=x_0+p(3);
        y_1=y_0+p(4);
        
        position=[x_0, x_1; y_0, y_1];        
end

%% 建立基準線與mask (旋轉前與旋轉後)
function [mask, L_n, new_L_n]=image_label (txt_path, mask_path, name_part_num, img, left, right, left_n, right_n)
    
    left_c=[(left(1, 1)+left(1, 2))*0.5; (left(2, 1)+left(2, 2))*0.5];
    right_c=[(right(1, 1)+right(1, 2))*0.5; (right(2, 1)+right(2, 2))*0.5];

    v1=[left(1, 2), left_c(2)]-[left_c(1), left_c(2)];
    v2=[right_c(1), right_c(2)]-[left_c(1), left_c(2)];
    theta=acosd(sum(v1.*v2)/(norm(v1)*norm(v2))); % 旋轉角度
    if theta>30
        % 當旋轉角度過大時，鼻孔長寬比例會產生偏差，故取長寬相加除以二
        new_left=[(left_n(1)+left_n(2))*0.5; (left_n(1)+left_n(2))*0.5];
        new_right=[(right_n(1)+right_n(2))*0.5; (right_n(1)+right_n(2))*0.5];
    else
        new_left=left_n;
        new_right=right_n;
    end
    
    L_c=[double((left_c(1)*new_right(1)+right_c(1)*new_left(1))/(new_left(1)+new_right(1))); double((left_c(2)*new_right(1)+right_c(2)*new_left(1))/(new_left(1)+new_right(1)));];
    L_n=[left_c, L_c, right_c];
          
    L_d=points_distance(L_n(:, 1), L_n(:, 2));
    R_d=points_distance(L_n(:, 3), L_n(:, 2));
    new_L_n=[(L_n(1, 2)-L_d), L_n(1, 2), (L_n(1, 2)+R_d); L_n(2, 2), L_n(2, 2), L_n(2, 2)];
    inverse_L_n=[L_n(1, :); L_n(2, 3), L_n(2, 2), L_n(2, 1)];
    
    gray_img = rgb2gray(img);
    rotated_img=image_rotate_inverse (gray_img, L_n, 0, 'uint8' );
                    
    new_txt=[txt_path 'rotated\'];

    if ~exist(new_txt, 'dir')
        % Folder does not exist so create it.
        mkdir(new_txt);
    end

    dlmwrite([new_txt 'rotated_position_' name_part_num '.txt'], new_L_n);

    [mask]=make_mask (mask_path, name_part_num, rotated_img, new_L_n, inverse_L_n, new_left, new_right);
    
end

% 狗鼻紋遮罩mask
function [mask]=make_mask (path, name_part_num, img, L_n, inverse_L_n, left_n, right_n)

    x_0=L_n(1, 1)-left_n(1)*0.5;
    y_0=L_n(2, 1)-left_n(2)*0.5;

    x_1=L_n(1, 1)+left_n(1)*0.5;
    y_1=L_n(2, 1)+left_n(2)*0.5;

    left=[x_0, x_1, x_1, x_0, x_0; y_0, y_0, y_1 y_1, y_0];

    x_3=L_n(1, 3)-right_n(1)*0.5;
    y_3=L_n(2, 3)-right_n(2)*0.5;

    x_4=L_n(1, 3)+right_n(1)*0.5;
    y_4=L_n(2, 3)+right_n(2)*0.5;

    right=[x_3, x_4, x_4, x_3, x_3; y_3, y_3, y_4 y_4, y_3];

    
    % 繪製左右鼻孔標記影像==============================================
    
% %     new=[ path, 'labeled_img\'];
% % 
% %     if ~exist(new, 'dir')
% %         % Folder does not exist so create it.
% %         mkdir(new);
% %     end
% % 
% %     figure(1);
% %     imshow(img, 'Border','tight');
% %     hold on
% %     plot(left(1, :), left(2, :), '-o','color','r','LineWidth',2);
% %     plot(right(1, :), right(2, :), '-o','color','g','LineWidth',2);
% %     plot(L_n(1, :), L_n(2, :), '-o','color','y','LineWidth',2);
% %     hold off
% % 
% %     f=getframe(gcf);
% %     imwrite(f.cdata,[new 'position_' num2str(name_part_num) '.bmp']); 
% % 
% %     close all;
    
        
    % ============================================================

    % 繪製旋轉後的mask影像==============================================
    
    [row, column]=size(img);    
    mask_tmp=zeros(row, column);
%     mask_view=img;
    
    radius_x_1=L_n(1, 2)*0.95;
    radius_x_2=(column-L_n(1, 2))*0.95;
    radius_y_1=L_n(2, 2)*0.95;
    radius_y_2=(row-L_n(2, 2))*0.95;
    
    % 鼻子輪廓(左上角)
    for i=1:L_n(1, 2)
        for j=1:L_n(2, 2)
            if ((((i-L_n(1, 2))^2)/(radius_x_1^2)) + (((j-L_n(2, 2))^2)/(radius_y_1^2))) < 1
                mask_tmp(j, i)=1;
            end
        end
    end
    % 鼻子輪廓(左下角)
    for i=1:L_n(1, 2)
        for j=round(L_n(2, 2)):row
            if ((((i-L_n(1, 2))^2)/(radius_x_1^2)) + (((j-L_n(2, 2))^2)/(radius_y_2^2))) < 1
                mask_tmp(j, i)=1;
            end
        end
    end
    % 鼻子輪廓(右上角)
    for i=round(L_n(1, 2)):column
        for j=1:L_n(2, 2)
            if ((((i-L_n(1, 2))^2)/(radius_x_2^2)) + (((j-L_n(2, 2))^2)/(radius_y_1^2))) < 1
                mask_tmp(j, i)=1;
            end
        end
    end
    % 鼻子輪廓(右下角)
    for i=round(L_n(1, 2)):column
        for j=round(L_n(2, 2)):row
            if ((((i-L_n(1, 2))^2)/(radius_x_2^2)) + (((j-L_n(2, 2))^2)/(radius_y_2^2))) < 1
                mask_tmp(j, i)=1;
            end
        end
    end
    
    % 左側鼻孔
    left_column=left(1, 3)-left(1, 1);
    left_row=left(2, 3)-left(2, 1);
    left_a = left_column*0.5*0.95 ; 
    left_b = left_row*0.5*0.95 ;
    for i=L_n(1, 1):left(1, 3)
        for j=left(2, 1):left(2, 3)
            if ((((i-L_n(1, 1))^2)/(left_a^2)) + (((j-L_n(2, 1))^2)/(left_b^2))) < 1 && i>=1
                mask_tmp(round(j), round(i))=0;
            end
        end
    end
    for i=left(1, 1):L_n(1, 1)
        for j=L_n(2, 1):left(2, 3)
            if ((((i-L_n(1, 1))^2)/(left_a^2)) + (((j-L_n(2, 1))^2)/(left_b^2))) < 1 && i>=1
                mask_tmp(round(j), round(i))=0;
            end
        end
    end
    left_r1=left_a*0.5;
    left_c1=[(L_n(1, 1)-left_r1), L_n(2, 1)];
    for i=left(1, 1):L_n(1, 1)
        for j=left(2, 1):left(2, 3)
            if ((((i-left_c1(1))^2)/(left_r1^2)) + (((j-left_c1(2))^2)/(left_r1^2))) < 1 && i>=1
                mask_tmp(round(j), round(i))=1;
            end
        end
    end
    left_r2=[left_b*0.3, left_b*0.5];
    left_c2=[L_n(1, 1), (L_n(2, 1)-left_r2(2))];
    for i=left(1, 1):left(1, 3)
        for j=left(2, 1):L_n(2, 1)
            if ((((i-left_c2(1))^2)/(left_r2(1)^2)) + (((j-left_c2(2))^2)/(left_r2(2)^2))) < 1 && i>=1
                mask_tmp(round(j), round(i))=0;
            end
        end
    end
    
    % 右側鼻孔
    right_column=right(1, 3)-right(1, 1);
    right_row=right(2, 3)-right(2, 1);
    right_a = right_column*0.5*0.95 ; 
    right_b = right_row*0.5*0.95 ;
    
    for i=right(1, 1):L_n(1, 3)
        for j=right(2, 1):right(2, 3)
            if ((((i-L_n(1, 3))^2)/(right_a^2)) + (((j-L_n(2, 3))^2)/(right_b^2))) < 1 && i<=column
                mask_tmp(round(j), round(i))=0;
            end
        end
    end
    for i=L_n(1, 3):right(1, 3)
        for j=L_n(2, 3):right(2, 3)
            if ((((i-L_n(1, 3))^2)/(right_a^2)) + (((j-L_n(2, 3))^2)/(right_b^2))) < 1 && i<=column
                mask_tmp(round(j), round(i))=0;
            end
        end
    end
    right_r1=right_a*0.5;
    right_c1=[(L_n(1, 3)+right_r1), L_n(2, 3)];
    for i=L_n(1, 3):right(1, 3)
        for j=right(2, 1):right(2, 3)
            if ((((i-right_c1(1))^2)/(right_r1^2)) + (((j-right_c1(2))^2)/(right_r1^2))) < 1 && i<=column
                mask_tmp(round(j), round(i))=1;
            end
        end
    end
    right_r2=[right_b*0.3, right_b*0.5];
    right_c2=[L_n(1, 3), (L_n(2, 3)-right_r2(2))];
    for i=right(1, 1):right(1, 3)
        for j=right(2, 1):L_n(2, 3)
            if ((((i-right_c2(1))^2)/(right_r2(1)^2)) + (((j-right_c2(2))^2)/(right_r2(2)^2))) < 1 && i<=column
                mask_tmp(round(j), round(i))=0;
            end
        end
    end

    rotated_mask=mask_tmp;

    % 旋轉前影像mask
    [mask]=image_rotate_inverse (rotated_mask, inverse_L_n, 0, 'logical');  % 背景為為黑色(0) 白色(1 or 255)

end

% 影像旋轉公式
function [rotation_img]=image_rotate_inverse (img, L_n, back_code, image_type)

    [row, column]=size(img);
    
    % 求基準點離四周的距離最大值 ==>為了把圖片擴張成正方形(正中心為基準點)，以免旋轉時跑到界線外
    expand_dist=max([L_n(1, 2), (column-L_n(1, 2)), L_n(2, 2), (row-L_n(2, 2))]); 

    % 影像長寬起始位置 (用來之後擴增後縮回到原始大小)
    column_start=1; column_end=column;
    row_start=1; row_end=row;
    
    % 四周各需要擴張的距離
    column_dist_a=round(expand_dist-L_n(1, 2));
    column_dist_b=round(expand_dist-(column-L_n(1, 2)));
    row_dist_a=round(expand_dist-L_n(2, 2));
    row_dist_b=round(expand_dist-(row-L_n(2, 2)));
    
    column_start=column_start+column_dist_a;    % 左側擴張==>X軸起始點改變，需要紀錄新的起始點
    column_end=column_end+column_dist_a;      % 左側擴張==>X軸結束點改變，需要紀錄新的結束點
    row_start=row_start+row_dist_a;                     % 上側擴張==>Y軸起始點改變，需要紀錄新的起始點
    row_end=row_end+row_dist_a;                       % 上側擴張==>Y軸結束點改變，需要紀錄新的結束點
        
    temp_1=[];
    temp_2=[];
    temp_3=[];
    temp_4=[];
    
    if column_dist_a>0
        temp_1=ones(row, column_dist_a)*back_code;
    end
    if column_dist_b>0
        temp_2=ones(row, column_dist_b)*back_code;
    end
    
    temp_c=[temp_1 img temp_2];
    temp_r=temp_c;
    
    if row_dist_a>0
        temp_3=ones(row_dist_a, size(temp_c,2))*back_code;
    end
    if row_dist_b>0
        temp_4=ones(row_dist_b, size(temp_c,2))*back_code;
    end
    
    temp_r=[temp_3; temp_r; temp_4];
    
    temp=temp_r;
    
%     figure(12)
%     imshow(temp);
    
    % 旋轉中心點座標
    c=size(temp, 2)/2; % column (x軸)
    r=size(temp, 1)/2; % row (y軸)
    
    v1=[L_n(1, 1), L_n(2, 1)]-[L_n(1, 3), L_n(2, 3)];
    v2=[L_n(1, 1), L_n(2, 2)]-[L_n(1, 3), L_n(2, 2)];
    theta=acosd(sum(v1.*v2)/(norm(v1)*norm(v2))); % 旋轉角度
    if L_n(2, 2)-L_n(2, 3)>0    % 判別是順時鐘(角度為負)旋轉還是逆時鐘(角度為正)
        theta=theta*(-1);
    end

    background=ones(size(temp))*back_code;
    
    % 因為是從新的圖像座標映射到原始圖像中的座標，故旋轉角度為相反
    theta_tmp=(-1)*theta;
    
    % 旋轉矩陣之反矩陣
    inverse_R=[cosd(theta_tmp) sind(theta_tmp); -sind(theta_tmp) cosd(theta_tmp)];
    
%     % R=[cosd(theta) -sind(theta);sind(theta) cosd(theta)];

    [inverse_row, inverse_column]=size(temp);

    for j=1:inverse_row             % y 軸
        for i=1:inverse_column   % x 軸
            a=[i-c; j-r];           % 新圖像之以中心點為(0, 0)的座標位置
            b=inverse_R*a;    % 映射到原始圖像中的位置(中心點為(0, 0))

            b(1)=b(1)+c; % 加上中心點座標==>表示在原始圖像中相應的座標位置(X軸)
            b(2)=b(2)+r; % 加上中心點座標==>表示在原始圖像中相應的座標位置(Y軸)

            if (b(1)>=1 && b(1)<=inverse_column && b(2)>=1 && b(2)<=inverse_row)
                background(j, i)=temp(round(b(2)), round(b(1)));    % 將原始圖像中的值映射到新圖像(座標取整數)
            end
            
        end
    end
    
    if strcmp(image_type, 'uint8')
        b=uint8(background);
    end
    if strcmp(image_type, 'logical')
        b=background;
    end 
    
%     figure(13)
%     imshow(b);

    rotation_img=b(row_start: row_end, column_start: column_end); % 位移到圖像原始大小
   
end

%% 鼻紋提取流程====================================================

% 中心點擴張win*win的範圍座標
function [x_idx, y_idx]=window_idx(i, j, column, row, win)
    low_width = i-win;      up_width = i+win;
    low_hight = j-win;       up_hight = j+win;

    if low_width < 1
        low_width = 1;
    end
    if up_width > column
        up_width = column;
    end
    if low_hight < 1
        low_hight = 1;
    end
    if up_hight > row
        up_hight = row;
    end
    
    x_idx=[low_width, up_width];
    y_idx=[low_hight, up_hight];
    
end

% 求座標p1 p2 之連線距離
function [dist]=points_distance(p1, p2)
    dist=sqrt(((p1(1)-p2(1))^2)+((p1(2)-p2(2))^2));
end

%% 1. gray_img
% 原始影像轉為灰階影像
function [gray_img]=Gray_img(img)

    gray_img = rgb2gray(img); %% 將rgb轉成灰階影像並存到變數img中

end

%% 2. enhance_img
% 增強黑白對比
function [enhance_img]=Enhance_img(img, EH_r)

        [row,lengths] = size(img);
        enhance_img = zeros(row, lengths);
        pixel_value = zeros(row, lengths);  %儲存每個pixel的值

        for i = 1:row
            for j = 1:lengths
                pixel_value(i, j) = img(i, j);
            end
        end

        win=3;

        for i = 1:row
            for j = 1:lengths
                win_temp=zeros(1,win*win);
                count=0;
                for a=fix(-win/2):fix(win/2)
                    for b=fix(-win/2):fix(win/2)
                        if ((i+a)>0 && (i+a)<=row && (j+b)>0 && (j+b)<=lengths)
                            count=count+1;
                            win_temp(1,count)=pixel_value(i+a, j+b);
                        end
                    end
                end

                total=sum(win_temp);
                result=((total/(count*255))^EH_r)*255;

                enhance_img(i, j)=round(result);

            end
        end
     
end

%% 3. histogram_img
% 直方圖等化
function [histogram_img]=Histogram_img(img)

        [row,lengths] = size(img); %%[r,c]=size(A),size函數將矩?的行數返回到第一個輸出變量r，將矩陣的列數返回到第二個輸出變量
        histogram_img = zeros(row,lengths);
        location = zeros(2,row*lengths);  %儲存每個pixel的座標
        pixel_value = zeros(1,row*lengths);  %儲存每個pixel的值
        count = 1;
        for i = 1:row  %將每個點的座標和灰階值存下來
            for j = 1:lengths
                location(1,count) = i;
                location(2,count) = j;
                pixel_value(1,count) = img(i,j);
                count = count + 1;
            end
        end
        [sort_pixel,index] = sort(pixel_value);  %將灰階值由小到大排序,[排序後的陣列,對應回原本陣列的位置]

        tmp_num = double(row*lengths / 256);  %計算出平均每個灰階值需要分配幾個pixel
        gray_level = 1;
        gray_count = 0;
        for a = 1:row*lengths
            if gray_count == 0
                %計算每個灰階值要分配的實際pixel個數(因為tmp_num可能無法整除,如果每個都用四捨五入的個數,會導致到越後面灰階值給的個數誤差越大,所已將每次應該分配到的總個數減去前面已用的個數來表示該灰階值要分配的pixel個數)
                num = round(double(gray_level*tmp_num) - (a-1));
            end
            gray_count = gray_count + 1;

            histogram_img(location(1,index(a)),location(2,index(a))) = gray_level - 1;
            if gray_count == num
                gray_count = 0;
                gray_level = gray_level + 1;
            end
        end
        
        histogram_img=uint8(histogram_img);

end

%% 4. local gamma
% 區域對比拉開，為了凸顯鱗狀區塊在影像中更為突出
function [local_gamma]=Local_gamma (img, LG_Gamma, LG_Wsize)


        [row, column] = size(img); % [Y軸, X軸]
        local_gamma = zeros(size(img)); % 轉換成不帶正負號的8位元的整數

        tmp_r=LG_Gamma; % gamma值
        tmp_R=1/tmp_r; % gamma倒數
        window_size=LG_Wsize;

%             window_size = 60; %window大小=121

         for i = 1:row              % Y軸
            for j = 1:column     % X軸
                [x_idx, y_idx]=window_idx(j, i, column, row, window_size);
                tempBlock=img(y_idx(1):y_idx(2), x_idx(1):x_idx(2));
                
                max_H = max(tempBlock(window_size+1,:));  %垂值最大
                max_V = max(tempBlock(:,window_size+1));  %水平最大
                min_H = min(tempBlock(window_size+1,:));  %垂值最小
                min_V = min(tempBlock(:,window_size+1));  %水平最小
                sum_H = sum(tempBlock(window_size+1,:));  %垂直總合
                sum_V = sum(tempBlock(:,window_size+1));  %水平總合
                length_H = length(tempBlock(window_size+1,:)); %垂直總個數
                length_V = length(tempBlock(:,window_size+1)); %水平總個數
                max_all = max(max_H,max_V);
                min_all = min(min_H,min_V);
                P = double(img(i,j));
                tmp_sum = sum_H+sum_V-P;
                tmp_length = length_H+length_V-1;
                mean_all = round(tmp_sum/tmp_length); %全部平均
                
                if img(i,j) < mean_all
                    local_gamma(i,j) = (((double(img(i,j)-min_all)/double(mean_all-min_all))^tmp_r)*255);
                else
                    local_gamma(i,j) = (((double(img(i,j)-mean_all)/double(max_all-mean_all))^tmp_R)*255);
                end

            end
         end
        
         local_gamma=uint8(local_gamma);

        clear mean_all;

end

%% 5. median_filter
% 去除影像雜訊
function [median_img]=Median_filter (img)
  
        [row,column] = size(img);
        tmp_result = ones(size(img));
        window_size = 1; % 3x3
        
        
         for i = 1:row              % Y軸
            for j = 1:column     % X軸

                [x_idx, y_idx]=window_idx(j, i, column, row, window_size);
                tempBlock=img(y_idx(1):y_idx(2), x_idx(1):x_idx(2));
                
                S = reshape (tempBlock, 1, numel(tempBlock)); % 轉成一維陣列
                
                sorting = sort(S); % 由小到大排序
                mid = median(sorting); % 取中位數
                
                tmp_result(i,j) = mid;
                
                %         end
            end
        end
        
        median_img=uint8(tmp_result);
        
end

%% 6. Minimum Replace Window Center
% 中心被window最小值取代
function [minmum_img]=Minimum (img, Min_Wsize)
       
        [row,lengths] = size(img); %%[r,c]=size(A),size函數將矩?的行數返回到第一個輸出變量r，將矩陣的列數返回到第二個輸出變量
        pixel_value = zeros(row, lengths);  %儲存每個pixel的值
        
        for i = 1:row  %將每個點的座標和灰階值存下來
            for j = 1:lengths
                pixel_value(i, j) = img(i, j);
            end
        end
        
        window=Min_Wsize*2+1;
        window_size=Min_Wsize;
        window_index=zeros(2, window*window);
        window_value=zeros(window*window, row*lengths);
        new_value=zeros(3, row*lengths);
        
        count=1;
        for i = 1:row  %將每個點加上周圍(九個點)pixel值記錄下來
            for j = 1:lengths
                w=1;
                for a=-window_size:window_size
                    for b=-window_size:window_size
                        window_index(1, w)=(i+a); 
                        window_index(2, w)=(j+b);
                        w=w+1;
                    end
                end
                
                temp=zeros(1, window*window);
                for k = 1:window*window
                    if (window_index(1, k)>0 && window_index(1, k)<=row && window_index(2, k)>0 && window_index(2, k)<=lengths)
                        window_value(k, count)=pixel_value(window_index(1, k), window_index(2, k));
                    else
                        window_value(k, count)=-1;
                    end
                    temp(k)=window_value(k, count);
                end
               
                new_value(1, count)=i;
                new_value(2, count)=j;
                new_value(3, count)=min(temp(temp>=0));
                
                count=count+1;
            end
        end
        
        new_img=zeros(row, lengths);
        for x=1:count-1
            new_img(new_value(1, x), new_value(2, x)) = round(new_value(3, x));
        end
        
        minmum_img = uint8(new_img);

end          

%% 6. run length
% 強化影像鱗狀區塊間的輪廓與邊界
function [runlength_img]=Run_Length (img, RL_Wsize)

        smoothGray = zeros(size(img));
        window_size=RL_Wsize;%% 9x9
        
        %%window = 8; %% 17x17
        
        for z = 1:size(img,1)
            for j = 1:size(img,2)
                top = z - window_size;
                if top < 1
                    top = 1;
                end
                down = z + window_size;
                if down > size(img,1)
                    down = size(img,1);
                end
                left = j-window_size;
                if left < 1
                    left = 1;
                end
                right = j+window_size;
                if right > size(img,2)
                    right = size(img,2);
                end
                block = img(top:down, left:right);
                
                xSi = ceil(size(block,1)/2);
                ySi = ceil(size(block,2)/2);
                line8 = double(zeros(8,2));
                step = [1,0,1;2,1,0.5;1,1,1;1,2,0.5;0,1,1;-1,2,0.5;-1,1,1;-2,1,0.5];
                for inneri = 1:8
                    for innerj = 1:2
                        nowx = 0;
                        nowy = 0;
                        while(true)
                            if nowx+xSi<1 || nowx+xSi>size(block,1) || nowy+ySi<1 || nowy+ySi>size(block,2)
                                break;
                            end
                            line8(inneri,1) = line8(inneri,1) + double(block(nowx+xSi,nowy+ySi));
                            line8(inneri,2) = line8(inneri,2) + 1;
                            nowx = nowx + step(inneri,1);
                            nowy = nowy + step(inneri,2);
                        end
                        step(inneri,1:2) = -step(inneri,1:2);
                    end
                end
                
                lineMean = line8(:,1) ./ line8(:,2);
                if img(z,j)>mean(block(:))
                    smoothGray(z,j) = max(lineMean(:));
                elseif img(z,j)<mean(block(:))
                    smoothGray(z,j) = min(lineMean(:));
                else
                    smoothGray(z,j) = img(z,j);
                end
            end
        end
        
        runlength_img = uint8(smoothGray);

end

%% 7. AVG (Minimum+RunLength)
% (Minimum+RunLength)/2
function [avg_img]=Avg (img_minimum, img_runlength)

        [row,lengths] = size(img_runlength); %%[r,c]=size(A),size函數將矩?的行數返回到第一個輸出變量r，將矩陣的列數返回到第二個輸出變量
        pixel_value_1 = zeros(row, lengths);  %儲存每個pixel的值
        pixel_value_2 = zeros(row, lengths);
        new_img=zeros(row, lengths);
        %         new_value=zeros(row, lengths);
        
        for i = 1:row  %將每個點的座標和灰階值存下來
            for j = 1:lengths
                pixel_value_1(i, j) = img_runlength(i, j);
                pixel_value_2(i, j) = img_minimum(i, j);
                new_img(i, j)=round((pixel_value_1(i, j)+pixel_value_2(i, j))/2);
            end
        end

        avg_img = uint8(new_img);
        
end

%% 8. multi-scale line
% 移除不相關的背景像素
function [multiscale_img]=Multiscale_line (img, mask, MS_Wsize, MS_step)
        
    img = im2double(img);
    mask = im2bw(mask); %#ok<IM2BW>

    img = 1-img; % 反白後
    % img = fakepad(img,mask);

    features = standardize(img,mask);
    Ls = 1:MS_step:MS_Wsize; % Ls為初始值為1 5 9 13....W的陣列
    for j = 1:numel(Ls)
        L = Ls(j);
        R = get_lineresponse(img,MS_Wsize,L);
        R = standardize(R,mask);
        features = features+R;
        %    disp(['L = ',num2str(L),' finished!']);
    end
    segmentedimg = features/(1+numel(Ls));
    multiscale_img=uint8(segmentedimg*255);

end

function [R] = get_lineresponse(img,W,L)
    % img: 反白
    % W: window size, L: line length
    % R: line detector response

    avgmask = fspecial('average',W);
    avgresponse = imfilter(img,avgmask,'replicate');

    % w = (W*2)-1;
    % M_avgmask = fspecial('average',w); % 大window之平均
    % M_avgresponse = imfilter(img,M_avgmask,'replicate');
    % 
    % L_avgmask = fspecial('average',L); % 長度L_window之平均
    % L_avgresponse = imfilter(img,L_avgmask,'replicate');

    % gamma = 1.0;
    % weight = ((L_avgresponse./M_avgresponse).^gamma); % 長度L_windowt除以大window，擇一再取gamma
    % weight = (L_avgresponse./M_avgresponse);

    maxlinestrength = -Inf*ones(size(img));
    for theta = 0:15:165
        linemask = get_linemask(theta,L);
        linemask = linemask/sum(linemask(:));
        imglinestrength = imfilter(img,linemask);    
        imglinestrength = imglinestrength - avgresponse;    
        maxlinestrength = max(maxlinestrength,imglinestrength);    
    end
    R = maxlinestrength;
    % R = maxlinestrength.*weight;
    % R = (maxlinestrength.^gamma).*weight;
end

function [linemask] = get_linemask(theta,masksize)
% (theta,masksize)
% Create a mask for line with angle theta
if theta > 90
   mask = getbasemask(180- theta,masksize);
   linemask = rotatex(mask);
else
   linemask = getbasemask(theta,masksize);
end
% imshow(linemask,'InitialMagnification','fit');
end

function [rotatedmask] = rotatex(mask)
[h,w] = size(mask);
rotatedmask = zeros(h,w);

for i = 1:h
    for j = 1:w
        rotatedmask(i,j) = mask(i,w-j+1);
    end
end
end

function [mask] = getbasemask(theta,masksize)

mask = zeros(masksize);

halfsize = (masksize-1)/2;

if theta == 0
    mask(halfsize+1,:) = 1;
elseif theta == 90
    mask(:,halfsize+1) = 1;
else
    x0 = -halfsize;
    y0 = round(x0*(sind(theta)/cosd(theta)));

    if y0 < -halfsize
        y0 = -halfsize;
        x0 = round(y0*(cosd(theta)/sind(theta)));
    end

    x1 = halfsize;
    y1 = round(x1*(sind(theta)/cosd(theta)));

    if y1 > halfsize
        y1 = halfsize;
        x1 = round(y1*(cosd(theta)/sind(theta)));
    end

    pt0 = [halfsize-y0+1 halfsize+x0+1];
    pt1 = [halfsize-y1+1 halfsize+x1+1];

    mask = drawline(pt0,pt1,mask);
end

end

function [img] = drawline(pt0,pt1,orgimg)
img = orgimg;
linepts = getlinepts(pt0,pt1);
for i = 1:size(linepts,1)
   img(linepts(i,1),linepts(i,2)) = 1; 
end

end

function [linepts] = getlinepts(pt0,pt1)
% Return the points along the straight line connecting pt1 and pt2
if pt0(2) < pt1(2)
    x0 = pt0(2);    y0 = pt0(1);
    x1 = pt1(2);    y1 = pt1(1);
else
    x0 = pt1(2);    y0 = pt1(1);
    x1 = pt0(2);    y1 = pt0(1);
end

dx = x1 - x0;   dy = y1 - y0;
ind = 1;
linepts = zeros(numel(x0:x1),2);
step = 1;
if dx == 0 
   x = x0;
   if dy < 0,   step = -1;  end
   for y = y0:step:y1
        linepts(ind,:) = [y,x];
        ind = ind + 1;
   end
elseif abs(dy) > abs(dx)
    if dy < 0,  step = -1;  end
    for y = y0:step:y1
       x = round((dx/dy)*(y - y0) + x0);
       linepts(ind,:) = [y,x];
       ind = ind + 1;
    end
else
    for x = x0:x1
        y = round((dy/dx)*(x - x0) + y0);
        linepts(ind,:) = [y, x]; 
        ind = ind + 1;
    end
end

end

function [simg] = standardize(img,mask,wsize)

if (nargin == 2 || wsize == 0)
    simg = globalstandardize(img,mask);  
else
    img(mask == 0) = 0;
    img_mean = nlfilter(img,[wsize, wsize],@getmean);
    img_std = nlfilter(img,[wsize, wsize],@getstd);

    simg = (img - img_mean)./img_std;
    simg(img_std == 0) = 0;
    simg(mask == 0) = 0;
end

end

function [simg] = globalstandardize(img,mask)
usedpixels = double(img(mask==1));
m = mean(usedpixels);
s = std(usedpixels);

simg = zeros(size(img));
simg(mask == 1) = (usedpixels - m)/s;
end

function [m] = getmean(x)
usedx = x(x ~= 0);
m = mean(usedx);
end

function [s] = getstd(x)
usedx = x(x ~= 0);
s = std(usedx);
end

%% 9. thredshold (雙門檻)
% 使用Otsu演算法去除陰影及切割出鱗片輪廓
function [thredshold_img]=Thredshold (img)
    
    [row,lengths] = size(img);
    result_img = [row,lengths];

    % 取雙門檻(Otsu && 整張平均+整張最小值除以2)
    Ostu = (graythresh(img))*255;

    tmp_min = min(img(:));
    Mean = (sum(img(:))/(row*lengths));
    new_Mean = (Mean+tmp_min)/2;

    M = max(Ostu,new_Mean);
    m = min(Ostu,new_Mean);

    for i = 1:row
        for j = 1:lengths

            if img(i,j) >= M  %% 大於Max則變白
                result_img(i,j) = 255;
            elseif img(i,j) <= m %% 小於mix則變黑
                result_img(i,j) = 0;
            else
                result_img(i,j) = 128; %% 其餘則變灰
            end

        end
    end

    thredshold_img = uint8(result_img);
        
end

%% 10. fillhole
% 填補小洞
function [fill_img]=Fillhole (img)
        
    result_fill = bwareaopen(~img,150,4); %小於150的區塊去除 使用4鄰居算法
    fill_img=(~result_fill)*255;
        
end

%% 11. closing
% 將斷開或漏洞的鱗狀區塊填補
function [closing_img]=Closing(img)

    [row,lengths] = size(img);
    con = 1;
    times = 0;

    while con > 0
        count = 0;

        for i = 2:row-1
            for j = 2:lengths-1

                if img(i,j)==128
                    if img(i-1,j-1)==255 || img(i-1,j)==255 || img(i-1,j+1)==255 || img(i,j-1)==255 || img(i,j+1)==255 || img(i+1,j-1)==255 || img(i+1,j)==255 || img(i+1,j+1)==255
                        img(i,j)=255;
                        count = count+1;
                    end
                end

            end
        end
        %   count
        times = times+1;
        y(1,times) = count;
        if count == 0 % count為0時,就停止
            con = 0;
        end
    end

    for i = 1:row
        for j = 1:lengths

            if img(i,j)==128
                img(i,j)=0;
            end

        end
    end

    %% closing

    k = 1; %% 做幾次
    se = strel('square',3); %% 方形(square) / 圓形(disk) / 3x3 / 5x5
    %%創建指定形狀的對應結構元素 SE = strel(shape,parameters)%%
    %%parameters控制大小%%

    for count = 1:1:k
        img = imdilate(img,se);
    end

    for count = 1:1:k
        result_img = imerode(img,se);
    end

    closing_img = uint8(result_img);

end

%% 12. thining
% 細線化
function [thin_img]=Thining (img)

    img = img / 255; % 轉為0~1之間
    [row,lengths] = size(img);
    tmp = zeros(row,lengths);
    con = 1;
    times = 0;

    %二值化
    for i = 1:row
        for j = 1:lengths

            if ( img(i,j) > 0 )
                tmp(i,j) = 1;

            end
        end
    end

    while con > 0

        count = 0;

        for i = 2:row-1
            for j = 2:lengths-1

                if(tmp(i,j)==1)

                    if ( tmp(i-1,j-1)==0&&tmp(i,j-1)==0&&tmp(i+1,j-1)==0&&tmp(i,j)==1&&tmp(i-1,j+1)==1&&tmp(i,j+1)==1&&tmp(i+1,j+1)==1 )
                        tmp(i,j)=0;
                        count = count+1;

                    elseif ( tmp(i-1,j-1)==1&&tmp(i+1,j-1)==0&&tmp(i-1,j)==1&&tmp(i,j)==1&&tmp(i+1,j)==0&&tmp(i-1,j+1)==1&&tmp(i+1,j+1)==0 )
                        tmp(i,j)=0;
                        count = count+1;

                    elseif ( tmp(i,j-1)==1&&tmp(i-1,j)==1&&tmp(i,j)==1&&tmp(i+1,j)==0&&tmp(i,j+1)==0&&tmp(i+1,j+1)==0 )
                        tmp(i,j)=0;
                        count = count+1;

                    elseif ( tmp(i,j-1)==0&&tmp(i+1,j-1)==0&&tmp(i-1,j)==1&&tmp(i,j)==1&&tmp(i+1,j)==0&&tmp(i,j+1)==1 )
                        tmp(i,j)=0;
                        count = count+1;

                    elseif ( tmp(i-1,j-1)==1&&tmp(i,j-1)==1&&tmp(i+1,j-1)==1&&tmp(i,j)==1&&tmp(i-1,j+1)==0&&tmp(i,j+1)==0&&tmp(i+1,j+1)==0 )
                        tmp(i,j)=0;
                        count = count+1;

                    elseif ( tmp(i-1,j-1)==0&&tmp(i+1,j-1)==1&&tmp(i-1,j)==0&&tmp(i,j)==1&&tmp(i+1,j)==1&&tmp(i-1,j+1)==0&&tmp(i+1,j+1)==1 )
                        tmp(i,j)=0;
                        count = count+1;

                    elseif ( tmp(i,j-1)==1&&tmp(i-1,j)==0&&tmp(i,j)==1&&tmp(i+1,j)==1&&tmp(i-1,j+1)==0&&tmp(i,j+1)==0 )
                        tmp(i,j)=0;
                        count = count+1;

                    elseif ( tmp(i-1,j-1)==0&&tmp(i,j-1)==0&&tmp(i-1,j)==0&&tmp(i,j)==1&&tmp(i+1,j)==1&&tmp(i,j+1)==1 )
                        tmp(i,j)=0;
                        count = count+1;
                    end
                end
            end

        end
        %   count
        times = times+1;
        y(1,times) = count;
        if count == 0 % count為0時,就停止
            con = 0;
        end
    end

    thin_img = tmp*255;

end

%% 13. connect
% 把所有破碎的鱗狀輪廓連線
function [connect]=Connect(img, connect_TH, connect_Maxlen)

%     figure(1)
%     imshow(img);
    img_tmp=img; 
    BW=img_tmp/255;
    [L, ~] = bwlabeln(BW, 4); % 線段編號
    BW2 = bwmorph(BW,'endpoints'); % 找端點
    [end_y, end_x]=find(BW2==1);
    
    win=10;
    for i=1:length(end_x)
        endpoints=[end_x(i), end_y(i)];
        [branchpoints, line_len]=find_branchpoints(BW, endpoints, connect_TH);
        
        p=[end_x(i)-win, end_x(i)+win; end_y(i)-win, end_y(i)+win];
        p(p<1)=1;
        if p(1, 1)>size(BW, 2)
            p(1, 1)=size(BW, 2);
        end
        if p(1, 2)>size(BW, 2)
            p(1, 2)=size(BW, 2);
        end
        if p(2, 1)>size(BW, 1)
            p(2, 1)=size(BW, 1);
        end
        if p(2, 2)>size(BW, 1)
            p(2, 2)=size(BW, 1);
        end
        range_x=p(1, 1):p(1, 2);
        range_y=p(2, 1):p(2, 2);
        BW2(end_y(i), end_x(i))=0;
        block=BW2(range_y, range_x);
        [tmp_y, tmp_x]=find(block==1);
        BW2(end_y(i), end_x(i))=1;
        BW(end_y(i), end_x(i))=2;
        block_tmp=BW(range_y, range_x);
        [c_y, c_x]=find(block_tmp==2);
        BW(end_y(i), end_x(i))=1;
        label=L(range_y, range_x);
        
        D=zeros(length(tmp_x), 1);
        for j=1:length(tmp_x)
            if label(tmp_y(j), tmp_x(j))~=L(end_y(i), end_x(i))
                D(j)=points_distance([c_x, c_y], [tmp_x(j), tmp_y(j)]);
            else
                D(j)=1000;
            end           
        end
        [Min, I]=min(D);
        
        theta=atan2d((end_y(i)-branchpoints(2)), (end_x(i)-branchpoints(1)));
        if Min~=1000 & Min~=0            
            connect_x=range_x(tmp_x(I));
            connect_y=range_y(tmp_y(I));
            
            connect_theta=atan2d((connect_y-end_y(i)), (connect_x-end_x(i)));
            
            if abs(connect_theta-theta)<30
                img_tmp=connect_points(img_tmp, [end_x(i), end_y(i)], [connect_x, connect_y]);
            end
        end
        if line_len>connect_TH
                ending=[end_x(i)+connect_Maxlen*cos(theta), end_y(i)+connect_Maxlen*sin(theta)];
%                 a=(end_y(i)-branchpoints(2))/(end_x(i)-branchpoints(1));
%                 b=end_y(i)-a*end_x(i);
%                 ending=[(end_x(i)+connect_TH), (a*(end_x(i)+connect_TH)+b)];
                img_tmp=connect_points(img_tmp, [end_x(i), end_y(i)], ending);
        end
    end
%     figure(2)
%     imshow(img_tmp);
    
    connect=img_tmp;

    
end

function [branchpoints, line_len]=find_branchpoints(BW, endpoints, connect_TH)
    s=2;
    start=endpoints;
    B=BW;
    line_len=1;
    while s==2 
        p=[start(1)-1, start(1)+1; start(2)-1, start(2)+1];
        p(p<1)=1;
        if p(1, 1)>size(BW, 2)
            p(1, 1)=size(BW, 2);
        end
        if p(1, 2)>size(BW, 2)
            p(1, 2)=size(BW, 2);
        end
        if p(2, 1)>size(BW, 1)
            p(2, 1)=size(BW, 1);
        end
        if p(2, 2)>size(BW, 1)
            p(2, 2)=size(BW, 1);
        end
        x_p=p(1, 1):p(1, 2); y_p=p(2, 1):p(2, 2); 
        block=B(y_p, x_p);
        total=sum(sum(block));
        B(start(2), start(1))=0;
        if total==2
            block=B(y_p, x_p);
            [y_b, x_b]=find(block==1);
            start=[x_p(x_b), y_p(y_b)];
            line_len=line_len+1;
            if line_len==connect_TH
                branch=start;
            end
        else
            branchpoints=start;
        end
        s=total;
    end
    if line_len>connect_TH
        branchpoints=branch;
    end
end

function [result]=connect_points(img, starting, ending)
    result=img;
    if (starting(1)-ending(1))~=0
        a=(ending(2)-starting(2))/(ending(1)-starting(1));
        b=starting(2)-a*starting(1);
        for n=min(starting(1), ending(1)):0.05:max(starting(1), ending(1))
            m=n*a+b;
            if round(m)<1
                m=1;
            end
            if  round(m)>size(img, 1)
                m=size(img, 2);
            end
            if round(n)<1
                n=1;
            end
            if round(n)>size(img, 2)
                n=size(img, 2);
            end
            result(round(m),round(n))=255;
        end
    else
        for m=min(starting(2), ending(2)):0.05:max(starting(2), ending(2))
            n=round(starting(1));
            if round(m)<1
                m=1;
            end
            if  round(m)>size(img, 1)
                m=size(img, 2);
            end
            if round(n)<1
                n=1;
            end
            if round(n)>size(img, 2)
                n=size(img, 2);
            end
            result(round(m),round(n))=255;
        end
    end   
end

%% 14. spur
% 去分岔
function [spur_img]=Spur(img)
    
%     name = [path '13.connection\final\connection_' name_part_num  '_EHr=' num2str(EH_r)  '_LG=' num2str(LG_Gamma) '_LW=' num2str(LG_Wsize) '_MW=' num2str(Min_Wsize) '_RW=' num2str(RL_Wsize)  '_MsW=' num2str(MS_Wsize) '_Mstep=' num2str(MS_step)  '_TH=' num2str(TH) '_Maxlen=' num2str(connect_Maxlen) '.bmp'];
% 
%     if (exist(name,'file'))
%         img = imread(name);
%     else
%         img=-1;
%     end
%     if (img~=-1)
        connect_img=img;
        [row,lengths] = size(img);

        img_fill = bwareaopen(~img,10,4); %填埔小洞，降低小區塊的數量

        img=~img_fill;

        result_area = bwmorph(img,'thin',Inf); %細線化
        result_area = bwmorph(result_area,'spur',Inf); %去分岔

        %                      (1)Determine the connected components.
        label_fill = bwlabeln(result_area,8);
        %                      (2)Compute the area of each component.
        S = regionprops(label_fill, 'Area');
        %                      (3)Remove small objects.
        result_area =ismember(label_fill, find([S.Area] >2));
        result_area=result_area*255;
        for x=2:1:row-1
            for y=2:1:lengths-1
                if(result_area(x,y)==255)

                    if connect_img(x,y)==120
                        result_area(x,y)=120;
                    end
                    if connect_img(x,y)==200
                        result_area(x,y)=100;
                    end

                end
            end
        end
        
        spur_img=result_area;
        
%         new=[path,'14.spur\'];
% 
%         if ~exist(new, 'dir')
%             % Folder does not exist so create it.
%             mkdir(new);
%         end
% 
%         result_name= [new 'newspur_' name_part_num  '_EHr=' num2str(EH_r)  '_LG=' num2str(LG_Gamma) '_LW=' num2str(LG_Wsize) '_MW=' num2str(Min_Wsize) '_RW=' num2str(RL_Wsize)  '_MsW=' num2str(MS_Wsize) '_Mstep=' num2str(MS_step)  '_TH=' num2str(TH) '_Maxlen=' num2str(connect_Maxlen) '.bmp'];
%         imwrite(spur_img, result_name);

end

%% 15. adjustment after segmentation
% 切割後的調整階段
function [chrom_roiarea, chrom_whitearea, result_img]=Adjust_seg(img, L_n, lower, upper)
%     name =[path '14.spur\newspur_' name_part_num  '_EHr=' num2str(EH_r)  '_LG=' num2str(LG_Gamma) '_LW=' num2str(LG_Wsize) '_MW=' num2str(Min_Wsize) '_RW=' num2str(RL_Wsize)  '_MsW=' num2str(MS_Wsize) '_Mstep=' num2str(MS_step)  '_TH=' num2str(TH) '_Maxlen=' num2str(connect_Maxlen) '.bmp']; %new data
%     mask_name = [mask_path name_part_num '.bmp'];

%     if (exist(name,'file'))
%         img = imread(name);
%     else
%         img=-1;
%     end

%     if (exist(mask_name,'file'))
%         mask = imread(mask_name);  %區分前景背景
%     else
%         mask=-1;
%     end
%     %             [row,lengths] = size(img);
%     if (mask~=-1)
        
%         figure(1)
%         imshow(img);
        
        [row, column]=size(img);
        
        img_fill = bwareaopen(~img,8,4); % 填補極小洞 % %作用：消除~img中小於10的面積區塊，使用4鄰居算法。目的：減少小區塊面積
        img = ~img_fill; %img為 黑底
        
%         figure(2)
%         imshow(~img);
        
        img_tmp=~img;
        
        % 把左右鼻孔中心點區塊與背景合併(消除區塊周圍線段)
        label_img = bwlabel(img_tmp, 4);
        label_leftnose=label_img(round(L_n(2, 1)), round(L_n(1, 1)));
        label_rightnose=label_img(round(L_n(2, 3)), round(L_n(1, 3)));
        if label_leftnose~=1
            [merge_idx]=region_merge(label_img, label_leftnose);
            if ~isempty(merge_idx)
            for r=1:length(merge_idx(:, 2))
                    img_tmp(merge_idx(r, 2), merge_idx(r, 1))=1;
            end
            end
        end
        if label_rightnose~=1
            [merge_idx]=region_merge(label_img, label_rightnose);
            if ~isempty(merge_idx)
            for r=1:length(merge_idx(:, 2))
                    img_tmp(merge_idx(r, 2), merge_idx(r, 1))=1;
            end
            end
        end
        
%         figure(3)
%         imshow(img_tmp);

        [L, num] = bwlabel(img_tmp, 4);        
        Area_label=regionprops(L,'Area');
        new_Area_label= Area_label(2:end, :);
        label_std=std([new_Area_label.Area]);%計算面積標準差
        label_mean=mean([new_Area_label.Area]); %計算面積總平均
        
        low_area=lower*label_std+label_mean;
                
        for i=2:num
            if Area_label(i).Area<low_area
                [merge_idx]=region_merge(L, i);
                if ~isempty(merge_idx)
                for r=1:length(merge_idx(:, 2))
                    img_tmp(merge_idx(r, 2), merge_idx(r, 1))=1;
                end 
                end
            end
        end
        
%         figure(9)
%         imshow(img_tmp);
        
        result_area = bwmorph(~img_tmp,'thin',Inf);%細線化
        result_area = bwmorph(result_area,'spur',Inf); %去分岔 %result_area為黑底
        result_area=~result_area;
        img_tmp= bwareaopen(result_area,2,4);

%         figure(10)
%         imshow(img_tmp); 
        chrom_roiarea=0;
        chrom_blackarea=0;
        chrom_whitearea=0;
        
        [L, num] = bwlabel(img_tmp, 4);
        Area_label=regionprops(L,'Area');
        new_Area_label= Area_label(2:end, :);
        label_std=std([new_Area_label.Area]);%計算面積標準差
        label_mean=mean([new_Area_label.Area]); %計算面積總平均

        up_area=upper*label_std+label_mean;

        D_n=points_distance(L_n(:, 1), L_n(:, 3));
        R_base=[L_n(1, 1), L_n(1, 3); (min(L_n(2, 1), L_n(2, 3))-(D_n*0.5)), max(L_n(2, 1), L_n(2, 3))]; % 基準範圍(R_base)
        
%         figure(10)
%         imshow(img_tmp);
%         hold on
%         plot([R_base(1,1), R_base(1,2), R_base(1,2), R_base(1,1), R_base(1,1)], [R_base(2, 2), R_base(2, 2), R_base(2, 1), R_base(2, 1), R_base(2, 2)], '-o','color','r','LineWidth',2);
%         plot(L_n(1,:), L_n(2,:), '-o','color','b','LineWidth',2);
%         hold off
        
%         roi=[];
%         black=[];
        
        for j=2:num
            
            [p, q]=find(L==j);
            n=numel(p);
            x=sum(q)/n;	%編號i 中心點X軸座標 (編號區域所有點X值加總/區域總個數)
            y=sum(p)/n;	%編號i 中心點Y軸座標 (編號區域所有點Y值加總/區域總個數)
            
            if (R_base(1, 1)<=x && R_base(1, 2)>=x && R_base(2, 1)<=y && R_base(2, 2)>=y)
                chrom_roiarea=chrom_roiarea+Area_label(j).Area;
                if Area_label(j).Area>up_area
                    img_tmp(find(L==j))=0;
                    chrom_blackarea=chrom_blackarea+Area_label(j).Area;
%                     black=[black, j];
                else
                    chrom_whitearea=chrom_whitearea+Area_label(j).Area;
                end
            else
                if Area_label(j).Area>up_area
                    img_tmp(find(L==j))=0;
                end
            end
                
%             figure(12)
%             imshow(img_tmp);
%             hold on
%             plot([R_base(1,1), R_base(1,2), R_base(1,2), R_base(1,1), R_base(1,1)], [R_base(2, 2), R_base(2, 2), R_base(2, 1), R_base(2, 1), R_base(2, 2)], '-o','color','r','LineWidth',2);
%             plot(L_n(1,:), L_n(2,:), '-o','color','b','LineWidth',2);
%             hold off
%             text(x, y, int2str(j), 'color', 'g');

        end
        
        img_tmp = ~bwmorph(~img_tmp,'clean',Inf); % 去除影像中孤立的亮點(pixel value=1 ==> 0, 故需反白)
        
        result_img=img_tmp;
        
%         new=[path,'15.adjust_seg\'];
%         
%         if ~exist(new, 'dir')
%             % Folder does not exist so create it.
%             mkdir(new);
%         end
%         
%         result_name= [new 'adjust_' name_part_num  '_EHr=' num2str(EH_r)  '_LG=' num2str(LG_Gamma) '_LW=' num2str(LG_Wsize) '_MW=' num2str(Min_Wsize) '_RW=' num2str(RL_Wsize)  '_MsW=' num2str(MS_Wsize) '_Mstep=' num2str(MS_step)  '_TH=' num2str(TH) '_Maxlen=' num2str(connect_Maxlen)  '_lower=' num2str(lower) '_upper=' num2str(upper) '.bmp'];
%         imwrite(result_img,result_name);
end

% 區塊合併 (消除鄰近區塊周圍線段，合併基準為相鄰區塊線段最長者)
function [merge_idx]=region_merge(L, i)
    merge_idx=[];
    [row, column]=size(L);
    BW=L;
    BW(BW~=i)=0;
    B = bwboundaries(BW);
    if ~isempty(B)
        bound=cell2mat(B(1,1)); % 邊緣座標(y, x)，最靠近黑線的部分


%     figure(5)
%     imshow(BW);
%     hold on;
%     plot(bound(:, 2), bound(:, 1), 'o', 'color', 'r', 'LineStyle', 'none');
%     hold off;
    
    bound_around=[]; % 欲合併區塊之周圍線段(黑線)座標

    for j=1:length(bound(:, 1))
        
        [x_idx, y_idx]=window_idx(bound(j, 2), bound(j, 1), column, row, 1);

        a=round(x_idx(1)):round(x_idx(2));
        b=round(y_idx(1)):round(y_idx(2));
        [q, p]=find(L(b, a)==0);

        tmp=[];
        for c=1:length(q)
            t=[a(p(c)), b(q(c))];
            tmp=[tmp; t];
        end
        bound_around=[bound_around; tmp];

    end

    bound_around = unique(bound_around,'rows');

    bound_labels=[];

    for k=1:length(bound_around(:, 1))
        
        [x_idx, y_idx]=window_idx(bound_around(k, 1), bound_around(k, 2), column, row, 1);
        tempBlock=L(round(y_idx(1)):round(y_idx(2)), round(x_idx(1)):round(x_idx(2)));

        u=unique(tempBlock);
        u=u(u~=0);  % label=0 為區塊周圍線段

        if length(u)==2 && ~isempty(find(u==i))
            index=u(u~=i);
            tmp=[bound_labels; bound_around(k, 1), bound_around(k, 2), index];
            bound_labels=tmp;
        end

    end

    l=unique(bound_labels(:, 3));
    count=[];
    for s=1:length(l)
        n=find(bound_labels(:, 3)==l(s));
        count=[count; l(s), length(n)];
    end
    [M]=sortrows(count, 2, 'descend');

    [idx]=find(bound_labels(:, 3)==M(1, 1));
    
    for r=1:length(idx)
        merge_idx(r, :)=[bound_labels(idx(r), 1), bound_labels(idx(r), 2)];
    end
     end
end

