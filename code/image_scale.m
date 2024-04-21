%% 調整影像尺寸
%% owner: yang shu chun
clear all;

path='E:\專題研究資料\楊舒淳\darknet_模板1_200409\darknet\darknet-master\build\darknet\x64\dognose_all\';

files = dir([path 'cut\*.bmp']);
for i=1:length(files)
    img_name=files(i).name;%獲得檔名字串
    img=imread([path 'cut\' img_name]);

    [row, column, ~]=size(img);
    
    if column>=2000
        c=column;
        img_temp=img;
        while c>2000
            img_temp = imgResize(img_temp, 'scale', 0.5);
            [r, c, ~]=size(img_temp);
        end
    elseif column<1000
        c=column;
        img_temp=img;
        while c<1000
            img_temp = imgResize(img_temp, 'scale', 2);
            [r, c, ~]=size(img_temp);
        end
    else
        img_temp=img;
    end
    
    new=[path 'resize\'];
    
    if ~exist(new, 'dir')
        % Folder does not exist so create it.
        mkdir(new);
    end
    
    result_name= [new img_name];
    imwrite(img_temp ,result_name);
    
end

%% 影像縮放(原理：雙線性插值法)
% input
% (1) img     : 輸入影像 (可為RGB或gray)
%      type   : 'scale', 縮放比例 
%      value : 縮放值
% (2)  img    : 輸入影像 (可為RGB或gray)
%       type   : 'size', 縮放長度
%       value : [輸出Y軸長度, 輸出X軸長度] (其一可為空值，輸入-1)
function [output]=imgResize(img, type, value)
    [~, ~, n]=size(img); % 判斷為彩色或是灰階影像
    if n==1
        output=imgResize_oneChannel(img, type, value);
    else
        red=img(:, :, 1);
        R=imgResize_oneChannel(red, type, value);
        green=img(:, :, 2);
        G=imgResize_oneChannel(green, type, value);
        blue=img(:, :, 3);
        B=imgResize_oneChannel(blue, type, value);
        output = cat(3, R, G, B);
    end
end

function [result]=imgResize_oneChannel(img, type, value)
    if strcmp(type, 'scale') % 比例縮放影像
        scale(1)=value(1);
        scale(2)=value(1); 
        result=zeros(floor(size(img, 1)*scale(1)), floor(size(img, 2)*scale(2)));
    end
    if strcmp(type, 'size') % 固定長寬影像
        if (value(1)==-1) % X軸長度固定，Y軸等比例長度
            scale(2)=value(2)/size(img, 2); % 縮放比例
            scale(1)=scale(2); 
            value(1)=floor(size(img, 1)*scale(1)); % 輸出影像Y軸長度
        elseif (value(2)==-1) % Y軸長度固定，X軸等比例長度
            scale(1)=value(1)/size(img, 1);
            scale(2)=scale(1);
            value(2)=floor(size(img, 2)*scale(2)); % 輸出影像X軸長度
        else
            scale(1)=value(1)/size(img, 1); % Y軸長度固定，計算縮放比例
            scale(2)=value(2)/size(img, 2); % X軸長度固定，計算縮放比例
        end
        result=zeros(value(1), value(2)); % 輸出影像長寬大小
    end
    
    % 為了避免四周矩陣在映射座標時無法找到對應位置，顧四周各擴張1pixel大小，值取周圍鄰邊帶入
    temp=zeros(size(img, 1)+2, size(img, 2)+2);
    temp(2:end-1, 2:end-1)=img;
    temp(1, 2:end-1)=img(1,1:end);
    temp(end, 2:end-1)=img(end,1:end);
    temp(2:end-1, 1)=img(1:end,1);
    temp(2:end-1, end)=img(1:end, end);
    temp(1, 1)=img(1, 1); temp(1, end)=img(1, end); 
    temp(end, 1)=img(end, 1); temp(end, end)=img(end, end);
    
    for y=1:size(result, 1)
        for x=1:size(result, 2)
            %  本方法是輸出圖像反映射回原始圖像座標
            x_p=x/scale(2)+1;   % 因為對應圖像(temp)為四周圍擴張1pixel，故X軸起始點+1
            y_p=y/scale(1)+1;   % 因為對應圖像(temp)為四周圍擴張1pixel，故Y軸起始點+1
            x_0=floor(x_p); x_1=min(x_0+1, (size(temp, 2))); % 避免映射超出範圍，上限最高為temp之X軸長度
            y_0=floor(y_p); y_1=min(y_0+1, (size(temp, 1))); % 避免映射超出範圍，上限最高為temp之Y軸長度
            
            f_x0y0=temp(y_0, x_0); f_x0y1=temp(y_1, x_0); % 上下限座標於原圖的像素值
            f_x1y0=temp(y_0, x_1); f_x1y1=temp(y_1, x_1);
            
            t1=((x_1-x_p)/(x_1-x_0))*(f_x0y0)+((x_p-x_0)/(x_1-x_0))*(f_x1y0); % X軸方向線性插值
            t2=((x_1-x_p)/(x_1-x_0))*(f_x0y1)+((x_p-x_0)/(x_1-x_0))*(f_x1y1);
            t=((y_1-y_p)/(y_1-y_0))*(t1)+((y_p-y_0)/(y_1-y_0))*(t2); % Y軸方向線性插值

            result(y, x)=t;
        end
    end
    result=uint8(result);
end
