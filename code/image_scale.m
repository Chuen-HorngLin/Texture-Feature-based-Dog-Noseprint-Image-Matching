%% �վ�v���ؤo
%% owner: yang shu chun
clear all;

path='E:\�M�D��s���\���βE\darknet_�ҪO1_200409\darknet\darknet-master\build\darknet\x64\dognose_all\';

files = dir([path 'cut\*.bmp']);
for i=1:length(files)
    img_name=files(i).name;%��o�ɦW�r��
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

%% �v���Y��(��z�G���u�ʴ��Ȫk)
% input
% (1) img     : ��J�v�� (�i��RGB��gray)
%      type   : 'scale', �Y���� 
%      value : �Y���
% (2)  img    : ��J�v�� (�i��RGB��gray)
%       type   : 'size', �Y�����
%       value : [��XY�b����, ��XX�b����] (��@�i���ŭȡA��J-1)
function [output]=imgResize(img, type, value)
    [~, ~, n]=size(img); % �P�_���m��άO�Ƕ��v��
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
    if strcmp(type, 'scale') % ����Y��v��
        scale(1)=value(1);
        scale(2)=value(1); 
        result=zeros(floor(size(img, 1)*scale(1)), floor(size(img, 2)*scale(2)));
    end
    if strcmp(type, 'size') % �T�w���e�v��
        if (value(1)==-1) % X�b���שT�w�AY�b����Ҫ���
            scale(2)=value(2)/size(img, 2); % �Y����
            scale(1)=scale(2); 
            value(1)=floor(size(img, 1)*scale(1)); % ��X�v��Y�b����
        elseif (value(2)==-1) % Y�b���שT�w�AX�b����Ҫ���
            scale(1)=value(1)/size(img, 1);
            scale(2)=scale(1);
            value(2)=floor(size(img, 2)*scale(2)); % ��X�v��X�b����
        else
            scale(1)=value(1)/size(img, 1); % Y�b���שT�w�A�p���Y����
            scale(2)=value(2)/size(img, 2); % X�b���שT�w�A�p���Y����
        end
        result=zeros(value(1), value(2)); % ��X�v�����e�j�p
    end
    
    % ���F�קK�|�P�x�}�b�M�g�y�ЮɵL�k��������m�A�U�|�P�U�X�i1pixel�j�p�A�Ȩ��P��F��a�J
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
            %  ����k�O��X�Ϲ��ϬM�g�^��l�Ϲ��y��
            x_p=x/scale(2)+1;   % �]�������Ϲ�(temp)���|�P���X�i1pixel�A�GX�b�_�l�I+1
            y_p=y/scale(1)+1;   % �]�������Ϲ�(temp)���|�P���X�i1pixel�A�GY�b�_�l�I+1
            x_0=floor(x_p); x_1=min(x_0+1, (size(temp, 2))); % �קK�M�g�W�X�d��A�W���̰���temp��X�b����
            y_0=floor(y_p); y_1=min(y_0+1, (size(temp, 1))); % �קK�M�g�W�X�d��A�W���̰���temp��Y�b����
            
            f_x0y0=temp(y_0, x_0); f_x0y1=temp(y_1, x_0); % �W�U���y�Щ��Ϫ�������
            f_x1y0=temp(y_0, x_1); f_x1y1=temp(y_1, x_1);
            
            t1=((x_1-x_p)/(x_1-x_0))*(f_x0y0)+((x_p-x_0)/(x_1-x_0))*(f_x1y0); % X�b��V�u�ʴ���
            t2=((x_1-x_p)/(x_1-x_0))*(f_x0y1)+((x_p-x_0)/(x_1-x_0))*(f_x1y1);
            t=((y_1-y_p)/(y_1-y_0))*(t1)+((y_p-y_0)/(y_1-y_0))*(t2); % Y�b��V�u�ʴ���

            result(y, x)=t;
        end
    end
    result=uint8(result);
end
