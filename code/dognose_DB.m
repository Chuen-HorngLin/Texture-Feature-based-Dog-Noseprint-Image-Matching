%% 建立狗鼻紋DB (提取影像資料、基準點、特徵值)
%% owner: yang shu chun
clear all;
warning('off','all')
%%
dir='base';

path=['..\data\' dir '\'];
result='output\15. result\15. result_';
txt='txt\rotated\rotated_position_';
img_code='_EHr=0.9_LG=4_LW=60_MW=1_RW=10_MsW=40_Mstep=2_TH=11_Maxlen=11_lower=-0.6_upper=1.45.bmp';

db_start = datestr(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
disp(['Start set DB: ', db_start])
get_dognose_DB(path, dir, result, txt, img_code);
db_end = datestr(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
disp(['End set DB: ', db_end])


function get_dognose_DB(path, dir_p, img_p, txt_p, img_code)

    files=dir([path img_p '*.bmp']);
    for k=1:length(files)

        temp=strrep(files(k).name, img_code, '');
        name_part_num=strrep(temp, '15. result_', '');

        img_name = [path img_p name_part_num img_code];
        if (exist(img_name,'file'))
            img = imread(img_name);
        else
            continue;
        end

        txt_name = [path txt_p name_part_num '.txt'];
        if (exist(txt_name,'file'))
            txt = fopen(txt_name, 'r');
        else
            continue;
        end

        L_n = cell2mat( textscan(txt, '%f%f%f', 'Delimiter', ',', 'collectoutput', 1) ); % 基準線座標

        D_n=points_distance(L_n(:, 1), L_n(:, 3)); % 基準線距離

        [~, column]=size(img);
        R_n=[1, column; 1, max(L_n(2, :))]; % 比較範圍(R_n)
        R_base=[L_n(1, 1), L_n(1, 3); (min(L_n(2, 1), L_n(2, 3))-(D_n*0.5)), max(L_n(2, 1), L_n(2, 3))]; % 基準範圍(R_base)

        [L, num] = bwlabel(img, 4); % 將每個區塊編號
        [base, region]=region_extract(L, num, L_n, D_n, R_base, R_n); % 基準範圍與比較範圍內鱗狀區塊之資料


        id=strsplit(name_part_num, '_');
        i=str2double(id(1)); j=str2double(id(2));
        name=[i, j];
        
        info.name=name;         % 名稱
        info.L=L;               % 原始影像
        info.L_n=L_n;           % 基準線(左邊鼻孔中心點，基準點，右邊鼻孔中心點)
        info.D_n=D_n;           % 基準線長度
        info.base=base;         % 基準範圍內鱗狀區塊(編號，與基準點距離，與基準線夾角，位於基準點之左[-1]右[1]方)
        info.region=region;     % 比較範圍內鱗狀區塊(編號，特徵值[長條程度R_S、不規則程度R_D、Moment特徵、角度方向比例R_0, R_45、R_90、R_135]，中心點座標)
        
        new=['database_original\' dir_p '\'];
        if ~exist(new, 'dir')
            mkdir(new);
        end

        save([new 'dognose_' dir_p '_db-' name_part_num '.mat'], '-struct', 'info');

        fclose('all');
    end
end

%%  求目標範圍(R)內各個編號區塊中心點(C)X, Y座標、與基準點距離、L_n與L_dist夾角角度
function [base, region]=region_extract(L, num, L_n, D_n, R_base, R_n)
    
    base=cell2table(cell(0, 4));
    base.Properties.VariableNames={'Label', 'Distance', 'Theta', 'Direct'};
    region=cell2table(cell(0, 3));
    region.Properties.VariableNames={'Label', 'Center', 'Feature'};
    
    num_1=1; num_2=1;
    for i=2:num
        [p,q]=find(L==i);	% 找到編號i 內所有座標[y, x]
        n=numel(p);
        x=sum(q)/n;	% 編號i 中心點X軸座標 (編號區域所有點X值加總/區域總個數)
        y=sum(p)/n;	% 編號i 中心點Y軸座標 (編號區域所有點Y值加總/區域總個數)

        if (R_base(1, 1)<=x && R_base(1, 2)>=x && R_base(2, 1)<=y && R_base(2, 2)>=y)
            dist=points_distance(L_n(:, 2), [x, y]);       % L_dist區域中心點與基準線中心點連線距離
            if (x<L_n(1, 2))
                v1=[x, y]-[L_n(1, 2), L_n(2, 2)];
                v2=[L_n(1, 1), L_n(2, 1)]-[L_n(1, 2), L_n(2, 2)];
                theta=acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
                direct=-1;
            else
                v1=[x, y]-[L_n(1, 2), L_n(2, 2)];
                v2=[L_n(1, 3), L_n(2, 3)]-[L_n(1, 2), L_n(2, 2)];
                theta=acosd(sum(v1.*v2)/(norm(v1)*norm(v2)));
                direct=1;
            end

            base(num_1, :)=table(i, dist, theta, direct);
            
            num_1=num_1+1;
        end

        if (R_n(1, 1)<=x && R_n(1, 2)>=x && R_n(2, 1)<=y && R_n(2, 2)>=y) % 判別區塊是否在比較範圍內
            dist=points_distance(L_n(:, 2), [x, y]);       % L_dist區域中心點與基準線中心點連線距離
            [R_S, R_D, Moment, R_0, R_45, R_90, R_135]=region_feature(L, i, D_n, dist);            
            
            region(num_2, :)=table(i, [x, y], [R_S, R_D, Moment, R_0, R_45, R_90, R_135]);
            
            num_2=num_2+1;
        end
    end
end

%% 求區塊特徵值: 長條程度R_S、不規則程度R_D、Moment特徵、角度方向比例R_0, R_45、R_90、R_135
function [R_S, R_D, Moment, R_0, R_45, R_90, R_135]=region_feature(L, label, D_n, L_dist)

    L(L~=label)=0;
    BW = imbinarize(L); % 目標區塊範圍

    item=regionprops(BW, 'Area', 'Perimeter');                  % 該區塊之面積、周長
    [D_c, D, slope]=greatest_diameter(BW);                      % 區塊最長距離之中心點座標、長度、斜率

    R_S=(4*item.Area)/(pi*(D^2));                               % 長條程度R_S
    R_D=(4*pi*item.Area)/(item.Perimeter^2);                    % 不規則程度R_D
    Moment=sqrt(item.Area*L_dist)/(D_n*2);                      % Moment特徵，用來正規化每一區塊

    [D_0, D_45, D_90, D_135]=degree_dist(slope, BW, D_c, D);    % 角度方向長度D_0、D_45、D_90、D_135

    if abs(D-D_0)<1
        D=D_0;
    end

    R_0=D_0/D;          % 角度方向比例R_0
    R_45=D_45/D;        % 角度方向比例R_45
    R_90=D_90/D;        % 角度方向比例R_90
    R_135=D_135/D;      % 角度方向比例R_135
    
end

%% 求區塊最遠距離
function [D_c, D, slope]=greatest_diameter(BW)

    B = bwboundaries(BW);
    bound=cell2mat(B(1,1));

    [dist]=pdist2(bound, bound);
    D=max(dist(:));
    [a, b] = find(dist==D);

    D_c=[bound(a(1), 2),((bound(a(1), 2)+bound(b(1), 2))*0.5), bound(b(1), 2);
        bound(a(1), 1),((bound(a(1), 1)+bound(b(1), 1))*0.5), bound(b(1), 1)];
    slope=(D_c(2, 3)-D_c(2, 1))/(D_c(1, 3)-D_c(1, 1));

end

%% 求夾角之線段距離(不包含黑色部分)
function [D_0, D_45, D_90, D_135]=degree_dist(slope, BW, D_c, D)
    
    angle = atan(slope);
    th = 0:pi/4:2*pi;   %以最長距離之中心點取夾角 0, 45, 90, 135 180 225 270 315 八個方向的斜率

    % 計算各方位的距離(不包含區塊外的部分，指的是如果有內凹，超出區塊外的不算)
    D_0=degree_dist_count((th(1)+angle), BW, D_c, D)+degree_dist_count((th(5)+angle), BW, D_c, D);      % 0, 180
    D_45=degree_dist_count((th(2)+angle), BW, D_c, D)+degree_dist_count((th(6)+angle), BW, D_c, D);    % 45, 225
    D_90=degree_dist_count((th(3)+angle), BW, D_c, D)+degree_dist_count((th(7)+angle), BW, D_c, D);    % 90, 270
    D_135=degree_dist_count((th(4)+angle), BW, D_c, D)+degree_dist_count((th(8)+angle), BW, D_c, D);   % 135, 315

end

%% 求中心點取夾角D之線段長度(不包含黑色部分)
function [dist]=degree_dist_count(degree, BW, D_c, D)

    % 以中心取長度為D/2，取角度degree的方向，涵蓋的所有座標(間距取 D*0.5)
    d=0:0.005:(D*0.5);
    x_p = d .* cos(degree) + D_c(1, 2);
    y_p = d .* sin(degree) + D_c(2, 2);

    [row, column]=size(BW);
    idx=1;
    for i=1:length(d)
        x=round(x_p(i)); y=round(y_p(i));
        if (x>=1 && x<=column) && (y>=1 && y<=row)
            line(idx)=BW(round(y), round(x)); %紀錄長度涵蓋的座標 值 (如果在範圍內為1 範圍外為0)
            idx=idx+1;
        end
    end
    line_temp=[line, 0];
    line_temp_1=[0, line];
    edge=line_temp-line_temp_1; % 找尋邊界點
    [start_p]=find(edge==1); % 找尋開始連線的地方 (如果有內凹，可能不只一個)
    [end_p]=find(edge==-1); % 找尋結束的地方
    dist=0;

    if (isempty(start_p)==0)
        for j=1:length(start_p)
            dist=dist+points_distance([x_p(start_p(j)), y_p(start_p(j))], [x_p(end_p(j)-1), y_p(end_p(j)-1)]);
        end
    end
    
end

%% 求座標p1 p2 之連線距離
function [dist]=points_distance(p1, p2)
    dist=sqrt(((p1(1)-p2(1))^2)+((p1(2)-p2(2))^2));
end
