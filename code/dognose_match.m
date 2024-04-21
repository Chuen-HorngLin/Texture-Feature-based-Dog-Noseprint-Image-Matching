%% 計算影像相似度 (取第一張、取前兩張、取前三張之相似度)
clear all;
warning('off','all')

%% dognose_match 參數
% group1
% w=[0.047297297	0.209459459	0.202702703	0.182432432	0.121621622	0.175675676	0.060810811];
% r=[0.950704225	0.792253521	0.411971831	0.728873239	0.792253521	0.443661972	0.38028169];
% R=9;
% image_threshold=0.35; 
% region_threshold=0.2;
% win=0.025;

% group2
% w=[0.123076923	0.076923077	0.261538462	0.123076923	0.153846154	0.092307692	0.169230769];
% r=[0.669291339	1.141732283	0.275590551	0.787401575	1.181102362	0.354330709	0.590551181];
% R=10;
% image_threshold=0.3; 
% region_threshold=0.225;
% win=0.025;

% group3
w=[0.047169811	0.188679245	0.264150943	0.122641509	0.20754717	0.075471698	0.094339623];
r=[0.373333333	0.266666667	0.64	0.746666667	0.826666667	0.64	0.506666667];
R=8;
image_threshold=0.3; 
region_threshold=0.2;
win=0.025;
%% folder path
group='group3';
DB_path=['C:\Users\User\Desktop\畢業論文\final\random_group\database_original\' group '\'];
base_files=dir([DB_path 'base\*.mat']);
train_files=dir([DB_path 'train\*.mat']);
test_files=dir([DB_path 'test\*.mat']);

tic
[train_fitness, train_fitness_result]=match_fitness(DB_path, base_files, 'train', train_files, w, r, region_threshold, image_threshold, win);
save(['match_result_' group '_train.mat'], 'train_fitness', 'train_fitness_result', '-v7.3');
toc
tic
[test_fitness, test_fitness_result]=match_fitness(DB_path, base_files, 'test', test_files, w, r, region_threshold, image_threshold, win);
save(['match_result_' group '_test.mat'], 'test_fitness', 'test_fitness_result', '-v7.3');
toc

%% 計算整體準確率
function [fitness, fitness_result]=match_fitness(DB_path, base_files, t, files, w, r, region_threshold, image_threshold, win)
    fitness_result=cell(length(files), 5);
    score=[0, 0, 0];

    for i=1:length(files)
        temp=strrep(files(i).name, 'dognose_test_db-', '');
        name_part_num=strrep(temp, '.mat', '');

        test_db=load([DB_path t '\' files(i).name]);
        
        acc=cell(length(base_files), 3);
        parfor j=1:length(base_files)
            base_db=load([DB_path 'base\' base_files(j).name]);
            [match_length_a, match_base]=image_match(test_db, base_db, win);    % 測試樣本之基準範圍內鱗狀區塊數量、輸出測試樣本對應到基準樣本的候補鱗狀區塊特徵 
            [match_length_b, match_test]=image_match(base_db, test_db, win);    % 基準樣本之基準範圍內鱗狀區塊數量、輸出基準樣本對應到測試樣本的候補鱗狀區塊特徵
            match_length=2*(1/((1/match_length_a)+(1/match_length_b)));         % 基準範圍內鱗狀區塊數量取調和參數
            
            match_b=get_match_region(match_base, w, r, region_threshold);       % 計算最相似的鱗狀區塊(測試樣本對應基準樣本)
            match_t=get_match_region(match_test, w, r, region_threshold);       % 計算最相似的鱗狀區塊(基準樣本對應測試樣本)
            s_1=0;
            for m=1:length(match_b(:, 1))
                index=find(match_t(:, 1)==match_b(m, 2));
                if (isempty(index)~=1 && length(index)==1)
                    if (match_t(index, 2)==match_b(m, 1)) % 兩張對應到相同編號
                        s_1=s_1+1;
                    else
                        match_b(m, 2)=-1;   % 編號不同則表示沒有對應到
                    end
                else
                    match_b(m, 2)=-1;
                end
            end
            match=s_1/match_length; % 兩張圖相似程度
            
            acc(j, :)={base_db.name(1), match, match_b};

            fclose('all');
        end
        [M]=sortrows(acc, 2, 'descend');
        
        % 鱗狀區塊相似度結果
        fitness_result(i, 1)={name_part_num}; % 測試樣本編號
        fitness_result(i, 2)={-1}; fitness_result(i, 3)={-1}; fitness_result(i, 4)={-1}; % 回應前三相似基準樣本
        fitness_result(i, 5)={M}; % 比對結果(與每一基準樣本鱗狀區塊相似度結果)
        
        m_id=[cell2mat(M(1, 1)), cell2mat(M(2, 1)), cell2mat(M(3, 1))];
        all_id=cell2mat(M(:, 1));
        match_id=find(all_id(:, 1)==test_db.name(1));
        
        % 相似度第一之基準樣本與測試樣本相同
        if match_id==1
            if (cell2mat(M(1, 2))>image_threshold)
                score(1)=score(1)+1;    score(2)=score(2)+1;    score(3)=score(3)+1;
                fitness_result(i, 2)={m_id(1)}; fitness_result(i, 3)={m_id(1)}; fitness_result(i, 4)={m_id(1)};
            end
        % 相似度第二之基準樣本與測試樣本相同
        elseif match_id==2
            if (cell2mat(M(2, 2))>image_threshold)
                score(2)=score(2)+1;    score(3)=score(3)+1;
                fitness_result(i, 2)={m_id(1)}; fitness_result(i, 3)={m_id(2)}; fitness_result(i, 4)={m_id(2)};
            end
        % 相似度第三之基準樣本與測試樣本相同
        elseif match_id==3
            if (cell2mat(M(3, 2))>image_threshold)
                score(3)=score(3)+1;
                fitness_result(i, 2)={m_id(1)}; fitness_result(i, 3)={m_id(1)}; fitness_result(i, 4)={m_id(3)};
            end
        % 相似度前三之基準樣本與測試樣本都不同
        else
            for k=1:3
                if (cell2mat(M(k, 2))>image_threshold)
                    fitness_result(i, k+1)={cell2mat(M(k, 1))};
                end
            end

        end
        
        
    end
    fitness=score/length(files);
end

%% 找尋對應鱗狀區塊
function [match_length, match_labels]=image_match(base_db, test_db, win)
    
    % 依 兩張影像基準線長度比 計算基準點延伸長度
    L_dist_p=double(test_db.D_n/base_db.D_n)*base_db.base.Distance;
    
    % 依 對應基準圖的角度及距離 對應出在比較圖上面的 對應點之X Y軸座標
    P_x = round(test_db.L_n(1,2)+(L_dist_p.*cosd(base_db.base.Theta).*base_db.base.Direct));
    P_y = round(test_db.L_n(2,2)-L_dist_p.*sind(base_db.base.Theta));   
    
    % 候選區域長度
    window=test_db.D_n*win;
    
    match_labels=struct([]);
    num=1;
    for i=1:length(base_db.base.Label)
        [shortlist, F]=region_shortlist(test_db.L, window, P_x(i), P_y(i), test_db.region);        
        if (isempty(shortlist)~=1)
            idx=find(base_db.region.Label==base_db.base.Label(i));
            F_base=base_db.region.Feature(idx, :);

            match_labels(num).label=base_db.base.Label(i);
            match_labels(num).shortlist=shortlist;
            match_labels(num).F_base=F_base;
            match_labels(num).F=F;
            num=num+1;

        else
            continue;
        end

    end
    match_length=length(base_db.base.Label);
end

%% 求候選區塊編號列表
function [shortlist, F]=region_shortlist(L, w, P_x, P_y, region)
    [row, column]=size(L);
        
        lowW = round(P_x-w);
        highW = round(P_x+w);
        lowH = round(P_y-w);
        highH = round(P_y+w);

        if lowH < 1
            lowH = 1;
        end
        if highH > row
            highH = row;
        end
        if lowW < 1
            lowW = 1;
        end
        if highW > column
            highW = column;
        end
        tempBlock = L(lowH:highH, lowW:highW);
        
        u = unique(tempBlock);
        u = u(u~=0 & u~=1);
        
        shortlist=[];
        F=[];
        if (isempty(u)~=1)
            l=1;
            for i=1:length(u)
                idx=find(region.Label==u(i));
                if (isempty(idx)~=1)
                    shortlist(l)=u(i);
                    F(l, :)=region.Feature(idx, :);
                    l=l+1;
                end
            end
        end

end

%% 找尋相對鱗狀區塊編號列表
function [match]=get_match_region(region, w, r, region_threshold)
    match=zeros(length(region), 3);
    for i=1:length(region)
        % 計算鱗狀區塊相似度
        [idx, score]=region_match(region(i).F_base, region(i).F, w, r, region_threshold);
        if (idx~=0) % 相似度最高之鱗狀區塊
            match(i, :)=[region(i).label, region(i).shortlist(idx), score];
        else        % 沒有相似的鱗狀區塊
            match(i, :)=[region(i).label, 0, score];
        end        
    end

end

%% 計算候補鱗狀區塊相似度
function [match, score]=region_match(F_base, F, w, r, threshold)
    i=1:length(F(:,1));
    j=1:length(F_base);
    t=w(j).*(abs(F_base(j)-F(i, j)).^r(j));
    m=sum(t, 2);
    [Min, I]=min(m);    % 取最相似(差異度最小)的鱗狀區塊

    match=I;
    score=Min;
    if (Min>threshold)  % 判斷是否低於閥值
        match=0;
        score=0;
    end
end
