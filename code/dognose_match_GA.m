%% dognose_match_GA 狗鼻紋比對參數之基因演算法
%% owner: yang shu chun
clc; clear all;
%% 預設基本設定
t=datestr(now,'yyyy-mm-dd HH-MM');
group='group3';
ga_name=['dognose_match_GA_' group '_'];
filename = [ga_name t '.xlsx'];      % excel名稱
header = {'chrom', 'ith', 'w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'R', 'image_threshold', 'region_threshold', 'fitness'};  % 欄位名稱
var=[5,5,5,5,5,5,5, 5,5,5,5,5,5,5 ,3,3,3];                              % 變數位元數
best_num=8;                           % 最佳染色體保留數量
mutate=0.1;                           % 染色體變異機率

match_db_path=['C:\Users\User\Desktop\畢業論文\final\random_group\match_db\' group '\train\'];

db_start = datestr(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
disp(['Start set DB: ', db_start])

%% 基因找尋最佳變數組
chrom_length=sum(var);      % 染色體長度(變數位元數總和)
% best_chrom=[];
best_chrom=[];                     % 最佳染色體(初始為空)
final_chrom=[];                     % 每條染色體運算結果(最佳排到最差)
times=1;                                % 迭代次數
times_upper=300;                % 迭代次數上限
best_fitness=0;
best_count=0;
% while (times<=times_upper && best_count<=50)                 % 迭代次數超過上限 或是 fitness到達標準時停止 
while (times<=times_upper)                 % 迭代次數超過上限 或是 fitness到達標準時停止 
    tic     % 開始計時
    chrom_generation=ga_fun(times, var, best_num, best_chrom, mutate);  % 產生30條染色體
    chrom_ith=zeros(30, 1);
    score=zeros(30, 1);
    var_dec=zeros(30, length(var));
    parfor ith=1:30
%     for ith=1:30
        chrom=chrom_generation(ith, :);
        chrom_idx=ith+(times-1)*30;        
        chrom_ith(ith)=chrom_idx;    
        dec=bintodec(chrom, var);   % 將染色體從二進位轉成十進位
        [w, r, R, image_threshold, region_threshold]=get_var_value(dec);        
        var_dec(ith, :)=cat(2, w, r, R, image_threshold, region_threshold);
        score(ith)=get_fitness(match_db_path, w, r, image_threshold, region_threshold);    % fitness
    end
    times_chrom=cat(2, chrom_generation, chrom_ith, var_dec, score);
    [out,idx] = sortrows(times_chrom, (chrom_length+length(var)+2), 'descend'); %times_chrom依照fitness排序
    best_chrom=out(1:best_num, 1:chrom_length);    %最佳的幾條染色體儲存到best_chrom，用來下一輪染色體運算
    final_chrom=cat(1, out, final_chrom);   % 把每一輪染色體運算結果儲存在final_chrom
        
    if out(1, end)>best_fitness
        best_fitness=out(1, end);
        best_count=1;
    else
        best_count=best_count+1;
    end
    
    disp('=========================================');
    disp(['Time_' num2str(times) ' , ' 'best_fitness: ' num2str(out(1, end)) ' , ' 'best_count: ' num2str(best_count)]);
    
    save([ga_name t '.mat'], 'final_chrom', '-v7.3');
        
    times=times+1;
    
    toc
end
chrom_str=chrom2str(final_chrom(:, 1:chrom_length));
matrix2excel(filename, header, chrom_str, final_chrom(:,(chrom_length+1):end));

db_end = datestr(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
disp(['End set DB: ', db_end])

%% 變數二進位轉成十進位
function [var_dec]=bintodec(chrom, var)
    var_dec=zeros(1, length(var));
    index=1;
    for i=1:length(var)  % 依照每個變數的位元數轉換成十進位數字
        bin=cellstr(num2str(chrom(index : (index+var(i)-1)), '%d'));
        var_dec(i)=bin2dec(bin);
        index=index+var(i);
    end
end

%% 將十位數轉換成變數值
function [w, r, R, image_threshold, region_threshold]=get_var_value(var)
    
    w=zeros(1, 7);
    r=zeros(1, 7);
    R=7+var(15);
    for i=1:7
        w(i)=var(i)/sum(var(1:7));
        r(i)=(var(i+7)/sum(var(8:14)))*R*0.5;
    end
    
    image_threshold=0.3+(var(16)*0.025);
    
    region_threshold=0.2+(var(17)*0.025);
    
%     win=(var(18)+1)*0.0125;
    
end

%% chromosome轉字串
function [chrom]=chrom2str(arr)
    [row, ~]=size(arr);
    chrom=cell(row, 1);
    for i=1:row
        chrom(i)=cellstr(mat2str(arr(i, :)));
    end
end

%% 運算結果儲存到excel
function matrix2excel(filename, header, chrom, data)
    [row, column]=size(data);
    xls=cell(row+1, column+1);
    xls(1, :)=header;    
    xls(2:end, 1)=chrom;
    xls(2:end, 2:end)=num2cell(data);
    xlswrite(filename,xls);
end

%% 計算fitness
function [fitness]=get_fitness(path, w, r, image_threshold, region_threshold)

    files=dir([path '*.mat']);
    score=0;
%     fitness_result=cell(length(files), 3);
    for i=1:length(files)
        match_db=[];
        load([path files(i).name], 'name', 'match_db');
        acc=cell(length(match_db), 3);
        parfor j=1:length(match_db)
            match_base=get_match_region(match_db(j).match_base, w, r, region_threshold);
            match_test=get_match_region(match_db(j).match_test, w, r, region_threshold);
            s=0;
            for m=1:length(match_base(:, 1))
                index=find(match_test(:, 1)==match_base(m, 2));
                if (isempty(index)~=1 && length(index)==1)
                    if (match_test(index, 2)==match_base(m, 1))
                        s=s+1;
                    else
                        match_base(m, 2)=-1;
                    end
                else
                    match_base(m, 2)=-1;
                end
            end
            match=s/(match_db(j).match_length); % 兩張圖相似程度        
            acc(j, :)={match_db(j).baseID, match, match_base};            
        end
        
        [M]=sortrows(acc, 2, 'descend');
%         match_idx=cell2mat(M(1, 1));
        if (cell2mat(M(1, 2))>image_threshold)
            if (cell2mat(M(1, 1))==name(1))
                score=score+1;                
            end
%         else
%             match_idx=-1;
        end
        
%         fitness_result(i, 1)={name};
%         fitness_result(i, 2)={match_idx};
%         fitness_result(i, 3)={M};
        
        clear match_db;
        fclose('all');
    end
    fitness=score/length(files);
end

%% 取得對應的鱗狀區塊
function [match]=get_match_region(region, w, r, region_threshold)
    match=zeros(length(region), 3);
    for i=1:length(region)
        [idx, score]=region_match(region(i).F_base, region(i).F, w, r, region_threshold);
        if (idx~=0)
            match(i, :)=[region(i).label, region(i).shortlist(idx), score];
        else
            match(i, :)=[region(i).label, 0, score];
        end        
    end

end

%% 計算鱗狀區塊候補的相似度
function [match, score]=region_match(F_base, F, w, r, threshold)
    i=1:length(F(:,1));
    j=1:length(F_base);
    %     d=abs(F_base(j)-F(i, j));
    t=w(j).*(abs(F_base(j)-F(i, j)).^r(j));
    m=sum(t, 2);
    [Min, I]=min(m);

    match=I;
    score=Min;
    if (Min>threshold)
        match=0;
        score=0;
    end
end