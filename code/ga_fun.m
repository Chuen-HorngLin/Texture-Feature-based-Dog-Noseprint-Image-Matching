%% GA_algo 基因演算法
%% owner: yang shu chun

function [chrom_generation] = ga_fun(times, var, best_num, best_chrom, mutate)
    chrom_length=sum(var);  % 染色體總長度
%% Step1：initial 初始化
    chrom_initial = zeros(10, chrom_length); 
    if (times==1)    % 如果為第一輪，則隨機產生染色體     
        for i=1:10     % 產生10條隨機染色體
            chrom_initial(i, :) = get_chrom_initial(chrom_length);
        end
    else                   % 如果非第一輪，則前幾條繼承best_chrom，剩下隨機產生
        chrom_initial(1:best_num, :) = best_chrom;
        for j=(best_num+1):10
            chrom_initial(j, :) = get_chrom_initial(chrom_length);
        end
    end
%% Step 2: mutation 突變
    chrom_mutation = zeros(10, chrom_length);
    for i=1:10      % 若是達到突變率則改變隨機個數位元，1轉0 0轉1
        chrom_mutation(i, :) = get_chrom_mutation(chrom_length, mutate, chrom_initial(i, :));
    end
%% Step3: crossover 交配
    chrom_crossover = zeros(10, chrom_length);
    i=1;
    while i<11      % 隨機選擇兩條染色體，把chrom_initial第A條取隨機半組變數組與chrom_mutate第B條另一半組變數組組合成新的染色體
        R1=randi(10);
        R2=randi(10);
        if (R1~=R2)
            R1_chrom = chrom_initial(R1, :);
            R2_chrom = chrom_initial(R2, :);
            chrom_crossover(i, :) = get_chrom_crossover(chrom_length, R1_chrom, R2_chrom);
            
            i=i+1;
        end
    end
%% 基因組產生
    chrom_generation=cat(1, chrom_initial, chrom_mutation, chrom_crossover);
end

%% initial 初始化
function [chrom_initial]=get_chrom_initial(chrom_length)
    chrom_initial = zeros(1, chrom_length);
    initial_count=randi([1, chrom_length], 1);  % 基因中為1的位置總個數
    initial_index=sort(randperm(chrom_length, initial_count));  % 基因為1的位置
    chrom_initial(initial_index)=1;
    
%     chrom_initial = check_nonzero(chrom_initial, var);  % 檢查每組變數組至少有一個1
    
end

%% mutation 突變
function [chrom_mutation]=get_chrom_mutation(chrom_length, mutate, chrom_initial)
    mask = zeros(1, chrom_length);  % 用來之後轉換0 1
    mutate_count=round(mutate*chrom_length);    % 突變的基因個數
    mutate_index=sort(randperm(chrom_length, mutate_count));    %突變的基因位置
    mask(mutate_index)=1;
    chrom_mutation = bitxor(chrom_initial, mask);   % 突變位置位元0轉1，1轉0
    
%     chrom_mutation = check_nonzero(chrom_mutation, var);  % 檢查每組變數組至少有一個1

end

%% crossover 交配
function [chrom_crossover]=get_chrom_crossover(chrom_length, R1_chrom, R2_chrom)
    chrom_crossover = R1_chrom;
%     cross_count=randi([1, chrom_length], 1);    % 交換的基因個數
%     cross_index=sort(randperm(chrom_length, cross_count)); % 交換的基因位置
    cross_index=sort(randperm(chrom_length, round(chrom_length*0.5))); % 交換的基因位置
    chrom_crossover(cross_index)=R2_chrom(cross_index);
    
%     chrom_crossover = check_nonzero(chrom_crossover, var);  % 檢查每組變數組至少有一個1

end