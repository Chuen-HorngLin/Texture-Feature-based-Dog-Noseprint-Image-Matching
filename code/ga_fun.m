%% GA_algo ��]�t��k
%% owner: yang shu chun

function [chrom_generation] = ga_fun(times, var, best_num, best_chrom, mutate)
    chrom_length=sum(var);  % �V�����`����
%% Step1�Ginitial ��l��
    chrom_initial = zeros(10, chrom_length); 
    if (times==1)    % �p�G���Ĥ@���A�h�H�����ͬV����     
        for i=1:10     % ����10���H���V����
            chrom_initial(i, :) = get_chrom_initial(chrom_length);
        end
    else                   % �p�G�D�Ĥ@���A�h�e�X���~��best_chrom�A�ѤU�H������
        chrom_initial(1:best_num, :) = best_chrom;
        for j=(best_num+1):10
            chrom_initial(j, :) = get_chrom_initial(chrom_length);
        end
    end
%% Step 2: mutation ����
    chrom_mutation = zeros(10, chrom_length);
    for i=1:10      % �Y�O�F����ܲv�h�����H���ӼƦ줸�A1��0 0��1
        chrom_mutation(i, :) = get_chrom_mutation(chrom_length, mutate, chrom_initial(i, :));
    end
%% Step3: crossover ��t
    chrom_crossover = zeros(10, chrom_length);
    i=1;
    while i<11      % �H����ܨ���V����A��chrom_initial��A�����H���b���ܼƲջPchrom_mutate��B���t�@�b���ܼƲղզX���s���V����
        R1=randi(10);
        R2=randi(10);
        if (R1~=R2)
            R1_chrom = chrom_initial(R1, :);
            R2_chrom = chrom_initial(R2, :);
            chrom_crossover(i, :) = get_chrom_crossover(chrom_length, R1_chrom, R2_chrom);
            
            i=i+1;
        end
    end
%% ��]�ղ���
    chrom_generation=cat(1, chrom_initial, chrom_mutation, chrom_crossover);
end

%% initial ��l��
function [chrom_initial]=get_chrom_initial(chrom_length)
    chrom_initial = zeros(1, chrom_length);
    initial_count=randi([1, chrom_length], 1);  % ��]����1����m�`�Ӽ�
    initial_index=sort(randperm(chrom_length, initial_count));  % ��]��1����m
    chrom_initial(initial_index)=1;
    
%     chrom_initial = check_nonzero(chrom_initial, var);  % �ˬd�C���ܼƲզܤ֦��@��1
    
end

%% mutation ����
function [chrom_mutation]=get_chrom_mutation(chrom_length, mutate, chrom_initial)
    mask = zeros(1, chrom_length);  % �ΨӤ����ഫ0 1
    mutate_count=round(mutate*chrom_length);    % ���ܪ���]�Ӽ�
    mutate_index=sort(randperm(chrom_length, mutate_count));    %���ܪ���]��m
    mask(mutate_index)=1;
    chrom_mutation = bitxor(chrom_initial, mask);   % ���ܦ�m�줸0��1�A1��0
    
%     chrom_mutation = check_nonzero(chrom_mutation, var);  % �ˬd�C���ܼƲզܤ֦��@��1

end

%% crossover ��t
function [chrom_crossover]=get_chrom_crossover(chrom_length, R1_chrom, R2_chrom)
    chrom_crossover = R1_chrom;
%     cross_count=randi([1, chrom_length], 1);    % �洫����]�Ӽ�
%     cross_index=sort(randperm(chrom_length, cross_count)); % �洫����]��m
    cross_index=sort(randperm(chrom_length, round(chrom_length*0.5))); % �洫����]��m
    chrom_crossover(cross_index)=R2_chrom(cross_index);
    
%     chrom_crossover = check_nonzero(chrom_crossover, var);  % �ˬd�C���ܼƲզܤ֦��@��1

end