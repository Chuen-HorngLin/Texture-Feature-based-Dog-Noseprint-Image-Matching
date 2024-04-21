%% �إ��쪬�϶��S�x����ƶ��A�Ω�GA��̨ΰѼ�
%% owner: yang shu chun
%%
clear all;
warning('off','all')
group='group3';
DB_path=['database_original\' group '\'];

win=0.025;

base_files=dir([DB_path 'base\*.mat']);

files=dir([DB_path 'test\*.mat']);

for i=1:length(files)
    temp=strrep(files(i).name, 'dognose_test_db-', '');
    name_part_num=strrep(temp, '.mat', '');
    
    test_db=load([DB_path 'test\' files(i).name]);

    match_db=struct([]);
    parfor j=1:length(base_files)
        base_db=load([DB_path 'base\' base_files(j).name]);
        [match_length_a, match_base]=image_match(test_db, base_db, win);    % ���ռ˥�����ǽd���쪬�϶��ƶq�B��X���ռ˥��������Ǽ˥����Ը��쪬�϶��S�x 
        [match_length_b, match_test]=image_match(base_db, test_db, win);    % ��Ǽ˥�����ǽd���쪬�϶��ƶq�B��X��Ǽ˥���������ռ˥����Ը��쪬�϶��S�x
        match_length=2*(1/((1/match_length_a)+(1/match_length_b)));         % ��ǽd���쪬�϶��ƶq���թM�Ѽ�
        
        match_db(j).baseID=base_db.name(1);
        match_db(j).match_base=match_base;
        match_db(j).match_test=match_test;
        match_db(j).match_length=match_length;

        fclose('all');
    end
    
    info.name=test_db.name;
    info.match_db=match_db;
    
    match_path=['match_db\' group '\'];
    if ~exist(match_path, 'dir')
        mkdir(match_path);
    end
    save([match_path 'match_db-' name_part_num '.mat'], '-struct', 'info');
end

function [match_length, match_labels]=image_match(base_db, test_db, win)

    L_dist_p=double(test_db.D_n/base_db.D_n)*base_db.base.Distance;
    
    %�H ������ǹϪ����פζZ�� �����X�b����ϤW���� �����I��X Y�b�y��
    P_x = round(test_db.L_n(1,2)+(L_dist_p.*cosd(base_db.base.Theta).*base_db.base.Direct));
    P_y = round(test_db.L_n(2,2)-L_dist_p.*sind(base_db.base.Theta));   
    
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

%% �D�Ըɰ϶��s���C��
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
