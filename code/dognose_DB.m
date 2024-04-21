%% �إߪ���DB (�����v����ơB����I�B�S�x��)
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

        L_n = cell2mat( textscan(txt, '%f%f%f', 'Delimiter', ',', 'collectoutput', 1) ); % ��ǽu�y��

        D_n=points_distance(L_n(:, 1), L_n(:, 3)); % ��ǽu�Z��

        [~, column]=size(img);
        R_n=[1, column; 1, max(L_n(2, :))]; % ����d��(R_n)
        R_base=[L_n(1, 1), L_n(1, 3); (min(L_n(2, 1), L_n(2, 3))-(D_n*0.5)), max(L_n(2, 1), L_n(2, 3))]; % ��ǽd��(R_base)

        [L, num] = bwlabel(img, 4); % �N�C�Ӱ϶��s��
        [base, region]=region_extract(L, num, L_n, D_n, R_base, R_n); % ��ǽd��P����d���쪬�϶������


        id=strsplit(name_part_num, '_');
        i=str2double(id(1)); j=str2double(id(2));
        name=[i, j];
        
        info.name=name;         % �W��
        info.L=L;               % ��l�v��
        info.L_n=L_n;           % ��ǽu(�����դ����I�A����I�A�k���դ����I)
        info.D_n=D_n;           % ��ǽu����
        info.base=base;         % ��ǽd���쪬�϶�(�s���A�P����I�Z���A�P��ǽu�����A������I����[-1]�k[1]��)
        info.region=region;     % ����d���쪬�϶�(�s���A�S�x��[�����{��R_S�B���W�h�{��R_D�BMoment�S�x�B���פ�V���R_0, R_45�BR_90�BR_135]�A�����I�y��)
        
        new=['database_original\' dir_p '\'];
        if ~exist(new, 'dir')
            mkdir(new);
        end

        save([new 'dognose_' dir_p '_db-' name_part_num '.mat'], '-struct', 'info');

        fclose('all');
    end
end

%%  �D�ؼнd��(R)���U�ӽs���϶������I(C)X, Y�y�СB�P����I�Z���BL_n�PL_dist��������
function [base, region]=region_extract(L, num, L_n, D_n, R_base, R_n)
    
    base=cell2table(cell(0, 4));
    base.Properties.VariableNames={'Label', 'Distance', 'Theta', 'Direct'};
    region=cell2table(cell(0, 3));
    region.Properties.VariableNames={'Label', 'Center', 'Feature'};
    
    num_1=1; num_2=1;
    for i=2:num
        [p,q]=find(L==i);	% ���s��i ���Ҧ��y��[y, x]
        n=numel(p);
        x=sum(q)/n;	% �s��i �����IX�b�y�� (�s���ϰ�Ҧ��IX�ȥ[�`/�ϰ��`�Ӽ�)
        y=sum(p)/n;	% �s��i �����IY�b�y�� (�s���ϰ�Ҧ��IY�ȥ[�`/�ϰ��`�Ӽ�)

        if (R_base(1, 1)<=x && R_base(1, 2)>=x && R_base(2, 1)<=y && R_base(2, 2)>=y)
            dist=points_distance(L_n(:, 2), [x, y]);       % L_dist�ϰ줤���I�P��ǽu�����I�s�u�Z��
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

        if (R_n(1, 1)<=x && R_n(1, 2)>=x && R_n(2, 1)<=y && R_n(2, 2)>=y) % �P�O�϶��O�_�b����d��
            dist=points_distance(L_n(:, 2), [x, y]);       % L_dist�ϰ줤���I�P��ǽu�����I�s�u�Z��
            [R_S, R_D, Moment, R_0, R_45, R_90, R_135]=region_feature(L, i, D_n, dist);            
            
            region(num_2, :)=table(i, [x, y], [R_S, R_D, Moment, R_0, R_45, R_90, R_135]);
            
            num_2=num_2+1;
        end
    end
end

%% �D�϶��S�x��: �����{��R_S�B���W�h�{��R_D�BMoment�S�x�B���פ�V���R_0, R_45�BR_90�BR_135
function [R_S, R_D, Moment, R_0, R_45, R_90, R_135]=region_feature(L, label, D_n, L_dist)

    L(L~=label)=0;
    BW = imbinarize(L); % �ؼа϶��d��

    item=regionprops(BW, 'Area', 'Perimeter');                  % �Ӱ϶������n�B�P��
    [D_c, D, slope]=greatest_diameter(BW);                      % �϶��̪��Z���������I�y�СB���סB�ײv

    R_S=(4*item.Area)/(pi*(D^2));                               % �����{��R_S
    R_D=(4*pi*item.Area)/(item.Perimeter^2);                    % ���W�h�{��R_D
    Moment=sqrt(item.Area*L_dist)/(D_n*2);                      % Moment�S�x�A�Ψӥ��W�ƨC�@�϶�

    [D_0, D_45, D_90, D_135]=degree_dist(slope, BW, D_c, D);    % ���פ�V����D_0�BD_45�BD_90�BD_135

    if abs(D-D_0)<1
        D=D_0;
    end

    R_0=D_0/D;          % ���פ�V���R_0
    R_45=D_45/D;        % ���פ�V���R_45
    R_90=D_90/D;        % ���פ�V���R_90
    R_135=D_135/D;      % ���פ�V���R_135
    
end

%% �D�϶��̻��Z��
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

%% �D�������u�q�Z��(���]�t�¦ⳡ��)
function [D_0, D_45, D_90, D_135]=degree_dist(slope, BW, D_c, D)
    
    angle = atan(slope);
    th = 0:pi/4:2*pi;   %�H�̪��Z���������I������ 0, 45, 90, 135 180 225 270 315 �K�Ӥ�V���ײv

    % �p��U��쪺�Z��(���]�t�϶��~�������A�����O�p�G�����W�A�W�X�϶��~������)
    D_0=degree_dist_count((th(1)+angle), BW, D_c, D)+degree_dist_count((th(5)+angle), BW, D_c, D);      % 0, 180
    D_45=degree_dist_count((th(2)+angle), BW, D_c, D)+degree_dist_count((th(6)+angle), BW, D_c, D);    % 45, 225
    D_90=degree_dist_count((th(3)+angle), BW, D_c, D)+degree_dist_count((th(7)+angle), BW, D_c, D);    % 90, 270
    D_135=degree_dist_count((th(4)+angle), BW, D_c, D)+degree_dist_count((th(8)+angle), BW, D_c, D);   % 135, 315

end

%% �D�����I������D���u�q����(���]�t�¦ⳡ��)
function [dist]=degree_dist_count(degree, BW, D_c, D)

    % �H���ߨ����׬�D/2�A������degree����V�A�[�\���Ҧ��y��(���Z�� D*0.5)
    d=0:0.005:(D*0.5);
    x_p = d .* cos(degree) + D_c(1, 2);
    y_p = d .* sin(degree) + D_c(2, 2);

    [row, column]=size(BW);
    idx=1;
    for i=1:length(d)
        x=round(x_p(i)); y=round(y_p(i));
        if (x>=1 && x<=column) && (y>=1 && y<=row)
            line(idx)=BW(round(y), round(x)); %�������ײ[�\���y�� �� (�p�G�b�d�򤺬�1 �d��~��0)
            idx=idx+1;
        end
    end
    line_temp=[line, 0];
    line_temp_1=[0, line];
    edge=line_temp-line_temp_1; % ��M����I
    [start_p]=find(edge==1); % ��M�}�l�s�u���a�� (�p�G�����W�A�i�ण�u�@��)
    [end_p]=find(edge==-1); % ��M�������a��
    dist=0;

    if (isempty(start_p)==0)
        for j=1:length(start_p)
            dist=dist+points_distance([x_p(start_p(j)), y_p(start_p(j))], [x_p(end_p(j)-1), y_p(end_p(j)-1)]);
        end
    end
    
end

%% �D�y��p1 p2 ���s�u�Z��
function [dist]=points_distance(p1, p2)
    dist=sqrt(((p1(1)-p2(1))^2)+((p1(2)-p2(2))^2));
end
