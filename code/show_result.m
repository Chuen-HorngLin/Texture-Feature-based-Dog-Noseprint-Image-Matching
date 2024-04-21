%% 繪製鱗狀區塊比對影像
%% 
group='group3';
load(['match_result_' group '_train.mat']);
DB_path=['database_original\' group '\'];
output_path=['show_result\train\' group '\'];

fitness_result=train_fitness_result;

for i=1:length(fitness_result(:, 1))
    test_name=fitness_result{i, 1};
    M=fitness_result{i, 5};
    test_db=load([DB_path 'train\dognose_test_db-' test_name '.mat']);
    if (cell2mat(fitness_result(i, 2))~=test_db.name(1))
        idx=find(cell2mat(M(:, 1))==test_db.name(1));
        plot_id=[1, idx];
        show_result_1(output_path, DB_path, test_db, plot_id, M);
    else
        plot_id=1;
        show_result_1(output_path, DB_path, test_db, plot_id, M);
    end
    
end
%%
function show_result_1(path, DB_path, test_db, plot_id, M)

    for k=1:length(plot_id)
        match=cell2mat(M(plot_id(k), 3));
        b_n=cell2mat(M(plot_id(k), 1));
        file=dir([DB_path 'base\dognose_base_db-' num2str(b_n) '_*.mat']);
        base_db=load([DB_path 'base\' file(1).name]);

        h(1)=figure;
        imshow(test_db.L);
        for i=1:length(match(:, 1))
            idx=find(test_db.region.Label==match(i, 1));
            if (isempty(idx)~=1 && match(i, 2)~=-1)
                figure(1);
                [y, x]=find(test_db.L==test_db.region.Label(idx));
                hold on;
                plot(x, y, 'color', 'c');
                hold off;
                text(test_db.region.Center(idx, 1), test_db.region.Center(idx, 2), int2str(match(i, 1)), 'color', 'b');
            end
        end

        h(2)=figure;
        imshow(base_db.L);
        for i=1:length(match(:, 1))
            idx=find(base_db.region.Label==match(i, 2));
            if (isempty(idx)~=1 && match(i, 2)~=-1)
                figure(2);
                [y, x]=find(base_db.L==base_db.region.Label(idx));
                hold on;
                plot(x, y, 'color', 'g');
                hold off;
                text(base_db.region.Center(idx, 1), base_db.region.Center(idx, 2), int2str(match(i, 1)), 'color', 'r');
            end
        end

        if ~exist(path, 'dir')
            Folder does not exist so create it.
            mkdir(path);
        end

        path_name=[path 'match_result_' num2str(test_db.name(1)) '-' num2str(test_db.name(2)) '=' num2str(base_db.name(1)) '.fig'];
        savefig(h ,path_name);

        close(h);

        test_image=cat(3, test_db.L, test_db.L, test_db.L);
%         test_image_color=cat(3, test_db.L, test_db.L, test_db.L);
        for i=1:length(match(:, 1))
            idx=find(test_db.region.Label==match(i, 1));
            if (isempty(idx)~=1 && match(i, 2)~=-1)
                [y, x]=find(test_db.L==test_db.region.Label(idx));
                for p=1:length(y)
                    test_image(y(p), x(p), 1)=0;
                    test_image(y(p), x(p), 2)=255;
                    test_image(y(p), x(p), 3)=255;

%                     test_image_color(y(p), x(p), 1)=0;
%                     test_image_color(y(p), x(p), 2)=255;
%                     test_image_color(y(p), x(p), 3)=255;
                end
                test_image = insertText(test_image,[(test_db.region.Center(idx, 1)), (test_db.region.Center(idx, 2))],int2str(match(i, 1)),'BoxOpacity',1, 'BoxColor','cyan', 'AnchorPoint', 'Center', 'FontSize',8,'TextColor','black');
            end
        end



        base_image=cat(3, base_db.L, base_db.L, base_db.L);
%         base_image_color=cat(3, base_db.L, base_db.L, base_db.L);
        for i=1:length(match(:, 2))
            idx=find(base_db.region.Label==match(i, 2));
            if (isempty(idx)~=1 && match(i, 2)~=-1)
                [y, x]=find(base_db.L==base_db.region.Label(idx));
                for p=1:length(y)
                    base_image(y(p), x(p), 1)=0;
                    base_image(y(p), x(p), 2)=255;
                    base_image(y(p), x(p), 3)=0;

%                     base_image_color(y(p), x(p), 1)=0;
%                     base_image_color(y(p), x(p), 2)=255;
%                     base_image_color(y(p), x(p), 3)=0;
                end
                base_image = insertText(base_image,[(base_db.region.Center(idx, 1)), (base_db.region.Center(idx, 2))],int2str(match(i, 1)),'BoxOpacity',1, 'BoxColor','green', 'AnchorPoint', 'Center', 'FontSize',8,'TextColor','black');
            end
        end


        result_img=[path 'img\'];
        if ~exist(result_img, 'dir')
            mkdir(result_img);
        end
        result_img_color=[path 'img_color\'];
        if ~exist(result_img_color, 'dir')
            mkdir(result_img_color);
        end

        imwrite(test_image, [result_img 'match_result_' num2str(test_db.name(1)) '-' num2str(test_db.name(2)) '=' num2str(base_db.name(1)) '_test.bmp'])
        imwrite(base_image, [result_img 'match_result_' num2str(test_db.name(1)) '-' num2str(test_db.name(2)) '=' num2str(base_db.name(1)) '_base.bmp'])

%         imwrite(test_image_color, [result_img_color 'match_result_' num2str(test_db.name(1)) '-' num2str(test_db.name(2)) '=' num2str(base_db.name(1)) '_test.bmp'])
%         imwrite(base_image_color, [result_img_color 'match_result_' num2str(test_db.name(1)) '-' num2str(test_db.name(2)) '=' num2str(base_db.name(1)) '_base.bmp'])


    end
end