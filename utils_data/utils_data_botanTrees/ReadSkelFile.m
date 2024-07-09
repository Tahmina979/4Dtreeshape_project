function [branch, branch_num] = ReadSkelFile( skelFile )


fp1 = fopen(skelFile, 'r');
total_levelNum = fscanf(fp1, '%d', [1,1]);

for i = 1:total_levelNum
    
    if i == 1
        layer_string = fscanf(fp1, '%s', [1,1]);
        layer_id = fscanf(fp1, '%d', [1,1]);                        %layer_id��ָ�ĵڼ���layer.
        branch_num(i) = fscanf(fp1, '%d', [1,1]);
       
        fatherBranch_id = fscanf(fp1, '%d', [1,1]);                       %��������һ���id
        fatherBranch_point_id = fscanf(fp1, '%d', [1,1]);                     %����֦���о���Ľڵ�
        branch_point_num = fscanf(fp1, '%d', [1,1]);                %��¼ÿ�����ӵĵ�ĸ���
        
        branch(i,1).father_branch_id = fatherBranch_id;
        branch(i,1).father_point_id = fatherBranch_point_id;
            
            for j = 1 : branch_point_num                     %����ÿ������branch��skeleton������ֵ�Ͷ�Ӧ�İ뾶
                branch(i,1).point(j).x = fscanf(fp1, '%f', [1,1]);
                branch(i,1).point(j).y = fscanf(fp1, '%f', [1,1]);
                branch(i,1).point(j).z = fscanf(fp1, '%f', [1,1]);
                branch(i,1).point(j).r = fscanf(fp1, '%f', [1,1]);
            end
    end
  
    if i~=1        
        layer_string = fscanf(fp1, '%s', [1,1]);
        layer_id = fscanf(fp1, '%d', [1,1]);                                    %layer_id��ָ�ĵڼ���layer.
        branch_num(i) = fscanf(fp1, '%d', [1,1]);

        for k = 1:branch_num(i)
            fatherBranch_id = fscanf(fp1, '%d', [1,1]);                              %��������һ���id
            fatherBranch_point_id = fscanf(fp1, '%d', [1,1]);                              %���branch����Ӧ�ĺ�����
            branch_point_num = fscanf(fp1,'%d',[1,1]);                        %��¼ÿ�����ӵĵ�ĸ���
            
            branch(i,k).father_branch_id = fatherBranch_id;
            branch(i,k).father_point_id = fatherBranch_point_id;
            
            for j = 1:branch_point_num                                       %����ÿ������branch��skeleton������ֵ�Ͷ�Ӧ�İ뾶
                branch(i,k).point(j).x = fscanf(fp1,'%f',[1,1]);
                branch(i,k).point(j).y = fscanf(fp1,'%f',[1,1]);
                branch(i,k).point(j).z = fscanf(fp1,'%f',[1,1]);
                branch(i,k).point(j).r = fscanf(fp1,'%f',[1,1]);
            end
        end      
    end
end

% ��branch_numת����������
branch_num = branch_num';
fclose(fp1);
    

end

