%%
% ****************************************************************
% *                                                                                             *
% *                          Finding outlets                                           *
% *                                                                                             *
% *                          Connection Checking                                 *
% *                                                                                             *
% *                            node and conduit processing                  *
% *                                                                                             *
% ****************************************************************


%% 读取文件并指定字段

clear all
watershed = shaperead('土地利用_edit01.shp');
tubeNodes=shaperead('检查井_edit01.shp');
line=shaperead('管网_edit01.shp');    % 包含长度
% outlet = shaperead('排水口_edit01.shp');
watershedMeanEle=importdata('meanelev_subcat_edit01.txt'); 
% watershedMeanEle_array = watershedMeanEle.data;
    %第二列是id（FID）；第五列是平均高程
id_col = 2 ;
elev_col = 5;
%centerPoint=shaperead('point.shp'); %watershed中X1和Y1好像就是中心点的坐标。

% 指定nodatavalue
nodataValue = 60;

% 统计数量
numW=length(watershed); % 共有numW个子汇水区
numT=length(tubeNodes); % 共有numT个管点


% 土地利用类型中分类的字段 -------------------------------------------
buildings_string = '房屋';    % 需要单独处理的类型字段
% 以遍历的方式统计上述字段（类型）中类型的数量，并附新字段type分类
underlyting_type_info_struct = struct('type_char',[],'type_num',[]);
for i=1:numW
    temp_char=watershed(i).x0xE70xB10xBB0xE50x9E0x8B;
    count=sum(find(strcmp({underlyting_type_info_struct.type_char},temp_char)));
    if count==0 %类型未存
        size_underlyting_type_info_struct=length([underlyting_type_info_struct.type_num]);
        underlyting_type_info_struct(size_underlyting_type_info_struct+1).type_char=temp_char;
        underlyting_type_info_struct(size_underlyting_type_info_struct+1).type_num=size_underlyting_type_info_struct+1;
    end

    % 储存类型对应的编号
    watershed(i).type_num=...
        underlyting_type_info_struct(find(strcmp({underlyting_type_info_struct.type_char},temp_char))).type_num;
    % newLine_river(num_riverlinks).fromNode = ...
    %     newPoint_river(find(strcmp({newPoint_river.node_id},newLine_river(num_riverlinks).us_node_id))).new_id;
    
    % 计算多边形中心
    watershed(i).centerX = (max(watershed(i).X) + min(watershed(i).X)) / 2;
    watershed(i).centerY = (max(watershed(i).Y) + min(watershed(i).Y)) / 2;

    % 将区域的平均面积赋值到watershed结构体中
    index = find(watershedMeanEle.data(:,id_col)==(i-1));    % 查找对应的watershed的高程统计是否存在在watershedMeanEle中
    if length(index) > 0
        watershed(i).mean_elev = watershedMeanEle.data(index,elev_col);
    elseif length(index) == 0
        watershed(i).mean_elev = nodataValue;
    end

    % 对子汇水区进行编号
    watershed(i).id = i;

end

buildings_typeNum=underlyting_type_info_struct(find(strcmp({underlyting_type_info_struct.type_char},buildings_string))).type_num;

% 对建筑区域的平均高程添加20m
for i=1:numW
    if watershed(i).type_num == buildings_typeNum
        watershed(i).mean_elev = watershed(i).mean_elev+20;
    end
end



%% 感觉其实要先确定管线和点数据的连通性，再算下面的内容；如果存在需要剔除的点，就需要修正或者删除

clear shp_points
clear shp_link

% 这个要不还是用原始编号吧，如果出问题了再说
shp_pipe_points = tubeNodes;
shp_pipe_link = line;

% ①处理newPoint数据
for num_points=1:length(shp_pipe_points)
    shp_points(num_points).Geometry = shp_pipe_points(num_points).Geometry;
    shp_points(num_points).X = shp_pipe_points(num_points).X;
    shp_points(num_points).Y = shp_pipe_points(num_points).Y;
    shp_points(num_points).monitor = shp_pipe_points(num_points).x0xE70x9B0x910xE60xB50x8B;
    shp_points(num_points).id = shp_pipe_points(num_points).x0xE70xBC0x960xE50x8F0xB7;
    shp_points(num_points).topElev = shp_pipe_points(num_points).x0xE90xA10xB60xE60x870xE90xAB0x98; 
    shp_points(num_points).bottomElev = shp_pipe_points(num_points).x0xE50xBA0x950xE60x870xE90xAB0x98;

    % newPoint(num_points).chambfloor = shp_pipe_points(num_points).chambfloor;
    % newPoint(num_points).chambroof = shp_pipe_points(num_points).chambroof;
    % newPoint(num_points).ground_lev = shp_pipe_points(num_points).ground_lev;
    % newPoint(num_points).node_type = shp_pipe_points(num_points).node_type;
    % % swmm中输入、需要计算出来的值
    % % 但不确定这个是不是对的……？
    % newPoint(num_points).maxdepth = newPoint(num_points).ground_lev-newPoint(num_points).chambfloor;
    % % 一些不是特别重要的设定值
    % newPoint(num_points).initDepth = 0;
    % newPoint(num_points).surdepth = 15;
    % newPoint(num_points).aponded = 1000;
    % 
    % % 本次制作point新数据的主要目的在这里
    % newPoint(num_points).node_id = shp_pipe_points(num_points).node_id;
    % newPoint(num_points).new_id = num_points+node_idx;

end



% line.
% ②处理newLine数据
for num_pipes=1:length(shp_pipe_link)
    shp_link(num_pipes).Geometry = shp_pipe_link(num_pipes).Geometry;
    shp_link(num_pipes).BoundingBox = shp_pipe_link(num_pipes).BoundingBox;
    shp_link(num_pipes).X = shp_pipe_link(num_pipes).X;
    shp_link(num_pipes).Y = shp_pipe_link(num_pipes).Y;
    shp_link(num_pipes).monitor = shp_pipe_link(num_pipes).x0xE70x9B0x910xE60xB50x8B;
    shp_link(num_pipes).id = shp_pipe_link(num_pipes).x0xE70xBC0x960xE50x8F0xB7;
    shp_link(num_pipes).fromNode = shp_pipe_link(num_pipes).x0xE40xB80x8A0xE60xB80xB80xE40xBA0x95;
    shp_link(num_pipes).toNode = shp_pipe_link(num_pipes).x0xE40xB80x8B0xE60xB80xB80xE40xBA0x95;
    shp_link(num_pipes).length = shp_pipe_link(num_pipes).x0xE90x950xBF0xE50xBA0xA6;
    shp_link(num_pipes).lineType = shp_pipe_link(num_pipes).x0xE70xAE0xA10xE70xB10xBB0xE50x9E0x8B;
    shp_link(num_pipes).width = shp_pipe_link(num_pipes).x0xE70x9B0xB40xE50xBE0x84_0xE50xAE0xBD;
    shp_link(num_pipes).height = shp_pipe_link(num_pipes).x0xE60xB80xE90xAB0x98;
    shp_link(num_pipes).in_height = shp_pipe_link(num_pipes).x0xE80xBF0x9B0xE70xAE0xA10xE50xBA0x95;
    shp_link(num_pipes).out_height = shp_pipe_link(num_pipes).x0xE50x870xBA0xE70xAE0xA10xE50xBA0x95;

end


%% B-01 线数据（管网河道）拓扑关系检查连通性检查
% ----------（content011）通过筛选条件取出流域出口汇水点（P）并编号-----------------------------------
%https://blog.csdn.net/u013555719/article/details/105547558
% [shp_points.Id]--取出的是数组
% {shp_points.type}--取出的是字符串元胞
% 应该还有别的转成数组（而不是元胞）的办法吧
% [shp_points.Id]==50 返回的是1×76 logical 数组
% ismember({shp_points.type},'link') 返回的是1×76 logical 数组
    % temp_findPoint=shp_points([shp_points.Id]==50);
    % temp_findPoint=shp_points(ismember({shp_points.type},'link'));
    % temp_findPoint=shp_points(([shp_points.Id]==50)&(ismember({shp_points.type},'link')));
%点的Id类型为数组，线的Id类型为字符串

clear newPoint
clear newLine

outlet_id='Out001';

% 不对，Out001是不存在于shp_points，最后也不会存在在newPoint中的点
% 所以先用随便哪个点替代，最后删掉
newPoint(1)=shp_points(1);
newPoint(1).id=outlet_id;
% newPoint(1)=...
%     shp_points(find(strcmp({shp_points.id},outlet_id)));
newPoint(1).markPoint=0;

%初始化循环的要素
k=1;
judge_loop=1;

% ----------（content011）--------------------------------------------------------------------------------------
while judge_loop>0 %(如果newPoint中还有markPoint==0的节点就接着循环)
    % ----------（content012）根据出口节点P编号/名称，搜索下游节点为P的（管道或河道），存到新数组-----------------
        % isempty(fieldnames(tempLine_new))  判断结构体是否为空
        %如果结构体为空的时候，不能用常规方法在数组后面添加新对象
    clear tempLine
    downPointName=newPoint(k).id;   % 字符串
    tempLine=shp_link(find(strcmp({shp_link.toNode},downPointName)));



    %检查tempLine中是否存在与newLine中已储存的项目重复的项目
    Id_in_tempLine={tempLine.id}; % 取出的是元胞
    num_of_tempLine=size(Id_in_tempLine,2);
    clear tempLine_new
    %tempLine_new=struct();
    if num_of_tempLine>0 %有的端点上游已经没有线和点了，直接跳过这个取出对应上游线和对应点的步骤
        for i=1:num_of_tempLine
            %第一次向newLine存储线元素的情况（单独处理）-------------
            if exist('newLine')==0
                if exist('tempLine_new')==0
                    %如果结构体为空的时候，不能用常规方法在数组后面添加新对象
                    %https://blog.csdn.net/u010247905/article/details/51356797
                    tempLine_new(1)=tempLine(i);
                else
                    %结构体不为空时可以用常规方法添加新对象
                    tempLine_new=[tempLine_new,tempLine(i)];
                end
                continue;
            end
            %--------------------------------------------------------------------
            % if sum(ismember({newLine.Id},Id_in_tempLine{i}))>0
            % if sum(ismember([newLine.new_id],Id_in_tempLine{i}))>0
            if sum(ismember({newLine.id},Id_in_tempLine{i}))>0   % 换成字符串以后不知道成不成立
                continue
            end
            if exist('tempLine_new')==0
                %如果结构体为空的时候，不能用常规方法在数组后面添加新对象
                %https://blog.csdn.net/u010247905/article/details/51356797
                tempLine_new(1)=tempLine(i);
            else
                %结构体不为空时可以用常规方法添加新对象
                tempLine_new=[tempLine_new,tempLine(i)];
            end
        end
        clear tempLine
        tempLine=tempLine_new;
        if exist('newLine')==0
            newLine=tempLine;
        else
            newLine=[newLine,tempLine];
        end

        % ----------（content012）----------------------------------------------------------------------------------------------

        % ----------（content013）提取所有新线元素对应的上游节点，存到临时数组upperPoints(i)，i = from 1 to n------------
        % 操作应该和上面相似，这个更需要判断和去重
        upperPoint_IdArray={tempLine.fromNode}; % 现在是字符串了
        %根据id取出所有的上游节点
        clear upperPoint_struct_temp
        for i=1:size(upperPoint_IdArray,2)
            %查找up id对应的point在shp_point的数据
            clear temp_point
            temp_point=shp_points(find(strcmp({shp_points.id},upperPoint_IdArray(i))));


            if exist('upperPoint_struct_temp')==0
                upperPoint_struct_temp=temp_point;
            else
                upperPoint_struct_temp=[upperPoint_struct_temp,temp_point];
            end
        end

         
        %逐个判断取出的upperPoint_struct_temp中与newPoint是否存在重复
        for i=1:size(upperPoint_struct_temp,2)
            judge_exist=length(find(strcmp({newPoint.id},upperPoint_struct_temp(i).id)));
            if judge_exist>0 %表示存在重复
                continue;
            else %judge_exist==0，不重复
                %newPoint有流域出流点了，不是空的
                %error-20230303001：upperPoint_struct_temp(i)没有markPoint的参数和属性，不能直接添加
                upperPoint_i_temp=upperPoint_struct_temp(i);
                upperPoint_i_temp(1).markPoint=0;
                newPoint=[newPoint,upperPoint_i_temp];
            end
        end
        % ----------（content013）-------------------------------------------------------------------------------------------------------
    end
    % 选取下一个没有遍历的Point，and mark the previous one as '1'
    newPoint(k).markPoint=1;
    k=k+1; %逐个遍历newPoint里的元素点
    judge_loop=sum([newPoint.markPoint]==0);
    
end



% find(strcmp({newPoint.id},'Y030799'))





% ----------（content014）-------------------------------------------------------------------------------------------------
%
%绘制shp图
% https://blog.csdn.net/weixin_44058353/article/details/123241109
% https://blog.csdn.net/slandarer/article/details/123098323

% Lat_min = min([newPoint.X]);
% Lat_max = max([newPoint.X]);
% Lon_min = min([newPoint.Y]);
% Lon_max = max([newPoint.Y]);
% worldmap([Lat_min Lat_max],[Lon_min Lon_max]); 
newLine_draw=newLine';
fig2=figure(2);
fig2.Position=[1 1 800 800];
% mapshow(saveWatershed,'FaceColor','#f2debd');
% alpha(0.4)
mapshow(newLine_draw);
hp2(1)=mapshow(newLine_draw,'LineWidth',1.5);

xlabel('X坐标','fontsize',15,'FontWeight','bold')
ylabel('Y坐标','fontsize',15,'FontWeight','bold')
legend(hp2,'管线',...
    'Fontsize',12,'Location','northwest')



%% 

% % 计算平均管线长度，作为搜索半径的依据
% mean_line_length = mean([line.x0xE90x950xBF0xE50xBA0xA6]);

bar = waitbar(0,'计算到各子汇水区距离...');

for n = 1:numW
    %取节点n的中心
    Xc = watershed(n).centerX;
    Yc = watershed(n).centerY;
    
    %--PT1--------------------------------------------------------------------------------------------
    % 中心点到节点（可以直接计算距离，取其最小值并暂时保存）
    % 分别计算各中心点到各管网节点的距离
    distanceT=nan(1,numT);
    for nodes = 1:numT
        X=tubeNodes(nodes).X;
        Y=tubeNodes(nodes).Y;
        distanceT(nodes)=sqrt((X-Xc)^2+(Y-Yc)^2);
    end
    clear X
    clear Y
    minDistanceTube=min(distanceT);
    minDistanceTubeIndex=find(distanceT==minDistanceTube); %最小值在tubeNodes中的位列
    
    

    %--PT2--------------------------------------------------------------------------------------------
    % 中心点到子汇水区
    % 1 - 分别计算各中心点i到其余各中心点的距离
    % 2 - 对距离进行排序，并记录对应的子汇水区编号(可略，后续用min并用nan替代已遍历值即可)
    % 3 - 先取中心距离最短的子汇水区k，计算其中所有边界转折点到中心i的距离
    % 4 - 取距离最短的两个点，算中心点到线段的最短距离
    % 5 - 比较k中心和中心i的高程平均值。
        %如果elevation k>elevation i，则记录k为i子汇水区水流可能流入的邻子汇水区；
    % 6 - 取下一个距离最短的k'，重复步骤3-5
        %如果存在elevation k > elevation i（已记录）
        %取记录中点到线段距离最小的邻子汇水区作为子汇水区可能值。
    


    %1
    distanceCenter=nan(1,numW);
    for c=1:numW
        distanceCenter(c)=sqrt((Xc-watershed(c).centerX)^2+(Yc-watershed(c).centerY)^2);
    end
    distanceCenter(1,n)=nan;
    
    minDistanceCenter=min(distanceCenter);
    minDistanceCenterIndex=find(distanceCenter==minDistanceCenter);
    watershedMeanEle_array = [watershed.mean_elev];
    MeanEle_Center=watershedMeanEle_array(n);
    if MeanEle_Center == nodataValue
        MeanEle_Center = 9999;
    end
    count=0;
    saveID=[]; %watershed.ID
    saveIndex=[]; %watershed Index
    saveDis=[]; %the distance
    
    while count<numW

        % 比对中心高程，如果不对就跳过取下一个
        MeanEle_k=watershedMeanEle_array(minDistanceCenterIndex);
        if MeanEle_k<MeanEle_Center
            %3 -取多边形周围的坐标
            x=watershed(minDistanceCenterIndex).X;
            y=watershed(minDistanceCenterIndex).Y;
            distanceBox=sqrt((Xc-x).^2+(Yc-y).^2);
            %4 -取距离最近的两个点，计算点到线段的最短距离
            indexA=find(distanceBox==min(distanceBox));
            if size(indexA,2)>1
                indexA=indexA(1);
            end
            A=[x(indexA),y(indexA)];
            distanceBox(1,indexA)=nan;
            indexB=find(distanceBox==min(distanceBox));
            if size(indexB,2)>1
                indexB=indexB(1);
            end
            B=[x(indexB),y(indexB)];
            distanceBox(1,indexB)=nan;
                        if B(1)==A(1) && B(2)==A(2)
                            indexB=find(distanceBox==min(distanceBox));
                            if size(indexB)>1
                                indexB=indexB(1);
                            end
                            B=[x(indexB),y(indexB)];
                            distanceBox(1,indexB)=nan;
                        end
            
            longAB=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);
            P=[Xc Yc];
            AP=[(P(1)-A(1)),(P(2)-A(2))];
            AB=[(B(1)-A(1)),(B(2)-A(2))];
            longAP=sqrt((A(1)-P(1))^2+(A(2)-P(2))^2);
            r=(dot(AP,AB))/(longAB)^2;   
            if r>=1
                d=sqrt((P(1)-B(1))^2+(P(2)-B(2))^2);
            elseif r<=0
                d=sqrt((P(1)-A(1))^2+(P(2)-A(2))^2);
            else
                %AB直线一般式的系数
                a=B(2)-A(2);
                b=A(1)-B(1);
                c=B(1)*A(2)-A(1)*B(2);
                d=abs(a*P(1)+b*P(2)+c)/sqrt(a^2+b^2);
            end
            count=count+1;
            saveID=[saveID,watershed(minDistanceCenterIndex).id]; %watershed.ID
            saveIndex=[saveIndex,minDistanceCenterIndex]; %watershed Index
            saveDis=[saveDis,d]; %the distance
            %将ID取出并将准备搜索的子汇水区换成下一个
            distanceCenter(distanceCenter==minDistanceCenter)=nan;
            minDistanceCenter=min(distanceCenter);
            minDistanceCenterIndex=find(distanceCenter==minDistanceCenter);

        else
            count=count+1;
            %将ID取出并将准备搜索的子汇水区换成下一个
            distanceCenter(distanceCenter==minDistanceCenter)=nan;
            minDistanceCenter=min(distanceCenter);
            minDistanceCenterIndex=find(distanceCenter==minDistanceCenter);
        end
        
    end
    
        
    %--PT3--------------------------------------------------------------------------------------------
    %比较到点最小值和到多边形最小值，取更小值和对应要素作为该子汇水区下一个流向的区域。
    minDistancePolygon=min(saveDis);
    minDistancePolygonIndex=find(saveDis==minDistancePolygon);
    if length(minDistancePolygonIndex)>1
        minDistancePolygonIndex = minDistancePolygonIndex(1,1);  % 不止一个就删掉
    end
    minDistancePolygonID=saveID(minDistancePolygonIndex);
    if size(minDistancePolygon,1)==0
        minDistancePolygon=99999;
    end
        % 但如果这是个建筑那就直接流到排水点
    if watershed(n).type_num == buildings_typeNum
        minDistancePolygon=99999;
    end

    
    if minDistancePolygon>=minDistanceTube
        %流向节点并优先取节点
        indexOfWatershed = find([watershed.id]==watershed(n).id);
        watershed(indexOfWatershed).IDofWatershedOrNodes = num2str(tubeNodes(minDistanceTubeIndex).x0xE70xBC0x960xE50x8F0xB7);
        watershed(indexOfWatershed).outType = 'node';
        watershed(indexOfWatershed).outDistance = minDistanceTube;
        % 
        % flowTo(n).watershedID=watershed(n).id;
        % flowTo(n).IDofWatershedOrNodes=tubeNodes(minDistanceTubeIndex).x0xE70xBC0x960xE50x8F0xB7;
        % flowTo(n).Type='node'; 
        % flowTo(n).Distance=minDistanceTube;
    else
        %流向子汇水区
        indexOfWatershed = find([watershed.id]==watershed(n).id);
        watershed(indexOfWatershed).IDofWatershedOrNodes = minDistancePolygonID;
        watershed(indexOfWatershed).outType = 'subcatchment';
        watershed(indexOfWatershed).outDistance = minDistancePolygon;

        % flowTo(n).watershedID=watershed(n).id;
        % flowTo(n).IDofWatershedOrNodes=minDistancePolygonID;
        % flowTo(n).Type='subcatchment';
        % flowTo(n).Distance=minDistancePolygon;
    end

    str=['寻找子汇水区出流去向...',num2str(100*n/numW),'%'];
    waitbar(n/numW,bar,str)                       % 更新进度条bar，配合bar使用


end
close(bar)




%% B-02 添加面要素的拓扑关系检查连通性检查
saveWatershed=watershed;
for i=1:numW
    saveWatershed(i).outletID=watershed(i).IDofWatershedOrNodes;
    saveWatershed(i).outletType=watershed(i).outType;
end


%
% 为每个子汇水区添加虚拟管道
% 以便将子汇水区算入连通性分析当中

% step01- 添加虚拟管
% step02- 将子汇水区的编号作为信息存储到shp_points中
virtualLink_initialID=90000;
num_originLine=length(newLine);
num_originPoint=length(newPoint);
num_watershed=length(saveWatershed);

shp_pointsAndSubca=newPoint;
shp_linkAndSubca=newLine;

for k=1:num_watershed
    shp_linkAndSubca(k+num_originLine).id=num2str(virtualLink_initialID+k);
    shp_linkAndSubca(k+num_originLine).fromNode=num2str(saveWatershed(k).id);
    shp_linkAndSubca(k+num_originLine).toNode=num2str(saveWatershed(k).outletID);
    shp_pointsAndSubca(k+num_originPoint).id=num2str(saveWatershed(k).id);
end
%%;
for k = 1:length(shp_pointsAndSubca)
    shp_pointsAndSubca(k).markPoint=0;
end




%
clear newPoint2
clear newLine2
clear newWatershed

newPoint2(1)=newPoint(1);
newPoint2(1).markPoint=0;

%初始化循环的要素
k=1;
judge_loop=1;
%newLine2=struct();
%test:
%newLine2(1)=shp_link(1);

% ----------（content011）--------------------------------------------------------------------------------------
while judge_loop>0 %(如果newPoint中还有markPoint==0的节点就接着循环)
    % ----------（content012）根据出口节点P编号/名称，搜索下游节点为P的（管道或河道），存到新数组-----------------
        % isempty(fieldnames(tempLine_new))  判断结构体是否为空
        %如果结构体为空的时候，不能用常规方法在数组后面添加新对象
    clear tempLine
    downPointName=newPoint2(k).id;
    tempLine=shp_linkAndSubca(find(strcmp({shp_linkAndSubca.toNode},downPointName)));


    %检查tempLine中是否存在与newLine中已储存的项目重复的项目
    Id_in_tempLine={tempLine.id}; % 取出的是元胞
    num_of_tempLine=size(Id_in_tempLine,2);
    clear tempLine_new
    %tempLine_new=struct();
    if num_of_tempLine>0 %有的端点上游已经没有线和点了，直接跳过这个取出对应上游线和对应点的步骤
        for i=1:num_of_tempLine
            %第一次向newLine存储线元素的情况（单独处理）-------------
            if exist('newLine2')==0
                if exist('tempLine_new')==0
                    %如果结构体为空的时候，不能用常规方法在数组后面添加新对象
                    %https://blog.csdn.net/u010247905/article/details/51356797
                    tempLine_new(1)=tempLine(i);
                else
                    %结构体不为空时可以用常规方法添加新对象
                    tempLine_new=[tempLine_new,tempLine(i)];
                end
                continue;
            end
            %--------------------------------------------------------------------
%             if sum(ismember({newLine2.Id},Id_in_tempLine{i}))>0
            if sum(ismember({newLine2.id},Id_in_tempLine{i}))>0
                continue
            end
            if exist('tempLine_new')==0
                %如果结构体为空的时候，不能用常规方法在数组后面添加新对象
                %https://blog.csdn.net/u010247905/article/details/51356797
                tempLine_new(1)=tempLine(i);
            else
                %结构体不为空时可以用常规方法添加新对象
                tempLine_new=[tempLine_new,tempLine(i)];
            end
        end
        clear tempLine
        tempLine=tempLine_new;
        if exist('newLine2')==0
            newLine2=tempLine;
        else
            newLine2=[newLine2,tempLine];
        end

        % ----------（content012）----------------------------------------------------------------------------------------------

        % ----------（content013）提取所有新线元素对应的上游节点，存到临时数组upperPoints(i)，i = from 1 to n------------
        % 操作应该和上面相似，这个更需要判断和去重
        upperPoint_IdArray={tempLine.fromNode}; % 现在是字符串了
        %根据id取出所有的上游节点
        clear upperPoint_struct_temp
        for i=1:size(upperPoint_IdArray,2)
            %查找up id对应的point在shp_point的数据
            clear temp_point
            temp_point=shp_pointsAndSubca(find(strcmp({shp_pointsAndSubca.id},upperPoint_IdArray(i))));
            if exist('upperPoint_struct_temp')==0
                upperPoint_struct_temp=temp_point;
            else
                upperPoint_struct_temp=[upperPoint_struct_temp,temp_point];
            end
        end
        %逐个判断取出的upperPoint_struct_temp中与newPoint是否存在重复
        for i=1:size(upperPoint_struct_temp,2)
            judge_exist=length(find(strcmp({newPoint2.id},upperPoint_struct_temp(i).id)));
            if judge_exist>0 %表示存在重复
                continue;
            else %judge_exist==0，不重复
                %newPoint有流域出流点了，不是空的
                %error-20230303001：upperPoint_struct_temp(i)没有markPoint的参数和属性，不能直接添加
                upperPoint_i_temp=upperPoint_struct_temp(i);
                upperPoint_i_temp(1).markPoint=0;
                newPoint2=[newPoint2,upperPoint_i_temp];
                % 判断这个点是不是子汇水区，是的话更新可连通的子汇水区的队列
                % if  length(find((strcmp({saveWatershed.id},upperPoint_i_temp.id))))>0
                if  sum(ismember([saveWatershed.id] , str2num(upperPoint_i_temp.id)))>0
                    if exist('newWatershed')==0
                        newWatershed=saveWatershed([saveWatershed.id]==str2num(upperPoint_i_temp.id));
                        % newWatershed=temp_point;
                    else
                        temp_watershed=saveWatershed([saveWatershed.id]==str2num(upperPoint_i_temp.id));
                        newWatershed=[newWatershed,temp_watershed];
                    end
                end

            end
        end
        % ----------（content013）-------------------------------------------------------------------------------------------------------
    end
    % 选取下一个没有遍历的Point，and mark the previous one as '1'
    newPoint2(k).markPoint=1;
    k=k+1; %逐个遍历newPoint里的元素点
    judge_loop=sum([newPoint2.markPoint]==0);
    
end

%%
% ----------（content014）-------------------------------------------------------------------------------------------------
%
%绘制shp图
% https://blog.csdn.net/weixin_44058353/article/details/123241109
% https://blog.csdn.net/slandarer/article/details/123098323

% Lat_min = min([newPoint2.X]);
% Lat_max = max([newPoint2.X]);
% Lon_min = min([newPoint2.Y]);
% Lon_max = max([newPoint2.Y]);
% worldmap([Lat_min Lat_max],[Lon_min Lon_max]); 
newLine_draw=newLine2';
fig1=figure(1);
fig1.Position=[1 1 800 800];
hp1(1)=mapshow(newWatershed,'FaceColor','#f2debd');
alpha(0.3)
hp1(2)=mapshow(newLine_draw,'LineWidth',1.5);

xlabel('X坐标','fontsize',15,'FontWeight','bold')
ylabel('Y坐标','fontsize',15,'FontWeight','bold')
legend(hp1,'子汇水区','管线',...
    'Fontsize',12,'Location','northwest')





%% 导出图像要素

% =========================================
% % 01-02-导出新的shp数据

% % 将原始id中数字的部分转换成str
% % 首先是newPoint
% for i=2:length(newPoint)    % 要将outlet排除
%     if isa(newPoint(i).node_id,'double')
%         newPoint(i).node_id=num2str(newPoint(i).node_id);
%     end
% 
% end
% 
% % 首先是newLine
% for i=1:length(newLine)
%     if isa(newLine(i).us_node_id,'double')
%         newLine(i).us_node_id=num2str(newLine(i).us_node_id);
%     end
%     if isa(newLine(i).ds_node_id,'double')
%         newLine(i).ds_node_id=num2str(newLine(i).ds_node_id);
%     end
% 
% end

%watershed
for i=1:length(newWatershed)
    if isa(newWatershed(i).IDofWatershedOrNodes,'double')
        newWatershed(i).IDofWatershedOrNodes=num2str(newWatershed(i).IDofWatershedOrNodes);
    end
    if isa(newWatershed(i).outletID,'double')
        newWatershed(i).outletID=num2str(newWatershed(i).outletID);
    end
    if isa(newWatershed(i).id,'double')
        newWatershed(i).id=num2str(newWatershed(i).id);
    end

end
%%
%
% 导出新的shp数据
clear newPoint_save
newPoint_save(1)=newPoint(2);
for i=2:length(newPoint)-1
    newPoint_save(i)=newPoint(i+1);
end
shapewrite(newPoint_save,"new_points.shp");
shapewrite(newLine,"new_lines.shp");    %生成shp,dbf, shx三个文件
shapewrite(newWatershed,'new_subcatchments.shp');    %生成shp,dbf, shx三个文件

% 导出新的excel表格
newPoint_table = struct2table(newPoint);
newLine_table = struct2table(newLine);
newSubca_table = struct2table(newWatershed);
writetable(newPoint_table, 'newPoint_table.xls')
writetable(newLine_table, 'newLine_table.xls')
writetable(newSubca_table, 'newWatershed_table.xlsx')





%% [coordinates]
% pipe points
clear coordinates_Pipe_point
for i=1:length(newPoint_save)
    coordinates_Pipe_point(i).id=newPoint_save(i).id;
    coordinates_Pipe_point(i).X=newPoint_save(i).X;
    coordinates_Pipe_point(i).Y=newPoint_save(i).Y;
end
data01=struct2table(coordinates_Pipe_point);
writetable(data01,'coordinates_Pipe_poins.txt','Delimiter',' ');

% % river points
% clear coordinates_River_point
% for i=1:size(shp_river_points,1)
%     coordinates_River_point(i).Id=shp_river_points(i).Id;
%     coordinates_River_point(i).X=shp_river_points(i).X;
%     coordinates_River_point(i).Y=shp_river_points(i).Y;
% end
% 
% data01=struct2table(coordinates_River_point);
% writetable(data01,'coordinates_River_point.txt','Delimiter',' ');

%% [verticts]
% pipe lines
clear verticts_pipe
verticts_pipe(1).Id=0;
verticts_pipe(1).X=0;
verticts_pipe(1).Y=0;
num_verticts_pipe=0;
for i=1:length(newLine)
    num_of_points=size(newLine(i).X,2);
    for j=1:num_of_points
        if ~isnan(newLine(i).X(1,j))
            num_verticts_pipe=num_verticts_pipe+1;
            verticts_pipe(num_verticts_pipe).Id=newLine(i).id;
            verticts_pipe(num_verticts_pipe).X=newLine(i).X(1,j);
            verticts_pipe(num_verticts_pipe).Y=newLine(i).Y(1,j);
        end
    end 
end
data02 = struct2table(verticts_pipe);
writetable(data02,'verticts_pipe.txt','Delimiter',' ');

% % river lines
% clear verticts_river
% verticts_river(1).Id=0;
% verticts_river(1).X=0;
% verticts_river(1).Y=0;
% num_verticts_river=0;
% for i=1:size(shp_river_link,1)
%     num_of_points=size(shp_river_link(i).X,2)-1;
%     for j=1:num_of_points
%         num_verticts_river=num_verticts_river+1;
%         verticts_river(num_verticts_river).Id=shp_river_link(i).Id;
%         verticts_river(num_verticts_river).X=shp_river_link(i).X(1,j);
%         verticts_river(num_verticts_river).Y=shp_river_link(i).Y(1,j);
%     end 
% end
% data02 = struct2table(verticts_river);
% writetable(data02,'verticts_river.txt','Delimiter',' ');

%%  [polygon]
clear polygon_subcatchments
polygon_subcatchments(1).Id=0;
polygon_subcatchments(1).X=0;
polygon_subcatchments(1).Y=0;
num_polygon_subcatchments=0;
for i=1:length(newWatershed)
    num_of_points=size(newWatershed(i).X,2)-1;
    for j=1:num_of_points
        if ~isnan(newWatershed(i).X(1,j))
            num_polygon_subcatchments=num_polygon_subcatchments+1;
            polygon_subcatchments(num_polygon_subcatchments).Id=newWatershed(i).id;
            polygon_subcatchments(num_polygon_subcatchments).X=newWatershed(i).X(1,j);
            polygon_subcatchments(num_polygon_subcatchments).Y=newWatershed(i).Y(1,j);
        end
    end 
end

data03 = struct2table(polygon_subcatchments);
writetable(data03,'polygon_subcatchments.txt','Delimiter',' ');









