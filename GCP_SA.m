function []=GCP_SA()
%用模拟退火求解：关键在于随机扰动的方式，以及扰动位点个数的选取
%程序作者：梁壹厅，2019
%route：着色方案，一个维度和顶点数相同的行向量
%name:测试样本文件名
%Ktarget:目标最少颜色数目
%N:顶点数
%Gm:邻接矩阵
%%
close all
clear all
names={'gc_50_5'};
Ktarget_all=[9];
conflicts=50;
for i=1:numel(names);
    name=names{i};Ktarget=Ktarget_all(i);
    while(conflicts>0);
    [route,conflicts]=gcp_sa(name,Ktarget);
    end
    disp(['颜色数目=',num2str(Ktarget),' 冲突次数=',num2str(conflicts)]);
    write_to_file(route,name);
end
%%
%子函数
function []=write_to_file(route,name)
string='';
for i=1:numel(route);
    string=[string,' ',num2str(route(i))];
end
savepath=['./result_of_',name,'.txt'];
fid=fopen(savepath,'w');
fprintf(fid,'%s',string);
disp(['结果存在：',savepath])
fclose(fid);
function [route_best,conflicts]=gcp_sa(name,Ktarget)
%求解图着色问题
%T;温度
%r:温度收缩因子
%route_best：最优的方案

%预处理
[N,~,G]=prepro(name);
calculateCost=@(X)fitness_tzs(X,G);

%进入模拟退火阶段
T0=1e6;r = 0.99;  Ts = 0.0001;  
route=randi(Ktarget,1,N);% 初始化
cost = calculateCost(route);
T = T0;
cost_min = cost;
route_best=route;
iter = 1;
cost_list=[];%储存代表序列
T_list=T;%储存温度序列
figure(1);
while(T > Ts)
    %重要改进1！！！！增强后期局部搜索能力
    if T<0.1; r=0.9999;end
    
    %随机扰动和适应度计算
    mode=3;%选择第三种随机扰动模式
    newRoute = createNeibor(route,mode,Ktarget);%随机扰动产生新的解
    newCost = calculateCost(newRoute);%计算适应度
    
    %Metropolis接受准则
    delta = newCost - cost;
    if(delta < 0)
        cost = newCost;
        route = newRoute;
    else
        p=exp(-delta/T);
        if rand <= p
            cost = newCost;
            route = newRoute;
        end
    end
    
    %退火操作（更新温度）
    T = T*r; %  annealing
    
    %储存信息和作图
    if cost < cost_min
        cost_min = cost;
        route_best=route;
    end
    cost_list=[cost_list,cost];%储存历史适应度值
    T_list=[T_list,T];
    iter = iter+1;
    if mod(iter,50)==1;    
        disp(['#' 'Iteration ' num2str(iter) ': 冲突数目 = ' num2str(cost) ' 温度T = ' num2str(T)]);
    end
    if mod(iter,100)==1;
        %适应度轨迹和温度轨迹图
        figure(1);
        subplot(2,1,1);plot(cost_list);xlabel('迭代次数');ylabel('适应度');
        subplot(2,1,2);plot(T_list);
        pause(0.1);
    end
    
    %终止条件之一
    if cost_min<=0;
        %也就是冲突数目为0的条件下，终止迭代
        break;
    end
end

%输出最优涂色方案的冲突数目
conflicts=cost_min;
function [N,Gm,G]=prepro(name)%预处理函数，负责读取数据，并处理成想要的数据格式
path='./data sets/';
G=textread([path,name]);
N=G(1,1);M=G(1,2);
G=G(2:end,:)+1;%matlab编号从1开始
Gm=zeros(N,N);%邻接矩阵
for i=1:M;
    Gm(G(i,1),G(i,2))=1;
end
function f=fitness_tzs(X,G)
%计算一个涂色方案X对应的冲突数目
[Np,~]=size(X);
max_memory=2e7;
M=size(G,1);
if M*Np<=max_memory;
    ii=G(:,1)';jj=G(:,2)';
    f=sum(X(:,ii)==X(:,jj),2);
else
     %如果样本数目太多，同时边的数目太大导致邻接矩阵超出matlab最大内存限制，则需要分片进行计算
    f=zeros(Np,1);
    m=floor(max_memory/M);
    k0=1;k=m;
    ii=G(:,1)';jj=G(:,2)';
    while k0<Np;
        ind=k0:k;
        f(ind)=sum(X(ind,ii)==X(ind,jj),2);
        k0=k+1;k=min(k+m,Np);
    end
end
function newRoute = createNeibor(route,mode,Ktarget)
switch mode
    case 1
        % 随机交换两个点的值
        newRoute = Swap(route);
    case 2
        % 随机逆转一段序列
        newRoute=Reversion(route);
    case 3 
        % 随机改变一个点的值
        newRoute=Change(route,Ktarget);
end
function newRoute = Change(route,Ktarget)
n=numel(route);
%重要改进2！！！随机扰动方式的选择
k=1;%
i1=randsample(n,k);%随机选取两个点改变其值 
newRoute=route;
newRoute(i1)=randi(Ktarget,1,numel(i1));%随机选取
function newRoute = Swap(route)
n = numel(route);
i = randsample(n,2);
i1 = i(1);
i2 = i(2);
newRoute = route;
newRoute([i1 i2]) = route([i2 i1]);
function newRoute = Reversion(route)
n=numel(route);
i=randsample(n,2);
i1=min(i(1),i(2));
i2=max(i(1),i(2));
newRoute=route;
newRoute(i1:i2)=route(i2:-1:i1);




