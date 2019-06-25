function []=GCP_AcoRLF()
%%使用极小极大蚁群算法结合RLF算法求解GCP问题%%
%参考论文：An aco algorithm for the graph coloring problem,2008
%程序作者：梁壹厅，2019
%核心思想：蚁群算法通过信息素机制，在RLF中选择进入第k种颜色的顶点时，提供全局信息，从而进一步提高RLF的性能
%route：着色方案，一个维度和顶点数相同的行向量
%name:测试样本文件名
%Ktarget:目标最少颜色数目
%Tau:蚁群算法中的信息素矩阵
%N:顶点数
%Gm:邻接矩阵
%Rho:蚁群算法的衰减因子
%Alpha:蚁群算法全局代价（信息素）权重
%Beta:蚁群算法局部代价权重
%%
%1，问题'gc_4_1','gc_20_1','gc_50_7','gc_100_5'
close all
clear all
names={'gc_4_1','gc_20_1','gc_50_7','gc_100_5'};
Ktarget_all=[2,3,14,17];
for i=1:numel(names);
    name=names{i};Ktarget=Ktarget_all(i);
    disp(['图',name,'求解结果:']);
    [Kbest,route]=gcp_acorlf(name,Ktarget);
    %disp(['最优颜色数目=',num2str(Kbest)]);
    write_to_file(route,name);
    disp(' ');
end
%%
%2,问题'gc_1000_5'
%较大规模的问题'gc_1000_5'用上面的蚁群算法求解结果不好，只有115多 
%反而直接使用RLF算法即可获得比要求的110更好的解（RLF本身在选取每种颜色的第一个节点的时候就有随机因素在里面），说明这里全局信息的作用很小
%为了让结果更优，这里直接进行随机搜索。
%目前最优的解是104种颜色
name='gc_1000_5';
Ktarget=104;
Kbest=Ktarget*1000;
maxiter=400;%最大随机搜索次数
for iter=1:maxiter;
    if iter==1;flag=1;else flag=0;end
    [Ki,route_i]=RLF(name,Ktarget,flag);%适用于大规模问题
    if Kbest>=Ki;
        Kbest=Ki;
        route=route_i;
    end
    if mod(iter,10)==0;
        %disp(['迭代轮次=',num2str(iter),',最优颜色数目=',num2str(Kbest)]);
    end
end
disp(['图',name,'求解结果:']);
disp(['迭代轮次=',num2str(iter),',最优颜色数目=',num2str(Kbest)]);
write_to_file(route,name);

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

function [Kbest,route_best]=gcp_acorlf(name,Ktarget)
%使用极小极大蚁群算法结合RLF算法求解GCP问题
[N,Gm,G]=prepro(name);%预处理输
Rho=0.990;Alpha=2;Beta=4;%蚁群参数
Nc=50;%蚂蚁数量
Kbest=N;%当前最优颜色数
ckbest=[];%储存每一轮的Kbest
for iter=1:100
    %信息素上下限
    Taumax=1/(Kbest-Ktarget+1)/Rho;%
    Taumin=Taumax/10;
    if iter==1;
        Tau=Taumax*ones(N,N);
    end
    
    Tau_new=Tau;
    for i=1:Nc
        %RLF操作
        [K1,route,Tau1]=RLF_aco(Tau,Ktarget,N,Gm,Rho,Alpha,Beta);
        if K1<=Kbest;
            %最小最大蚁群算法（MMAS),用最优的那个更新信息素
            Tau_new=Tau1;
            Kbest=K1;
            route_best=route;
        end
    end
    Tau=min(Taumax,max(Taumin,Tau_new));
    
    %作图和显示
    ckbest=[ckbest,Kbest];
    %disp(['迭代轮次=',num2str(iter),',最少颜色数目=',num2str(Kbest)]);
    figure(1);plot(ckbest); xlabel('迭代次数');ylabel('颜色数');
    pause(0.1);
    
    %终止条件
    if Kbest<=Ktarget;
        %如果已经达到目标颜色数目就停止迭代
        break;
    end
end
disp(['迭代轮次=',num2str(iter),',最优颜色数目=',num2str(Kbest)]);

function [N,Gm,G]=prepro(name)%预处理函数，负责读取数据，并处理成想要的数据格式
path='./data sets/';%'C:\Users\liang\Desktop\MATLAB\图着色问题\data sets\';
G=textread([path,name]);
N=G(1,1);M=G(1,2);
G=G(2:end,:)+1;%matlab编号从1开始
Gm=zeros(N,N);%邻接矩阵
for i=1:M;
    Gm(G(i,1),G(i,2))=1;
end
function [Kbest,route,Tau,dcc]=RLF_aco(Tau,Ktarget,N,Gm,Rho,Alpha,Beta)
%图着色问题的RLF寻优算法，比另一种直接寻优算法DSATUR好
Kmax=Ktarget*2; 
unvisited=ones(1,N);
route=zeros(1,N);
CC={};
for k=1:Kmax;
    Ck=[];
    for i=1:N;
        if numel(Ck)>0;
            li=sum(Gm(Ck,:),1)+sum(Gm(:,Ck),2)';
            Wk=find((li==0).*(unvisited==1));%允许集合:1,不能和Ck中的顶点有连接，2，不能是已经被涂过的
        else
            Wk=find(unvisited==1);
        end
        %%
        %基于信息素进行路径选择
        if numel(Wk)<1;
            break;
        end
        if numel(Ck)==0;
            i_chosen=randi(numel(Wk));%随机选点
            v=Wk(i_chosen);
        else
            %ACORLF的核心
            %用信息素辅助选择进入k颜色的下一个节点
            %原理是计算新节点与k颜色中已有元素的相同关系信息素之和，统计
            degWk=sum(Gm(Wk,Wk),1)+sum(Gm(Wk,Wk),2)';
            eta=numel(Wk)-degWk;
            tau=(sum(Tau(Ck,Wk),1)+sum(Tau(Wk,Ck),2)')/numel(Ck)/2;%计算允许集合中元素与Ck中元素的信息累积量
            v=Chosen_aco(Wk,tau,eta,Alpha,Beta);%蚁群算法中的路径选择算法
        end
        %%
        route(v)=k;Ck=[Ck,v];
        unvisited(v)=0;
    end
    CC{k}=Ck;
    if sum(unvisited)==0;
        break;
    end
end
Kbest=numel(unique(route));  
dcc=0;
for i=1:numel(CC);
    dcc=dcc+numel(CC{i})^2;
end
Tau=Pheromone(Tau,CC,Kbest,Ktarget,Rho,dcc);%蚁群算法信息素更新
function Tau=Pheromone(Tau,CC,Kbest,Ktarget,Rho,dcc)
%蚁群算法信息素更新模块
for k=1:numel(CC);
    Ck=CC{k};
    Tau(Ck,Ck)=(1-Rho)*Tau(Ck,Ck)+1/(Kbest-Ktarget+1);
end
function to_visit=Chosen_aco(V_next,tau,eta,Alpha,Beta)
%蚁群算法中的路径选择算法
P=(tau/(mean(tau)).^Alpha).*(eta.^Beta);
P=P/(sum(P));
%按概率原则，轮盘赌方法选取下一个节点
Pcum=cumsum(P);     %cumsum，元素累加即求和
Select=find(Pcum>=rand); %若计算的概率大于原来的就选择这条路线
to_visit=V_next(Select(1));
function [Kbest,route]=RLF(name,Ktarget,flag)
%用RLF算法求解最优的着色方案
%这里把数据预处理也写进来了
%数据预处理
if flag==1
    [N,Gm,~]=prepro(name);
    save('rfl_preprocess.mat');
else
    load('rfl_preprocess.mat');
end

%RLF算法部分
Kmax=Ktarget*2;
unvisited=ones(1,N);%未访问点
route=zeros(1,N);%
CC={};
for k=1:Kmax;
    Ck=[];
    for i=1:N;
        if numel(Ck)>0;
            li=sum(Gm(Ck,:),1)+sum(Gm(:,Ck),2)';
            Wk=find((li==0).*(unvisited==1));%允许集合:1,不能和Ck中的顶点有连接，2，不能是已经被涂过的
        else
            Wk=find(unvisited==1);
        end
        if numel(Wk)<1;
            break;
        end
        if numel(Ck)==0;
            i_chosen=randi(numel(Wk));%随机选点
        else
            %RLF算法核心：选择进入颜色k集合的顶点
            degWk=sum(Gm(Wk,Wk),1)+sum(Gm(Wk,Wk),2)';%局部优化条件
            [~,i_chosen]=max(numel(Wk)-degWk);%其实是找在剩余可行集合中度最小的那个
        end
        v=Wk(i_chosen);
        route(v)=k;Ck=[Ck,v];
        unvisited(v)=0;
    end
    CC{k}=Ck;
    if sum(unvisited)==0;
        break;
    end
end
Kbest=numel(unique(route));
