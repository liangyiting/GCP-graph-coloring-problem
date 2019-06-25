function []=GCP_SA()
%��ģ���˻���⣺�ؼ���������Ŷ��ķ�ʽ���Լ��Ŷ�λ�������ѡȡ
%�������ߣ���Ҽ����2019
%route����ɫ������һ��ά�ȺͶ�������ͬ��������
%name:���������ļ���
%Ktarget:Ŀ��������ɫ��Ŀ
%N:������
%Gm:�ڽӾ���
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
    disp(['��ɫ��Ŀ=',num2str(Ktarget),' ��ͻ����=',num2str(conflicts)]);
    write_to_file(route,name);
end
%%
%�Ӻ���
function []=write_to_file(route,name)
string='';
for i=1:numel(route);
    string=[string,' ',num2str(route(i))];
end
savepath=['./result_of_',name,'.txt'];
fid=fopen(savepath,'w');
fprintf(fid,'%s',string);
disp(['������ڣ�',savepath])
fclose(fid);
function [route_best,conflicts]=gcp_sa(name,Ktarget)
%���ͼ��ɫ����
%T;�¶�
%r:�¶���������
%route_best�����ŵķ���

%Ԥ����
[N,~,G]=prepro(name);
calculateCost=@(X)fitness_tzs(X,G);

%����ģ���˻�׶�
T0=1e6;r = 0.99;  Ts = 0.0001;  
route=randi(Ktarget,1,N);% ��ʼ��
cost = calculateCost(route);
T = T0;
cost_min = cost;
route_best=route;
iter = 1;
cost_list=[];%�����������
T_list=T;%�����¶�����
figure(1);
while(T > Ts)
    %��Ҫ�Ľ�1����������ǿ���ھֲ���������
    if T<0.1; r=0.9999;end
    
    %����Ŷ�����Ӧ�ȼ���
    mode=3;%ѡ�����������Ŷ�ģʽ
    newRoute = createNeibor(route,mode,Ktarget);%����Ŷ������µĽ�
    newCost = calculateCost(newRoute);%������Ӧ��
    
    %Metropolis����׼��
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
    
    %�˻�����������¶ȣ�
    T = T*r; %  annealing
    
    %������Ϣ����ͼ
    if cost < cost_min
        cost_min = cost;
        route_best=route;
    end
    cost_list=[cost_list,cost];%������ʷ��Ӧ��ֵ
    T_list=[T_list,T];
    iter = iter+1;
    if mod(iter,50)==1;    
        disp(['#' 'Iteration ' num2str(iter) ': ��ͻ��Ŀ = ' num2str(cost) ' �¶�T = ' num2str(T)]);
    end
    if mod(iter,100)==1;
        %��Ӧ�ȹ켣���¶ȹ켣ͼ
        figure(1);
        subplot(2,1,1);plot(cost_list);xlabel('��������');ylabel('��Ӧ��');
        subplot(2,1,2);plot(T_list);
        pause(0.1);
    end
    
    %��ֹ����֮һ
    if cost_min<=0;
        %Ҳ���ǳ�ͻ��ĿΪ0�������£���ֹ����
        break;
    end
end

%�������Ϳɫ�����ĳ�ͻ��Ŀ
conflicts=cost_min;
function [N,Gm,G]=prepro(name)%Ԥ�������������ȡ���ݣ����������Ҫ�����ݸ�ʽ
path='./data sets/';
G=textread([path,name]);
N=G(1,1);M=G(1,2);
G=G(2:end,:)+1;%matlab��Ŵ�1��ʼ
Gm=zeros(N,N);%�ڽӾ���
for i=1:M;
    Gm(G(i,1),G(i,2))=1;
end
function f=fitness_tzs(X,G)
%����һ��Ϳɫ����X��Ӧ�ĳ�ͻ��Ŀ
[Np,~]=size(X);
max_memory=2e7;
M=size(G,1);
if M*Np<=max_memory;
    ii=G(:,1)';jj=G(:,2)';
    f=sum(X(:,ii)==X(:,jj),2);
else
     %���������Ŀ̫�࣬ͬʱ�ߵ���Ŀ̫�����ڽӾ��󳬳�matlab����ڴ����ƣ�����Ҫ��Ƭ���м���
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
        % ��������������ֵ
        newRoute = Swap(route);
    case 2
        % �����תһ������
        newRoute=Reversion(route);
    case 3 
        % ����ı�һ�����ֵ
        newRoute=Change(route,Ktarget);
end
function newRoute = Change(route,Ktarget)
n=numel(route);
%��Ҫ�Ľ�2����������Ŷ���ʽ��ѡ��
k=1;%
i1=randsample(n,k);%���ѡȡ������ı���ֵ 
newRoute=route;
newRoute(i1)=randi(Ktarget,1,numel(i1));%���ѡȡ
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




