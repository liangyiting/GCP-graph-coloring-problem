function []=GCP_AcoRLF()
%%ʹ�ü�С������Ⱥ�㷨���RLF�㷨���GCP����%%
%�ο����ģ�An aco algorithm for the graph coloring problem,2008
%�������ߣ���Ҽ����2019
%����˼�룺��Ⱥ�㷨ͨ����Ϣ�ػ��ƣ���RLF��ѡ������k����ɫ�Ķ���ʱ���ṩȫ����Ϣ���Ӷ���һ�����RLF������
%route����ɫ������һ��ά�ȺͶ�������ͬ��������
%name:���������ļ���
%Ktarget:Ŀ��������ɫ��Ŀ
%Tau:��Ⱥ�㷨�е���Ϣ�ؾ���
%N:������
%Gm:�ڽӾ���
%Rho:��Ⱥ�㷨��˥������
%Alpha:��Ⱥ�㷨ȫ�ִ��ۣ���Ϣ�أ�Ȩ��
%Beta:��Ⱥ�㷨�ֲ�����Ȩ��
%%
%1������'gc_4_1','gc_20_1','gc_50_7','gc_100_5'
close all
clear all
names={'gc_4_1','gc_20_1','gc_50_7','gc_100_5'};
Ktarget_all=[2,3,14,17];
for i=1:numel(names);
    name=names{i};Ktarget=Ktarget_all(i);
    disp(['ͼ',name,'�����:']);
    [Kbest,route]=gcp_acorlf(name,Ktarget);
    %disp(['������ɫ��Ŀ=',num2str(Kbest)]);
    write_to_file(route,name);
    disp(' ');
end
%%
%2,����'gc_1000_5'
%�ϴ��ģ������'gc_1000_5'���������Ⱥ�㷨��������ã�ֻ��115�� 
%����ֱ��ʹ��RLF�㷨���ɻ�ñ�Ҫ���110���õĽ⣨RLF������ѡȡÿ����ɫ�ĵ�һ���ڵ��ʱ�����������������棩��˵������ȫ����Ϣ�����ú�С
%Ϊ���ý�����ţ�����ֱ�ӽ������������
%Ŀǰ���ŵĽ���104����ɫ
name='gc_1000_5';
Ktarget=104;
Kbest=Ktarget*1000;
maxiter=400;%��������������
for iter=1:maxiter;
    if iter==1;flag=1;else flag=0;end
    [Ki,route_i]=RLF(name,Ktarget,flag);%�����ڴ��ģ����
    if Kbest>=Ki;
        Kbest=Ki;
        route=route_i;
    end
    if mod(iter,10)==0;
        %disp(['�����ִ�=',num2str(iter),',������ɫ��Ŀ=',num2str(Kbest)]);
    end
end
disp(['ͼ',name,'�����:']);
disp(['�����ִ�=',num2str(iter),',������ɫ��Ŀ=',num2str(Kbest)]);
write_to_file(route,name);

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

function [Kbest,route_best]=gcp_acorlf(name,Ktarget)
%ʹ�ü�С������Ⱥ�㷨���RLF�㷨���GCP����
[N,Gm,G]=prepro(name);%Ԥ������
Rho=0.990;Alpha=2;Beta=4;%��Ⱥ����
Nc=50;%��������
Kbest=N;%��ǰ������ɫ��
ckbest=[];%����ÿһ�ֵ�Kbest
for iter=1:100
    %��Ϣ��������
    Taumax=1/(Kbest-Ktarget+1)/Rho;%
    Taumin=Taumax/10;
    if iter==1;
        Tau=Taumax*ones(N,N);
    end
    
    Tau_new=Tau;
    for i=1:Nc
        %RLF����
        [K1,route,Tau1]=RLF_aco(Tau,Ktarget,N,Gm,Rho,Alpha,Beta);
        if K1<=Kbest;
            %��С�����Ⱥ�㷨��MMAS),�����ŵ��Ǹ�������Ϣ��
            Tau_new=Tau1;
            Kbest=K1;
            route_best=route;
        end
    end
    Tau=min(Taumax,max(Taumin,Tau_new));
    
    %��ͼ����ʾ
    ckbest=[ckbest,Kbest];
    %disp(['�����ִ�=',num2str(iter),',������ɫ��Ŀ=',num2str(Kbest)]);
    figure(1);plot(ckbest); xlabel('��������');ylabel('��ɫ��');
    pause(0.1);
    
    %��ֹ����
    if Kbest<=Ktarget;
        %����Ѿ��ﵽĿ����ɫ��Ŀ��ֹͣ����
        break;
    end
end
disp(['�����ִ�=',num2str(iter),',������ɫ��Ŀ=',num2str(Kbest)]);

function [N,Gm,G]=prepro(name)%Ԥ�������������ȡ���ݣ����������Ҫ�����ݸ�ʽ
path='./data sets/';%'C:\Users\liang\Desktop\MATLAB\ͼ��ɫ����\data sets\';
G=textread([path,name]);
N=G(1,1);M=G(1,2);
G=G(2:end,:)+1;%matlab��Ŵ�1��ʼ
Gm=zeros(N,N);%�ڽӾ���
for i=1:M;
    Gm(G(i,1),G(i,2))=1;
end
function [Kbest,route,Tau,dcc]=RLF_aco(Tau,Ktarget,N,Gm,Rho,Alpha,Beta)
%ͼ��ɫ�����RLFѰ���㷨������һ��ֱ��Ѱ���㷨DSATUR��
Kmax=Ktarget*2; 
unvisited=ones(1,N);
route=zeros(1,N);
CC={};
for k=1:Kmax;
    Ck=[];
    for i=1:N;
        if numel(Ck)>0;
            li=sum(Gm(Ck,:),1)+sum(Gm(:,Ck),2)';
            Wk=find((li==0).*(unvisited==1));%������:1,���ܺ�Ck�еĶ��������ӣ�2���������Ѿ���Ϳ����
        else
            Wk=find(unvisited==1);
        end
        %%
        %������Ϣ�ؽ���·��ѡ��
        if numel(Wk)<1;
            break;
        end
        if numel(Ck)==0;
            i_chosen=randi(numel(Wk));%���ѡ��
            v=Wk(i_chosen);
        else
            %ACORLF�ĺ���
            %����Ϣ�ظ���ѡ�����k��ɫ����һ���ڵ�
            %ԭ���Ǽ����½ڵ���k��ɫ������Ԫ�ص���ͬ��ϵ��Ϣ��֮�ͣ�ͳ��
            degWk=sum(Gm(Wk,Wk),1)+sum(Gm(Wk,Wk),2)';
            eta=numel(Wk)-degWk;
            tau=(sum(Tau(Ck,Wk),1)+sum(Tau(Wk,Ck),2)')/numel(Ck)/2;%������������Ԫ����Ck��Ԫ�ص���Ϣ�ۻ���
            v=Chosen_aco(Wk,tau,eta,Alpha,Beta);%��Ⱥ�㷨�е�·��ѡ���㷨
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
Tau=Pheromone(Tau,CC,Kbest,Ktarget,Rho,dcc);%��Ⱥ�㷨��Ϣ�ظ���
function Tau=Pheromone(Tau,CC,Kbest,Ktarget,Rho,dcc)
%��Ⱥ�㷨��Ϣ�ظ���ģ��
for k=1:numel(CC);
    Ck=CC{k};
    Tau(Ck,Ck)=(1-Rho)*Tau(Ck,Ck)+1/(Kbest-Ktarget+1);
end
function to_visit=Chosen_aco(V_next,tau,eta,Alpha,Beta)
%��Ⱥ�㷨�е�·��ѡ���㷨
P=(tau/(mean(tau)).^Alpha).*(eta.^Beta);
P=P/(sum(P));
%������ԭ�����̶ķ���ѡȡ��һ���ڵ�
Pcum=cumsum(P);     %cumsum��Ԫ���ۼӼ����
Select=find(Pcum>=rand); %������ĸ��ʴ���ԭ���ľ�ѡ������·��
to_visit=V_next(Select(1));
function [Kbest,route]=RLF(name,Ktarget,flag)
%��RLF�㷨������ŵ���ɫ����
%���������Ԥ����Ҳд������
%����Ԥ����
if flag==1
    [N,Gm,~]=prepro(name);
    save('rfl_preprocess.mat');
else
    load('rfl_preprocess.mat');
end

%RLF�㷨����
Kmax=Ktarget*2;
unvisited=ones(1,N);%δ���ʵ�
route=zeros(1,N);%
CC={};
for k=1:Kmax;
    Ck=[];
    for i=1:N;
        if numel(Ck)>0;
            li=sum(Gm(Ck,:),1)+sum(Gm(:,Ck),2)';
            Wk=find((li==0).*(unvisited==1));%������:1,���ܺ�Ck�еĶ��������ӣ�2���������Ѿ���Ϳ����
        else
            Wk=find(unvisited==1);
        end
        if numel(Wk)<1;
            break;
        end
        if numel(Ck)==0;
            i_chosen=randi(numel(Wk));%���ѡ��
        else
            %RLF�㷨���ģ�ѡ�������ɫk���ϵĶ���
            degWk=sum(Gm(Wk,Wk),1)+sum(Gm(Wk,Wk),2)';%�ֲ��Ż�����
            [~,i_chosen]=max(numel(Wk)-degWk);%��ʵ������ʣ����м����ж���С���Ǹ�
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
