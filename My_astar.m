%这是自己编纂的代码 不是大华写的
clc,clear;
tic
load("all.mat");%矩阵C仔细看一眼，分清楚前后置集关系,nm是P的前置集，np是P的后置集
C=-C;
M0(2)=4;
%%初始化矩阵信息，定义终态标识Mf（只需要满足仓库库所的零件值为20，故为P(11)==20）
open_list=[M0];%记录待检验的标识的集合
close_list=[];%与open_list相对
Mn=[];
M_all=[M0];
M_temp=[];
parent=zeros(1,size(M0,2));%注意这个parent，起点是没有parent的
fire_t=[0];
%初始化矩阵信息完毕
G_baocun=[];
G_close=[];
H_close=[];
F_close=[];
%%初始化启发式算法
G(1)=0;
Topen=[];%用于保存时间变化信息
H(1)=0;
T1=[8 9 9 5 5 4 8 4 8 12 12 8 4 0 5 8 9];%不变的T，作为标准方便返回
Topen=T1;%变化的T，具体保存为矩阵
F(1)=G(1)+Hsis(M0,100000);
manhat=[];
avg_weight=0;
spend=[];
Mf=M0;
Mf(11)=M0(2);
Mf(2)=0;

%启发式算法初始化完毕

%%主函数部分
while isempty(open_list)==0
    enable_t=[];
    index=find(F==min(F),1);
    G_current=G(index);
    G_baocun=[G_baocun,G_current];
    Mn=open_list(index,:);%从open_list中选取优先级最高的节点n
    T=Topen(index,:);
    if Mn(11)==Mf(11)
        break%如果结点n为终点，则说明算法结束，退出主函数部分
    end
    %%将节点n从open_list中删除，并加入close_list中
        close_list=[close_list;open_list(index,:)];
        open_list(index,:)=[];
         G_close=[G_close,G(index)];
        H_close=[H_close,H(index)];
        F_close=[F_close,F(index)];
       
     enable_t=[];
    %遍历邻近标识
    for i=1:size(np,2)
        if max(Mn,np(:,i)')==Mn%可触发的t
            enable_t=[enable_t,i];
        end
    end
    for i=1:size(np,2)
        if max(Mn,np(:,i)')==Mn%可触发的t
            M_temp=Mn+C(:,i)';%得到邻接结点
             T=Topen(index,:);
            if ismember(M_temp,close_list,'rows')==1%已经在关闭表中，但是未修改父节点
                G_newcost=G_current+T(i);%!!!还米想好每轮G的计算%计算邻接矩阵的时间花费
                manhat=[manhat,norm(Mn - M_temp, 'fro')];
                %时间Petri本意是对的，掺和起来了
               spend=[spend,T(i)];
                avg_weight=mean(spend./manhat);
                if T(i)==0
                    T(i)=T1(i);
                else
                T(enable_t)=T(enable_t)-T(i);
                T(find(T<0))=0;
                end
                [Lia,Locb]=ismember(M_temp,M_all,'rows');
                [Lia,Locb2]=ismember(M_temp,close_list,'rows');
                if G_newcost<G_close(Locb2)
                    parent(Locb,:)=Mn;%如果临近结点gcost更小，则修改父节点
                    fire_t(Locb)=i;
                end
                open_list=[open_list;M_temp];
                G=[G,G_newcost];
                Topen=[Topen;T];
                H=[H,Hsis(M_temp,avg_weight)];
                F=[F,G(end)+H(end)];
                close_list(Locb2,:)=[];

                continue
            elseif ismember(M_temp,open_list,'rows')==1
                G_newcost=G_current+T(i);%!!!还米想好每轮G的计算%计算邻接矩阵的时间花费
                 manhat=[manhat,norm(Mn - M_temp, 'fro')];
                %时间Petri本意是对的，掺和起来了
               spend=[spend,T(i)];
                avg_weight=mean(spend./manhat);
                 if T(i)==0
                    T(i)=T1(i);
                else
                T(enable_t)=T(enable_t)-T(i);
                T(find(T<0))=0;
                end
                [Lia,Locb]=ismember(M_temp,M_all,'rows');
                [Lia,Locb2]=ismember(M_temp,open_list,'rows');
                if G_newcost<G(Locb2)
                    parent(Locb,:)=Mn;%如果临近结点gcost更小，则修改父节点
                    fire_t(Locb)=i;
                    Topen(Locb2,:)=T;
                end
            elseif ismember(M_temp,open_list,'rows')==0
                open_list=[open_list;M_temp];
                G_newcost=G_current+T(i);%!!!还米想好每轮G的计算%计算邻接矩阵的时间花费
                 manhat=[manhat,norm(Mn - M_temp, 'fro')];
                %时间Petri本意是对的，掺和起来了
               spend=[spend,T(i)];
                avg_weight=mean(spend./manhat);
                 if T(i)==0
                    T(i)=T1(i);
                else
                T(enable_t)=T(enable_t)-T(i);
                T(find(T<0))=0;
                 end
                 Topen=[Topen;T];
                G=[G,G_newcost];
                H=[H,Hsis(M_temp,avg_weight)];
                F=[F,G(end)+H(end)];

            end
            if ismember(M_temp,M_all,'rows')==0
                M_all=[M_all;M_temp];
                parent=[parent;Mn];
                fire_t=[fire_t,i];
            end
        end
    end
    Topen(index,:)=[];
        G(index)=[];
        H(index)=[];
        F(index)=[];
end
final_route=[];%保存t的发射序列
final_M=[Mn];%保存所有经过的结点
find_parent=Mn;%用到的临时变量
final_time=G_current;
while isequal(find_parent,zeros(1,size(M0,2)))==0%当搜寻到头结点时，停止
    [Lia,Locb]=ismember(find_parent,M_all,"rows");%寻找在总图中的序号，以找到自己的父节点和前置发射t
    final_route=[final_route,fire_t(Locb)];%记录前置发射t
    final_M=[final_M;parent(Locb,:)];%记录父节点
    find_parent=parent(Locb,:);%更新当前临时变量
end
final_M=fliplr(final_M);
final_route=fliplr(final_route);
toc
    

    