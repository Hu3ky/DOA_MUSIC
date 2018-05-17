clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%常量%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
derad = pi/180;       
radeg = 180/pi;
twpi = 2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%信号参数设定%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fc=2.2*10^9;%载波频率
Fs=5*10^9;%采样速率
Rb=1*10^5;%码元速率
C=3*10^8;%光速
L=C/Fc;%波长
SNR=5;%信噪比
Len=500;%快拍数
Theta=45;%入射仰角角度
Psi=45;%入射方位角角度
NUM=1;%来波信号个数
K=twpi/L;%波数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%角度扫描%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Resolution=0.1;%%分辨率
Theta_start=0;
Theta_end=90;
Theta_test=Theta_start:Resolution:Theta_end;%仰角
Psi_start=0;
Psi_end=360;
Psi_test=Psi_start:Resolution:Psi_end;%方位角
Theta_len=length(Theta_test);
Psi_len=length(Psi_test);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%阵列大小%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3;%N行
M=3;%M列
Element_Num=N*M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%子阵间隔%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_x=0.5*L;
d_y=0.5*L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%阵列流型矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(Element_Num,1);
for count_1=1:M    
    for count_2=1:N    
        A(count_1+(count_2-1)*M)=exp(-1i*K*((count_1-1)*d_x*sin(Theta*derad)*cos(Psi*derad)+(count_2-1)*d_y*sin(Theta*derad)*sin(Psi*derad)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%所有的a%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_all=cell(Theta_len,Psi_len);
for count_1=1:Theta_len
    for count_2=1:Psi_len
        a=zeros(Element_Num,1);
        for count_3=1:M    
            for count_4=1:N    
                a(count_3+(count_4-1)*M)=exp(-1i*K*((count_3-1)*d_x*sin(Theta_test(count_1)*derad)*cos(Psi_test(count_2)*derad)+(count_4-1)*d_y*sin(Theta_test(count_1)*derad)*sin(Psi_test(count_2)*derad)));
            end
        end
        a_all{count_1,count_2}=a;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%产生发射信号（BPSK）%%%%%%%%%%%%%%%%%%%%%%%%
S=PskMod(2,Fs,Rb,Fc,Len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%接收信号%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=zeros(Element_Num,Len);
for count_3=1:(Element_Num)
    X(count_3,:)=awgn(A(count_3).*S,SNR,'measured');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%实际的相关矩阵R%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=X*X'/Len;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%对R进行特征值分解，并构造G矩阵%%%%%%%%%%%%%%
[V,D]=eig(R);%%其中D为特征值构成的对角阵，每个特征值对应V矩阵中列向量；
EVA=diag(D)';%%取出所有的特征值
[EVA,I]=sort(EVA);%%进行从小到大排序，并记录索引数组
EVA=fliplr(EVA);%%进行翻转（特征值从大到小排序）
V=fliplr(V(:,I));%%排序完特征值对应的特征向量
G=V(:,(NUM+1):Element_Num);%%构成G矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%求得谱函数矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(Theta_len,Psi_len);
for count_1=1:Theta_len
    for count_2=1:Psi_len
        a=a_all{count_1,count_2};
        P(count_1,count_2)=1./abs(a'*G*G'*a);
    end
end
[signal_1,signal_2]=find(P==max(max(P)));
Theta_find=Theta_start+(signal_1-1)*Resolution;
Psi_find=Psi_start+(signal_2-1)*Resolution;

            



