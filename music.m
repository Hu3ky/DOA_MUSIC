clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
derad = pi/180;       
radeg = 180/pi;
twpi = 2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%�źŲ����趨%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fc=2.2*10^9;%�ز�Ƶ��
Fs=5*10^9;%��������
Rb=1*10^5;%��Ԫ����
C=3*10^8;%����
L=C/Fc;%����
SNR=5;%�����
Len=500;%������
Theta=45;%�������ǽǶ�
Psi=45;%���䷽λ�ǽǶ�
NUM=1;%�����źŸ���
K=twpi/L;%����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ƕ�ɨ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Resolution=0.1;%%�ֱ���
Theta_start=0;
Theta_end=90;
Theta_test=Theta_start:Resolution:Theta_end;%����
Psi_start=0;
Psi_end=360;
Psi_test=Psi_start:Resolution:Psi_end;%��λ��
Theta_len=length(Theta_test);
Psi_len=length(Psi_test);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���д�С%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3;%N��
M=3;%M��
Element_Num=N*M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_x=0.5*L;
d_y=0.5*L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�������;���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(Element_Num,1);
for count_1=1:M    
    for count_2=1:N    
        A(count_1+(count_2-1)*M)=exp(-1i*K*((count_1-1)*d_x*sin(Theta*derad)*cos(Psi*derad)+(count_2-1)*d_y*sin(Theta*derad)*sin(Psi*derad)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���е�a%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���������źţ�BPSK��%%%%%%%%%%%%%%%%%%%%%%%%
S=PskMod(2,Fs,Rb,Fc,Len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�����ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=zeros(Element_Num,Len);
for count_3=1:(Element_Num)
    X(count_3,:)=awgn(A(count_3).*S,SNR,'measured');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ʵ�ʵ���ؾ���R%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=X*X'/Len;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��R��������ֵ�ֽ⣬������G����%%%%%%%%%%%%%%
[V,D]=eig(R);%%����DΪ����ֵ���ɵĶԽ���ÿ������ֵ��ӦV��������������
EVA=diag(D)';%%ȡ�����е�����ֵ
[EVA,I]=sort(EVA);%%���д�С�������򣬲���¼��������
EVA=fliplr(EVA);%%���з�ת������ֵ�Ӵ�С����
V=fliplr(V(:,I));%%����������ֵ��Ӧ����������
G=V(:,(NUM+1):Element_Num);%%����G����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����׺�������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

            



