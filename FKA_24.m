%�������
clc;
clear all;
close all;
% �������� ������
G0=400; % ��� �������
kn=2.4; % ������� �����������
Pp=3500; % ��, ������� �������� ������ �������� ���
sigma=10; % �^2, ��� ����
D=0.95; % ����������� ����������� �����������
F=1*10^-8; % ����������� ������ �������
X=210*10^3; % �, ������ ������������ ��������
Y=120*10^3; % �, ������ ������������ ��������
k=1.37*10^-23; %��/(��*�), ���������� ���������
T0=290; %�, ������� ����������� �������
Ta = 24; 
tau_p = 3 * 10^(-6); 
f = 170 * 10^6; 
lambda = 3*10^8/f; %c/f
del_F = (0.8 / tau_p);
Tn = Ta + T0 * (kn - 1);

%������ ������� ������������ �������
r_max = sqrt((X/2).^2+(Y/2).^2); %����������� ��������� �����
Etr = Pp * tau_p; %������� ����������� �������
Eres = (Etr*(G0.^2)*sigma*(lambda.^2))/((4*pi).^3*(r_max.^4));%������� ������������ �������
%������ ����������� ����������� �����������
Pn=k*del_F*Tn ;%��������� �������� ����� ���������
N0=2*Pn/del_F ;%��� ����
q=sqrt(2*Eres/N0); %���
Drasch=F.^(1/(1+0.5*q.^2)); %����������� �����������


%------------------------�
%���������� ������ �����������
%������ ������������ �����
Lx=-X/2:100:X/2; %�� ��� X
Ly=-Y/2:100:Y/2; %�� ��� Y

%�������� ������ �� ����� ��� (x)
%rr=0:1:100000;
%qq=(G0*lambda./(4*pi.^1.5*rr.^2))*sqrt((Etr*sigma)/(2*N0));
%DraschTest=F.^(1./(1+0.5*qq.^2));
%plot(rr,DraschTest);

%������ ������������ ������������ �� ������������ �����
for i = 1:length(Ly)
for j = 1:length(Lx)
r_odin(i,j)=sqrt(Lx(j)^2+Ly(i)^2); %����� ����������
q_odin(i,j)=(G0*lambda/(4*pi^1.5*r_odin(i,j)^2))*sqrt((Etr*sigma)/(2*N0)); %��� � �������� �����
Dr_odin(i,j)=F^(1/(1+0.5*q_odin(i,j)^2)); %����������� ��
end
end
%����� ������������ ������
figure('Name','������ ����������� ����� ���');
contour(Lx,Ly,Dr_odin,'ShowText','on');
title('������ ����������� ����� ���');
xlabel('X,�'); ylabel('Y,�');

%/*****************************/
counter = 0; 
[sz_x, sz_y] = size(Dr_odin); 
max_val = sz_x * sz_y; %All of elements

for i = 1:sz_x 
   for j = 1:sz_y
      if Dr_odin(i, j) >= 0.95
         counter = counter + 1;  
      end 
   end
end

disp("�������� ������ �������� ������� � ������������ ����� 95%")
disp(counter * 100 / max_val)
%------------------------�

%������ ��� ��� ������ ���
nX=2; % ���������� ��� �� ��� X
nY=3; % ���������� ��� �� ��� Y
M_rls=nX*nY; % ����� ���-�� ���
krit = 2; %�������� ��� ����� ��������
Rast=6000; % ���������� ����� ���
%������� ������������ ����� � ������������� ���� ��� �� �����
LL=-Rast/2:100:Rast/2;
for i = 1:length(LL)
if LL(i)==0
LL(i)= 1;
else
LL(i) = 0;
end
end
LLX=repmat(LL,1,nX);
LLY=repmat(LL,1,nY);

LLx=0:100:X;
LLy=0:100:Y;

%�������� ������� �������� ���-�� ��� (���������� ��� ���������� ������ �� �������)
if mod(nX,2) == 1
LLLx=[zeros(1,(length(LLx)-length(LLX))/2) LLX zeros(1,(length(LLx)-length(LLX))/2)];
else
LLLx=[zeros(1,(length(LLx)-length(LLX)-1)/2) LLX zeros(1,(length(LLx)-length(LLX)+1)/2)];
end
if mod(nY,2) == 1
LLLy=[zeros(1,(length(LLy)-length(LLY))/2) LLY zeros(1,(length(LLy)-length(LLY))/2)];
else
LLLy=[zeros(1,(length(LLy)-length(LLY)-1)/2) LLY zeros(1,(length(LLy)-length(LLY)+1)/2)];
end
boolmat=zeros(length(LLLx),length(LLLy));

% ������� ������� ������� � ������� ������� ���� ���
for i = 1:length(LLLx)
for j = 1:length(LLLy)
if LLLx(i) == 1
if LLLy(j) == 1
boolmat(i,j)=1;
end
end
end
end
%figure('Name','COORD');
%(LLLy,LLLx,boolmat);

%������� ��������� ������� ��� "�" � "���"
OR=ones(length(LLy),length(LLx));
AND=ones(length(LLy),length(LLx));
ITOG=zeros(length(LLy),length(LLx));
ITOGF=zeros(length(LLy),length(LLx));
ITOGF1=zeros(length(LLy),length(LLx));
%�������� ������ ��� ��� ������ ���
figure();
hold on;
schet=0;
for k = 1:length(LLLx)
for h = 1:length(LLLy)
if boolmat (k,h) == 1
%������� ����� ������������ ����� ��� ������ ���
RastX = (-100*(k-1)):100:(100*(length(LLLx)-k));
RastY = (-100*(h-1)):100:(100*(length(LLLy)-h));
%�������� ����������� �������
for i = 1:length(RastY)
for j = 1:length(RastX)
RRR(i,j)=sqrt(RastX(j)^2+RastY(i)^2);
QQ1(i,j)=(G0*lambda/(4*pi^1.5*RRR(i,j)^2))*sqrt((Etr*sigma)/(2*N0));
DraschTest11(i,j)=F^(1/(1+0.5*QQ1(i,j)^2));
end
end
%������� ������� �������� ���� ���
contour(LLx,LLy,DraschTest11,'ShowText','on');
%������ ������������ �� ����� � � ��� ��� ���� ���
OR=OR.*(1-DraschTest11);
AND=AND.*DraschTest11;
%������� ��� ��������
 
schet=schet+1;
Dbr(:,:,schet)=DraschTest11;
Qbr(:,:,schet)=1-DraschTest11;
end
end
end
title('������ ����������� ��� ���� ���');
xlabel('X,�'); ylabel('Y,�');

%����� ��������
%����������� ������ ������������
P(:,:,1)=Qbr(:,:,1);
for eprst=2:schet
P(:,:,eprst)=Qbr(:,:,eprst).*P(:,:,eprst-1);
end
PP(:,:,1,1)=Dbr(:,:,1);
for i=2:M_rls
for j=1:i
if i==j
PP(:,:,j,i)=PP(:,:,j-1,i-1).*Dbr(:,:,i);
elseif j==1
PP(:,:,j,i)=Qbr(:,:,i).*PP(:,:,j,i-1)+Dbr(:,:,i).*P(:,:,i-1);
else
PP(:,:,j,i)=Qbr(:,:,i).*PP(:,:,j,i-1)+Dbr(:,:,i).*PP(:,:,j-1,i-1);
end
end
end
%���
%q=1-F
%F1f2f3f4 = F
Pf(1)=1-F;
for eprst=2:schet
    Pf(eprst)=(1-F)*Pf(eprst-1);
end

PPf(1,1)=F;
for i=2:M_rls
for j=1:i
if i==j
    PPf(j,i)=PPf(j-1,i-1).*F;
elseif j==1
            PPf(j,i)=(1-F).*PPf(j,i-1)+F.*Pf(i-1);
else
            PPf(j,i)=(1-F).*PPf(j,i-1)+F.*PPf(j-1,i-1);
end
end
end


for j = 1: size(PPf,1)
    Pf_sum(j) = sum(PPf(j:size(PPf,1),end));
end

if Pf_sum(1)<F
    disp('����������� ��������')
     a = ['0','/' num2str(schet)];
     disp(a)
end
optimal = 0;
for i = 2:length(Pf_sum)
    if (Pf_sum(i)<F)&&(Pf_sum(i-1)>F)
        optimal = Pf_sum(i);
        save_ind = i; 
    end
end

otobr1 = ['����������� �������� :',num2str(save_ind), '/',num2str(schet),...
    '    ��� ��� ����������', num2str(optimal)];
disp(otobr1)

for i = 1:length(Pf_sum)
    
    otobr2=['��������� ��������: ', num2str(i), '/',num2str(schet), ...
        '    �� ��� ����������', num2str(Pf_sum(i))];
    disp(otobr2)
    
end


for i=schet:-1:1
ITOG=ITOG+PP(:,:,i,schet);
figure('Name',['AND_Branner',num2str(i)]);
contour(LLx,LLy,ITOG,'ShowText','on');
title(['��������  ',num2str(i),'/',num2str(schet)]);
xlabel('X,�'); ylabel('Y,�');
end
ITOG=zeros(length(LLy),length(LLx));
for i=schet:-1:save_ind
    ITOG=ITOG+PP(:,:,i,schet);
end
figure('Name',['AND_Branner',num2str(save_ind),'optimal']);
contour(LLx,LLy,ITOG,'ShowText','on');
title(['�������� ����������� ',num2str(save_ind),'/',num2str(schet)]);
xlabel('X,�'); ylabel('Y,�');


%����� ��������
figure('Name','AND');
contour(LLx,LLy,AND,'ShowText','on');
xlabel('X,�'); ylabel('Y,�');
title('������ ���� ��� �� ����� "�"')
figure('Name','OR');
contour(LLx,LLy,1-OR,'ShowText','on');
xlabel('X,�'); ylabel('Y,�');
title('������ ���� ��� �� ����� "���"')