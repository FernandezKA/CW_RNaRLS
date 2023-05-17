%Очистка
clc;
clear all;
close all;
% Исходные данные
G0=400; % КНД антенны
kn=2.4; % Шумовой коэффициент
Pp=3500; % Вт, средняя мощность одного импульса РЛС
sigma=10; % м^2, ЭПР цели
D=0.95; % Вероятность правильного обнаружения
F=1*10^-8; % Вероятность ложной тревоги
X=210*10^3; % м, Длинна покрываемого квадрата
Y=120*10^3; % м, Ширина покрываемого квадрата
k=1.37*10^-23; %Вт/(Гц*К), постоянная Больцмана
T0=290; %К, шумовая температура антенны
Ta = 24; 
tau_p = 3 * 10^(-6); 
f = 170 * 10^6; 
lambda = 3*10^8/f; %c/f
del_F = (0.8 / tau_p);
Tn = Ta + T0 * (kn - 1);

%Расчет энергии принимаемого сигнала
r_max = sqrt((X/2).^2+(Y/2).^2); %Максимально удаленная точка
Etr = Pp * tau_p; %Энергия излученного сигнала
Eres = (Etr*(G0.^2)*sigma*(lambda.^2))/((4*pi).^3*(r_max.^4));%Энергия принимаемого сигнала
%Расчет вероятности правильного обнаружения
Pn=k*del_F*Tn ;%Суммарная мощность шумов приемника
N0=2*Pn/del_F ;%СПМ шума
q=sqrt(2*Eres/N0); %ОСШ
Drasch=F.^(1/(1+0.5*q.^2)); %Вероятность обнаружения


%------------------------—
%Построение кривой обнаружения
%Задаем координатную сетку
Lx=-X/2:100:X/2; %по оси X
Ly=-Y/2:100:Y/2; %по оси Y

%Тестовый расчет по одной оси (x)
%rr=0:1:100000;
%qq=(G0*lambda./(4*pi.^1.5*rr.^2))*sqrt((Etr*sigma)/(2*N0));
%DraschTest=F.^(1./(1+0.5*qq.^2));
%plot(rr,DraschTest);

%Расчет вероятностей обьнаружения по координатной сетке
for i = 1:length(Ly)
for j = 1:length(Lx)
r_odin(i,j)=sqrt(Lx(j)^2+Ly(i)^2); %Поиск расстояния
q_odin(i,j)=(G0*lambda/(4*pi^1.5*r_odin(i,j)^2))*sqrt((Etr*sigma)/(2*N0)); %ОСШ в заданной точке
Dr_odin(i,j)=F^(1/(1+0.5*q_odin(i,j)^2)); %Вероятность ПО
end
end
%Вывод получившейся кривой
figure('Name','Кривая обнаружения одной РЛС');
contour(Lx,Ly,Dr_odin,'ShowText','on');
title('Кривая обнаружения одной РЛС');
xlabel('X,м'); ylabel('Y,м');

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

disp("Значение оценки процента площади с обнаружением более 95%")
disp(counter * 100 / max_val)
%------------------------—

%Расчет ВПО для группы РЛС
nX=2; % Количество РЛС по оси X
nY=3; % Количество РЛС по оси Y
M_rls=nX*nY; % Общее кол-во РЛС
krit = 2; %Критерий для схемы Браннера
Rast=6000; % Расстояние между РЛС
%Создаем координатную сетку с расположением всех РЛС на сетке
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

%Проверка условий четности кол-ва РЛС (необходимо для совпадения матриц по размеру)
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

% Создаем булевую матрицу в которой еденицы есть РЛС
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

%Создаем еденичные матрицы для "и" и "или"
OR=ones(length(LLy),length(LLx));
AND=ones(length(LLy),length(LLx));
ITOG=zeros(length(LLy),length(LLx));
ITOGF=zeros(length(LLy),length(LLx));
ITOGF1=zeros(length(LLy),length(LLx));
%Проводим расчет ВПО для каждой РЛС
figure();
hold on;
schet=0;
for k = 1:length(LLLx)
for h = 1:length(LLLy)
if boolmat (k,h) == 1
%Создаем новую координатную сетку для КАЖДОЙ РЛС
RastX = (-100*(k-1)):100:(100*(length(LLLx)-k));
RastY = (-100*(h-1)):100:(100*(length(LLLy)-h));
%Проводим необходимые расчеты
for i = 1:length(RastY)
for j = 1:length(RastX)
RRR(i,j)=sqrt(RastX(j)^2+RastY(i)^2);
QQ1(i,j)=(G0*lambda/(4*pi^1.5*RRR(i,j)^2))*sqrt((Etr*sigma)/(2*N0));
DraschTest11(i,j)=F^(1/(1+0.5*QQ1(i,j)^2));
end
end
%Выводим области действия всех РЛС
contour(LLx,LLy,DraschTest11,'ShowText','on');
%РАсчет вероятностей по схеме И и ИЛИ для всех РЛС
OR=OR.*(1-DraschTest11);
AND=AND.*DraschTest11;
%Матрицы для Браннера
 
schet=schet+1;
Dbr(:,:,schet)=DraschTest11;
Qbr(:,:,schet)=1-DraschTest11;
end
end
end
title('Кривые обнаружения для всех РЛС');
xlabel('X,м'); ylabel('Y,м');

%Схема браннера
%Определение матриц вероятностей
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
%ВЛТ
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
    disp('оптимальный критерий')
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

otobr1 = ['оптимальный критерий :',num2str(save_ind), '/',num2str(schet),...
    '    его ВЛТ составляет', num2str(optimal)];
disp(otobr1)

for i = 1:length(Pf_sum)
    
    otobr2=['остальные критерии: ', num2str(i), '/',num2str(schet), ...
        '    их ВЛТ составляет', num2str(Pf_sum(i))];
    disp(otobr2)
    
end


for i=schet:-1:1
ITOG=ITOG+PP(:,:,i,schet);
figure('Name',['AND_Branner',num2str(i)]);
contour(LLx,LLy,ITOG,'ShowText','on');
title(['Критерий  ',num2str(i),'/',num2str(schet)]);
xlabel('X,м'); ylabel('Y,м');
end
ITOG=zeros(length(LLy),length(LLx));
for i=schet:-1:save_ind
    ITOG=ITOG+PP(:,:,i,schet);
end
figure('Name',['AND_Branner',num2str(save_ind),'optimal']);
contour(LLx,LLy,ITOG,'ShowText','on');
title(['Критерий оптимальный ',num2str(save_ind),'/',num2str(schet)]);
xlabel('X,м'); ylabel('Y,м');


%Вывод графиков
figure('Name','AND');
contour(LLx,LLy,AND,'ShowText','on');
xlabel('X,м'); ylabel('Y,м');
title('Расчет всех РЛС по схеме "И"')
figure('Name','OR');
contour(LLx,LLy,1-OR,'ShowText','on');
xlabel('X,м'); ylabel('Y,м');
title('Расчет всех РЛС по схеме "ИЛИ"')