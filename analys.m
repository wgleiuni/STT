step=2500;
nx=zeros(step,20000);
figure
set(gcf,'color','w')
hold on
for i=1:step
    eval(['load data/data27/nx',num2str(i),'.mat']);
    temp=eval(['nx',num2str(i)]);
    if ~isempty(temp)
        nx(i,:)=temp;
        y=findpeaks(temp(end-10000:end));
        x=ones(length(y),1)*(-1+2/step*i);
        scatter(x,y,2,'b')
        eval(['clear nx',num2str(i)])
    end
end
%%
ny=zeros(step,20000);

for i=1:step
    eval(['load data/data27/ny',num2str(i),'.mat']);
    temp=eval(['ny',num2str(i)]);
    if ~isempty(temp)
        ny(i,:)=temp;
        eval(['clear ny',num2str(i)])
    end
end
%%
nz=zeros(step,20000);

for i=1:step
    eval(['load data/data27/nz',num2str(i),'.mat']);
    temp=eval(['nz',num2str(i)]);
    if ~isempty(temp)
        nz(i,:)=temp;
        eval(['clear nz',num2str(i)])
    end
end

%%
sx=zeros(step,1000000);

for i=1:step
    eval(['load data/data24/s1',num2str(i),'.txt']);
    temp=eval(['s1',num2str(i)]);
    if ~isempty(temp)
        sx(i,:)=temp(:,2);
        eval(['clear s1',num2str(i)])
    end
end
%%
rule=1;
omega=10;
tau=2*pi/(omega);
dt=0.01;
numtau=floor(tau/dt);
cut=1000;
%%
figure
set(gcf,'color','w')
hold on
if rule==1
    for i=1:step
        y=findpeaks(nx(i,end-cut:end));
        x=ones(length(y),1)*i;
        scatter(x,y,2,'b')
    end
elseif rule==2
    for i=1:step
        y=nx(i,end-cut:numtau:end);
        x=ones(length(y),1)*i;
        scatter(x,y,2,'b')
    end
end

%%
%figure
%set(gcf,'color','w')
%hold on
if rule==1
    for i=1:step
        y=findpeaks(ny(i,end-cut:end));
        x=ones(length(y),1)*i;
        scatter(x,y,2,'r')
    end
elseif rule==2
    for i=1:step
        y=ny(i,end-cut:numtau:end);
        x=ones(length(y),1)*i;
        scatter(x,y,2,'r')
    end
end

%%
% figure
% set(gcf,'color','w')
% hold on
if rule==1
    for i=1:step
        y=findpeaks(nz(i,end-cut:end));
        x=ones(length(y),1)*i;
        scatter(x,y,2,'m')
    end
elseif rule==2
    for i=1:step
        y=nz(i,end-cut:numtau:end);
        x=ones(length(y),1)*i;
        scatter(x,y,2,'m')
    end
end

%%
figure
set(gcf,'color','w')
hold on
if rule==1
    for i=1:step
        y=findpeaks(sx(i,end-cut:end));
        x=ones(length(y),1)*i;
        scatter(x,y,2,'b')
    end
elseif rule==2
    for i=1:step
        y=nz(i,end-cut:numtau:end);
        x=ones(length(y),1)*i;
        scatter(x,y,2,'m')
    end
end

%% calculate sigma for some time
path(path,'~/Documents/work/STT/')
j=847;
dc=-1+2/step*j;
ac=1;
len=length(nx(1,:));
if 1==2
    sigma1=zeros(step,cut);
    sigma2=zeros(step,cut);
    for j=1:step
        for i=1:cut
            sigma1(j,i)=integral(@(x)STT_sigma1(x,nx(j,end-cut+i),ny(j,end-cut+i),nz(j,end-cut+i)),-pi/2,pi/2);
            sigma2(j,i)=integral(@(x)STT_sigma2(x,nx(j,end-cut+i),ny(j,end-cut+i),nz(j,end-cut+i)),-pi/2,pi/2);
        end
    end
elseif 1==1
    sigma1=zeros(1,cut);
    sigma2=zeros(1,cut);
    inte1=zeros(1,cut);
    inte2=zeros(1,cut);
    for i=1:cut
        inte1(i)=integral(@(x)STT_sigma1(x,nx(j,end-cut+i),ny(j,end-cut+i),nz(j,end-cut+i)),-pi/2,pi/2);
        inte2(i)=integral(@(x)STT_sigma2(x,nx(j,end-cut+i),ny(j,end-cut+i),nz(j,end-cut+i)),-pi/2,pi/2);
        sigma1(i)=(dc+ac*cos(omega*(len-cut+i)*dt))*inte1(i);
        sigma2(i)=(dc+ac*cos(omega*(len-cut+i)*dt))*inte2(i);
    end
end

%%
j=851;
i=10540;
theta=-pi/2:pi/1000:pi/2;
y1=STT_sigma1(theta,nx(j,end-i+1),ny(j,end-i+1),nz(j,end-i+1))./sin(theta);
y2=STT_sigma2(theta,nx(j,end-i+1),ny(j,end-i+1),nz(j,end-i+1))./cos(theta);
figure
plot(theta,y1,'b',theta,y2,'r')

%%
theta=0:2*pi/1000:2*pi;
nz1=0.9;
r=sqrt(1-nz1^2);
nx1=r*cos(theta);
ny1=r*sin(theta);
yy=zeros(length(theta),1);
for i=1:length(theta)
    yy(i)=integral(@(x)STT_sigma1(x,nx1(i),ny1(i),nz1),-pi/2,pi/2);
end

%% 3D for phase space
x1=nx(:,1);
y1=ny(:,1);
z1=nz(:,1);
t=mean(ny(:,end-cut:end),2);
X1=reshape(x1,50,50);
Y1=reshape(y1,50,50);
Z1=reshape(z1,50,50);
C=reshape(t,50,50);
surf(X1,Y1,Z1,C)
set(gcf,'color','w')
