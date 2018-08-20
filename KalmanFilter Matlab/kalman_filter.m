% v�ltoz�k deklar�l�sa
dt=0.1;
A=[1 dt 0.5*dt^2; % �llapott�r A m�trixa, csak ez van, nincs szab�lyoz� tag
    0 1 dt;
    0 0 1];
X0=[0; 20; 5]; %becs�lt  �llapotv�ltoz�k (kezd�snek ismert) : [hely, sebess�g,gyorsul�s]
Q=[0.01; 0.01; 0.01]; %folyamat hiba m�trix
P=[Q(1)^2     Q(1)*Q(2)   Q(1)*Q(3); 
    Q(2)*Q(1)    Q(2)^2     Q(2)*Q(3); 
    Q(3)*Q(1)  Q(3)*Q(2)      Q(3)^2]; % kovariancia m�trix
R=[0.5 0 0; 
    0 0.5 0; 
    0 0 0.5]; % m�r�si hiba m�trix
noisesigma=5; %zaj szigm�ja
Measured=[0 ;0 ;0]; % m�rt �llapotok 
X=A*X0;
K=zeros(3);
log=[0 0];

for i=1:length(Measurement)
    Pos_previous=Measured(1);
    Measured(1)=Measurement(i,2)+normrnd(0,noisesigma); % zaj r�vitele a m�r�si jelre

    log_pos(i)=Measured(1);
    Vel_previous=Measured(2);
    Measured(2)=(Measured(1)-Pos_previous)/dt; %m�rt sebess�g kisz�m�t�sa
    Measured(3)=(Measured(2)-Vel_previous)/dt; % m�rt gyorsul�s kisz�m�t�sa 

    P=A*P*A'+Q;
    
    K=P./(P+R);
    
    K(1,2)=0;
    K(1,3)=0;
    K(2,1)=0;
    K(2,3)=0;
    K(3,1)=0;
    K(3,2)=0;
    
    
    
    X=X+K*(Measured-X);

    
    P=(eye(3)-K)*P;

    
    X=A*X;
    log(i,1)=X(1);
    log(i,2)=X(1);
end
for i=1:length(Measurement)
    y1(i)=Measurement(i,2);
end
for i=1:length(log)
    y2(i)=log(i,2);
end

plot(log_pos);
hold on
plot(y2);
legend('M�rt zajos','Sz�rt illesztett')