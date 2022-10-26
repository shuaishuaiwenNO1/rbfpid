function [sys,x0,str,ts,simStateCompliance] = PIDRBFs(t,x,u,flag,~)
% 网络结构3-6-1

% 采样时间
Ts=0.001;
switch flag
    case 0
        [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(Ts);
    case 2
        sys=mdlUpdate(u);
    case 3
        sys=mdlOutputs(t,x,u);
    case {1,4,9}
        sys=[];
    otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(Ts)
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 3;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   

sys = simsizes(sizes);
x0  = [0;0;0];
str = [];
ts  = [Ts 0];
simStateCompliance = 'UnknownSimState';

function sys=mdlUpdate(u)

sys=[u(1)-u(2);u(1);u(1)-2*u(2)+u(3)];

function sys=mdlOutputs(t,x,u)
persistent w w_1 w_2 h ci ci_1 ci_2 bi bi_1 bi_2 iq kp0 ki0 kd0 ci_3 bi_3 w_3
xite = 0.5;% RBF网络学习效率以及pid参数学习效率
alfa = 0.05;% RBF网络动量因子
xitekp=0.25;
xitekd=0.25;
xiteki=0.25;%pid学习效率

if t==0
% PID参数初值
kp0=0.01;ki0=0.01;kd0=0.01;
 % 高斯基函数中心矩阵，取值范围在输入信号范围内
 ci=zeros(3,6);
 bi=3*ones(6,1);    %第i个节点的基宽向量B
 w=0.10*ones(6,1);   %权值向量W
 h= zeros(6,1);
ci_1=ci;
ci_3=ci_1;
ci_2=ci_1;    %c(i-1),c(i-3),c(i-2)，迭代用
bi_1=bi;
bi_2=bi_1; 
bi_3=bi_2;     %b(i-1),b(i-2),b(i-3)
w_1=w;
w_2=w_1;  
w_3=w_1;      %w(i-1),w(i-2),w(i-3)
iq=0;                    %u(k-1),y(k-1)，迭代用
                  
end

iq=iq+(kp0*(u(1)-u(2))+ki0*u(1)+kd0*(u(1)-2*u(2)+u(3)));%初始iq

RBF_input = [iq u(4) u(5)]';

for j = 1:1:6
    h(j) = exp(-norm(RBF_input - ci(:,j))^2/(2*bi(j)*bi(j)));
end
 ym=w'*h;                      %辨识器网络输出ym
 
%%%%%%%输出权重更新计算%%%%%%%

 d_w=0*w;
   for j=1:1:6
      d_w(j)=xite*(u(4)-ym)*h(j);%计算deta_w供迭代计算下次w，权值更新
   end
   w=w_1+d_w+alfa*(w_1-w_2);        %权值向量w调整
   
  d_bi=0*bi;
   for j=1:1:6
      d_bi(j)=xite*(u(4)-ym)*w(j)*h(j)*(bi(j)^-3)*norm(RBF_input-ci(:,j))^2;%计算deta_b供迭代计算下次第i个节点的基宽向量
   end
   bi=bi_1+ d_bi+alfa*(bi_1-bi_2);  %基宽向量B调整
   
   for j=1:1:6
     for i=1:1:3
      d_ci(i,j)=xite*(u(4)-ym)*w(j)*h(j)*(RBF_input(i)-ci(i,j))*(bi(j)^-2);%计算deta_ci供迭代计算下次第i个节点的中心矢量ci
     end
   end
   ci=ci_1+d_ci+alfa*(ci_1-ci_2);   %中心矢量调整

 %%%%%%%系统雅可比矩阵计算%%%%%%%  
yu=0;
  for j=1:1:6
      yu=yu+w(j)*h(j)*(-RBF_input(1)+ci(1,j))/bi(j)^2;
  end
dyout=yu; %被控对象的Jacobian信息

   %用梯度下降法调整P,I,D
   kp_1=kp0+xitekp*u(1)*dyout*x(1);
   ki_1=ki0+xiteki*u(1)*dyout*x(2);
   kd_1=kd0+xitekd*u(1)*dyout*x(3);
   if kp_1<0
      kp_1<0;
   end
   if kd_1<0
      kd_1=0;
   end
   if ki_1<0
      ki_1=0;
   end
   
  
   
   deta_u=kp_1*(u(1)-u(2))+ki_1*u(1)+kd_1*(u(1)-2*u(2)+u(3));%PID输出的调整量deta_u
   uu=iq+deta_u;%计算下次的输出
   ci_3=ci_2;
   ci_2=ci_1;
   ci_1=ci;
   
   bi_3=bi_2;
   bi_2=bi_1;
   bi_1=bi;
   
   w_3=w_2;
   w_2=w_1;
   w_1=w;
   
   kp0=kp_1;
   kd0=kd_1;
   ki0=ki_1;          

sys = [uu;ym];

