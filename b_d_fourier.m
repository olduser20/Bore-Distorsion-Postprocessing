tic

clc
clearvars

scale=500;


% LinerNo=2;
FrOrd=6;
SecNo_CMM=11;
SecNo_FEM=1;
strSecVec=['01' '02' '03' '04' '05' '06' '07' '08' '09' '10',...
           '11' '12' '13' '14' '15' '16' '17' '18' '19' '20',...
           '21' '22' '23' '24' '25' '26' '27' '28' '29' '30',...
           '31' '32' '33' '34' '35' '36' '37' '38' '39' '40',...
           '41' '42' '43' '44' '45' '46' '47' '48' '49' '50'];


% displacement_names_CMM=['LINER0' int2str(LinerNo) 'SEC' strSecVec(SecNo_FEM*2-1:SecNo_FEM*2) '.txt'];
% displacement_names_FEM=['STEP1LINER0' int2str(LinerNo) 'SEC' strSecVec(SecNo_FEM*2-1:SecNo_FEM*2) '.txt'];


indx=1;
indy=2;
% CMM Open Deck with Torque Plate
% CMM_OD_WT=cell(11,4);
% CMM Open Deck without Torque Plate
% CMM_OD_WO=cell(11,4);
% CMM Close Deck with Torque Plate
CMM_CD_WT=cell(11,4);
% CMM Close Deck without Torque Plate
CMM_CD_WO=cell(11,4);

% FEM Open Deck
% FEM_OD_WT=cell(41,4);
% FEM Close Deck
FEM_CD_WT=cell(41,4);


%Data=cell(3,1);
%Datan=cell(3,1);
%CMM_OD_WT=cell(3,1);
%CMM_OD_WO=cell(3,1);

n=4;

% diameters=[0.18 0.05 0.12 0.16
%            0.26 0.07 0.15 0.24];

%diameters=[0.08 0.01 0.00 0.05
%           0.12 0.02 0.02 0.10];

% address_CMM_OD_WT='CMM_OD_WT\';
% address_CMM_OD_WO='CMM_OD_WO\';
% address_CMM_CD_WT='CMM_CD_WT\';
% address_CMM_CD_WO='CMM_CD_WO\';

% address_FEM_OD_WT='FEM_OD_WT\';
address_FEM_CD_WT='FEM_CD_WT\';
DIS=[];

for LinerNo=1:n
    
    displacement_names_CMM=['LINER0' int2str(LinerNo) 'SEC' strSecVec(SecNo_CMM*2-1:SecNo_CMM*2) '.txt'];
    displacement_names_FEM=['STEP1LINER0' int2str(LinerNo) 'SEC' strSecVec(SecNo_FEM*2-1:SecNo_FEM*2) '.txt'];

    
    %nodes=importdata([address nodes_names(i,:)]);
%     displacement_CMM_CD_WT=importdata([address_CMM_CD_WT displacement_names_CMM(1,:)]);
%     displacement_CMM_CD_WO=importdata([address_CMM_CD_WO displacement_names_CMM(1,:)]);
%     displacement_CMM_OD_WT=importdata([address_CMM_OD_WT displacement_names_CMM(1,:)]);
%     displacement_CMM_OD_WO=importdata([address_CMM_OD_WO displacement_names_CMM(1,:)]);
    
    
%     nodes_FEM_OD_WT=importdata([address_FEM_OD_WT 'nodes.txt']);
    nodes_FEM_CD_WT=importdata([address_FEM_CD_WT 'nodes.txt']);
%     displacement_FEM_OD_WT=importdata([address_FEM_OD_WT displacement_names_FEM(1,:)]);
    displacement_FEM_CD_WT=importdata([address_FEM_CD_WT displacement_names_FEM(1,:)]);
    
%     temp_OD=node_extract(nodes_FEM_OD_WT,displacement_FEM_OD_WT);
    temp_CD=node_extract(nodes_FEM_CD_WT,displacement_FEM_CD_WT);
    
    
%     N_CMM_CD_WT=length(displacement_CMM_CD_WT(:,1));
%     N_CMM_CD_WO=length(displacement_CMM_CD_WO(:,1));
%     N_CMM_OD_WT=length(displacement_CMM_OD_WT(:,1));
%     N_CMM_OD_WO=length(displacement_CMM_OD_WO(:,1));
    
%     N_FEM_OD_WT=length(displacement_FEM_OD_WT(:,1));
    N_FEM_CD_WT=length(displacement_FEM_CD_WT(:,1));
    
    
%     diag_vec=[round(N/4+N/2) ; round(N/4)
%                 round(N/2)   ;  round(N) 
%             round(3*N/8+N/2) ; round(3*N/8)
%               round(N/8+N/2) ; round(N/8)];
    
    % Mirroring points w.r.t y direction
    %nodes(:,[4])=-nodes(:,[4]);
    %displacement(:,[4])=-displacement(:,[4]);

%     C_CMM_CD_WT=mean(displacement_CMM_CD_WT(:,1:3));
%     C_CMM_CD_WO=mean(displacement_CMM_CD_WO(:,1:3));
%     C_CMM_OD_WT=mean(displacement_CMM_OD_WT(:,1:3));
%     C_CMM_OD_WO=mean(displacement_CMM_OD_WO(:,1:3));
    
%     C_FEM_OD_WT=mean(temp_OD(:,2:4));
    C_FEM_CD_WT=mean(temp_CD(:,2:4));
    
    
%     R_CMM_CD_WT=norm(displacement_CMM_CD_WT(1,1:3)-C_CMM_CD_WT);
%     R_CMM_CD_WO=norm(displacement_CMM_CD_WO(1,1:3)-C_CMM_CD_WO);
%     R_CMM_OD_WT=norm(displacement_CMM_OD_WT(1,1:3)-C_CMM_OD_WT);
%     R_CMM_OD_WO=norm(displacement_CMM_OD_WO(1,1:3)-C_CMM_OD_WO);
    
%     R_FEM_OD_WT=norm(temp_OD(1,2:4)-C_FEM_OD_WT);
    R_FEM_CD_WT=norm(temp_CD(1,2:4)-C_FEM_CD_WT);
    
    
    % nodes(:,2:4)=nodes(:,2:4)-C.*ones(length(nodes(:,1)),3);
%     [theta_CMM_CD_WT,r_CMM_CD_WT]=cart2pol(displacement_CMM_CD_WT(:,indx),displacement_CMM_CD_WT(:,indy));
%     [theta_CMM_CD_WO,r_CMM_CD_WO]=cart2pol(displacement_CMM_CD_WO(:,indx),displacement_CMM_CD_WO(:,indy));
%     [theta_CMM_OD_WT,r_CMM_OD_WT]=cart2pol(displacement_CMM_OD_WT(:,indx),displacement_CMM_OD_WT(:,indy));
%     [theta_CMM_OD_WO,r_CMM_OD_WO]=cart2pol(displacement_CMM_OD_WO(:,indx),displacement_CMM_OD_WO(:,indy));
    
%     [theta_FEM_OD_WT,r_FEM_OD_WT]=cart2pol(temp_OD(:,indx),temp_OD(:,indy));
    [theta_FEM_CD_WT,r_FEM_CD_WT]=cart2pol(temp_CD(:,indx),temp_CD(:,indy));
    
    
    %displacement2=displacement(:,1:3)+scale*displacement(:,2:4);
    %[theta2,r2]=cart2pol(nodes2(:,indx-1),nodes2(:,indy-1));
    
    %dtheta=theta2-theta;
    %dr=r2-r;
    
    

    %nodes=[nodes theta r];
    
%     dr_CMM_CD_WT=r_CMM_CD_WT-39.3;
%     dr_CMM_CD_WO=r_CMM_CD_WO-39.3;
%     dr_CMM_OD_WT=r_CMM_OD_WT-39.3;
%     dr_CMM_OD_WO=r_CMM_OD_WO-39.3;
    
%     temp_CMM_CD_WT=[displacement_CMM_CD_WT(:,1:3) theta_CMM_CD_WT r_CMM_CD_WT dr_CMM_CD_WT];
%     temp_CMM_CD_WO=[displacement_CMM_CD_WO(:,1:3) theta_CMM_CD_WO r_CMM_CD_WO dr_CMM_CD_WO];
%     temp_CMM_OD_WT=[displacement_CMM_OD_WT(:,1:3) theta_CMM_OD_WT r_CMM_OD_WT dr_CMM_OD_WT];
%     temp_CMM_OD_WO=[displacement_CMM_OD_WO(:,1:3) theta_CMM_OD_WO r_CMM_OD_WO dr_CMM_OD_WO];
    
    
%     temp_CMM_CD_WT=sortrows(temp_CMM_CD_WT,4);
%     temp_CMM_CD_WTn=[temp_CMM_CD_WT;temp_CMM_CD_WT(1,:)];
%     temp_CMM_CD_WO=sortrows(temp_CMM_CD_WO,4);
%     temp_CMM_CD_WOn=[temp_CMM_CD_WO;temp_CMM_CD_WO(1,:)];
%     temp_CMM_OD_WT=sortrows(temp_CMM_OD_WT,4);
%     temp_CMM_OD_WTn=[temp_CMM_OD_WT;temp_CMM_OD_WT(1,:)];
%     temp_CMM_OD_WO=sortrows(temp_CMM_OD_WO,4);
%     temp_CMM_OD_WOn=[temp_CMM_OD_WO;temp_CMM_OD_WO(1,:)];
    
    
    
%     CMM_CD_WT{SecNo_CMM,LinerNo}=temp_CMM_CD_WT;
%     CMM_CD_WTn{SecNo_CMM,LinerNo}=temp_CMM_CD_WTn;
%     CMM_CD_WO{SecNo_CMM,LinerNo}=temp_CMM_CD_WO;
%     CMM_CD_WOn{SecNo_CMM,LinerNo}=temp_CMM_CD_WOn;
%     CMM_OD_WT{SecNo_CMM,LinerNo}=temp_CMM_OD_WT;
%     CMM_OD_WTn{SecNo_CMM,LinerNo}=temp_CMM_OD_WTn;
%     CMM_OD_WO{SecNo_CMM,LinerNo}=temp_CMM_OD_WO;
%     CMM_OD_WOn{SecNo_CMM,LinerNo}=temp_CMM_OD_WOn;
    
%     FEM_OD_WT{SecNo_FEM,LinerNo}=temp_OD;
    FEM_CD_WT{SecNo_FEM,LinerNo}=temp_CD;
    
    
    % Fourier series expansion


T0=2*pi;
w0=2*pi/T0;

y=zeros(length(r_FEM_CD_WT),1);
z=zeros(length(r_FEM_CD_WT),1);

for n=1:FrOrd
    y=0;

    for jj=1:length(r_FEM_CD_WT)
        y(jj,1)= r_FEM_CD_WT(jj) * cos((n-1)*w0*theta_FEM_CD_WT(jj));
        z(jj,1)= r_FEM_CD_WT(jj) * cos((n-1)*w0*theta_FEM_CD_WT(jj));
    end

    A(n,1) = (2/T0) * trapz (theta_FEM_CD_WT,y);
    B(n,1) = (2/T0) * trapz (theta_FEM_CD_WT,z);
    U(n,1) = 2000 * sqrt ((A(n)^2)+ (B(n)^2));

end
    
%     plot3(datan(:,2),datan(:,3),datan(:,4));hold on
%     plot(datan(:,2),datan(:,3),'g');hold on

    %polar(datan(:,5),datan(:,6),'y');hold on
%     polar()
    clr=['b' 'g' 'r' 'm'];
    
%     figure
    subplot(1,4,5-LinerNo)
%     polar(CMM_CD_WT{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_CD_WT{SecNo_CMM,LinerNo}(:,6),'m-');hold on
%     polar(CMM_CD_WT{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_CD_WT{SecNo_CMM,LinerNo}(:,6),'m');hold on
%     polar(CMM_CD_WT{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_CD_WT{SecNo_CMM,LinerNo}(:,6),'c');hold on
%     polar(CMM_OD_WT{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_OD_WT{SecNo_CMM,LinerNo}(:,6),'b');hold on
%     polar(CMM_OD_WO{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_OD_WO{SecNo_CMM,LinerNo}(:,6),'r')
    axis equal
    
%     figure
%     subplot(1,4,5-LinerNo)
    polar(FEM_CD_WT{SecNo_FEM,LinerNo}(:,2),FEM_CD_WT{SecNo_FEM,LinerNo}(:,3)+scale*FEM_CD_WT{SecNo_FEM,LinerNo}(:,3),'m');hold on
%     polar(FEM_OD_WT{SecNo_FEM,LinerNo}(:,4),FEM_OD_WT{SecNo_FEM,LinerNo}(:,5)+scale*FEM_OD_WT{SecNo_FEM,LinerNo}(:,6),'b');hold on
%     axis equal
%     
%     polar(CMM_CD_WTn{SecNo_CMM,LinerNo}(:,4),CMM_CD_WTn{SecNo_CMM,LinerNo}(:,5),'--');hold on
% %     text(datan(:,2),datan(:,3),num2str(datan(:,10)*180/pi))
%     plot3(C(1),C(2),C(3),'*m')
%     plot(C(1),C(2),'*m');hold on
    
%     polarplot(Datan{i}(diag_vec(:),10),Datan{i}(diag_vec(:),11),...
%         'marker','o','color',clr(i),'linestyle','none')
%     plot3(datan(:,2)+scale*datan(:,7),...
%           datan(:,3)+scale*datan(:,8),...
%           datan(:,4)+scale*datan(:,9),...
%           'color','m');
%     plot(datan(:,2)+scale*datan(:,7),...
%          datan(:,3)+scale*datan(:,8),...
%          'color','b');
%         axis equal
        grid on
%     
%     d_vertical=norm([Data{i}(round(N/4+N/2),7)-Data{i}(round(N/4),7),...
%                  Data{i}(round(N/4+N/2),8)-Data{i}(round(N/4),8),...
%                  Data{i}(round(N/4+N/2),9)-Data{i}(round(N/4),9)])*2
% 
% 
%     d_horizontal=norm([Data{i}(round(N/2),7)-Data{i}(N,7),...
%                        Data{i}(round(N/2),8)-Data{i}(N,8),...
%                        Data{i}(round(N/2),9)-Data{i}(N,9)])*2
% 
%     d_q24=norm([Data{i}(round(3*N/8+N/2),7)-Data{i}(round(3*N/8),7),...
%                 Data{i}(round(3*N/8+N/2),8)-Data{i}(round(3*N/8),8),...
%                 Data{i}(round(3*N/8+N/2),9)-Data{i}(round(3*N/8),9)])*2
% 
%     d_q13=norm([Data{i}(round(N/8+N/2),7)-Data{i}(round(N/8),7),...
%                 Data{i}(round(N/8+N/2),8)-Data{i}(round(N/8),8),...
%                 Data{i}(round(N/8+N/2),9)-Data{i}(round(N/8),9)])*2
%             
%     diameters=[d_vertical d_horizontal d_q24 d_q13];
%     dr_CMM_CD_WT
%     dr_CMM_CD_WO
%     dr_CMM_OD_WT
%     dr_CMM_OD_WO
%         DIS=[DIS;max(dr_CMM_CD_WO) max(dr_CMM_CD_WT) max(dr_CMM_OD_WO) max(dr_CMM_OD_WT)]

end


% figure(8)
%     stem(nodes(:,1)-displacement(:,1))
% figure
% stem(nodes(:,1)-displacement(:,1))

% figure(3)

% d_vertical=norm([Data{n}(round(N/4+N/2),7)-Data{n}(round(N/4),7),...
%                  Data{n}(round(N/4+N/2),8)-Data{n}(round(N/4),8),...
%                  Data{n}(round(N/4+N/2),9)-Data{n}(round(N/4),9)])*2
% 
% 
% d_horizontal=norm([Data{n}(round(N/2),7)-Data{n}(N,7),...
%                    Data{n}(round(N/2),8)-Data{n}(N,8),...
%                    Data{n}(round(N/2),9)-Data{n}(N,9)])*2
%                
% d_q24=norm([Data{n}(round(3*N/8+N/2),7)-Data{n}(round(3*N/8),7),...
%             Data{n}(round(3*N/8+N/2),8)-Data{n}(round(3*N/8),8),...
%             Data{n}(round(3*N/8+N/2),9)-Data{n}(round(3*N/8),9)])*2
%         
% d_q13=norm([Data{n}(round(N/8+N/2),7)-Data{n}(round(N/8),7),...
%             Data{n}(round(N/8+N/2),8)-Data{n}(round(N/8),8),...
%             Data{n}(round(N/8+N/2),9)-Data{n}(round(N/8),9)])*2

% C_orig=Datan{n}(:,6).*exp(Datan{n}(:,5)*1i);
%C_def=CMM_OD_WTn{i,1}(:,4).*exp(CMM_OD_WTn{i,1}(:,5)*1i);
%C_def_WO=Datan_WO{n}(:,4).*exp(Datan_WO{n}(:,5)*1i);
% plot(Data{n}(106:115,5)*180/pi,Data{n}(106:115,6),'.',Data{n}(106:115,10)*180/pi,Data{n}(106:115,11),'o')

% diameters_manual=[0.013 0.005 0.011 0.008;
%                   0.012 0.006 0.010 0.010];

% for q=2:2
% %     plot(diameters(q,:));hold on
%     plot(diameters_manual(q,:));hold on
% end
% ylim([0 0.15]);


toc