tic

clc
clearvars

scale=0;


% LinerNo=2;
FrOrd=6;
% SecNo_FEM=30;
strSecVec=['01' '02' '03' '04' '05' '06' '07' '08' '09' '10',...
           '11' '12' '13' '14' '15' '16' '17' '18' '19' '20',...
           '21' '22' '23' '24' '25' '26' '27' '28' '29' '30',...
           '31' '32' '33' '34' '35' '36' '37' '38' '39' '40',...
           '41' '42' '43' '44' '45' '46' '47' '48' '49' '50'];


indx=1;
indy=2;


% FEM Close Deck
FEM_CD_WT=cell(41,4);
FEM_CD_WT_Fourier=cell(41,4);

% Number of liners
n=4;

address_FEM_CD_WT='FEM_CD_WT\';

for LinerNo=1:n
    for SecNo_FEM=1:41
    
        displacement_names_FEM=['STEP1LINER0' int2str(LinerNo) 'SEC' strSecVec(SecNo_FEM*2-1:SecNo_FEM*2) '.txt'];
        nodes_FEM_CD_WT=importdata([address_FEM_CD_WT 'nodes.txt']);


        displacement_FEM_CD_WT=importdata([address_FEM_CD_WT displacement_names_FEM(1,:)]);

        temp_FEM_CD_WT=node_extract(nodes_FEM_CD_WT,displacement_FEM_CD_WT);
        temp_FEM_CD_WT=[temp_FEM_CD_WT displacement_FEM_CD_WT(:,2:4)];

        % Number of nodes on the section
        N_FEM_CD_WT=length(displacement_FEM_CD_WT(:,1));



        % Mirroring points w.r.t y direction
        %nodes(:,[4])=-nodes(:,[4]);
        %displacement(:,[4])=-displacement(:,[4]);

        % Calculating center of the circular section (exact geometry)
        C_FEM_CD_WT=mean(temp_FEM_CD_WT(:,2:4));

        % Calculating radius of the circular section (exact geometry)
        % For EF7 engine should be 39.300
        R_FEM_CD_WT=norm(temp_FEM_CD_WT(1,2:4)-C_FEM_CD_WT);


        % Translating the center of cordinate system to the section center 
        temp_FEM_CD_WT(:,(indx+1):(indy+1))=temp_FEM_CD_WT(:,(indx+1):(indy+1))-C_FEM_CD_WT(:,indx:indy);


        % Transforming to polar coordinate system
        [theta_FEM_CD_WT,r_FEM_CD_WT]=cart2pol(temp_FEM_CD_WT(:,indx+1),temp_FEM_CD_WT(:,indy+1));

        % Calculating the magnitude of displacement of each node
        dr_FEM_CD_WT=sqrt(temp_FEM_CD_WT(:,indx+4).^2+temp_FEM_CD_WT(:,indy+4).^2);

        temp_FEM_CD_WT=[temp_FEM_CD_WT theta_FEM_CD_WT r_FEM_CD_WT dr_FEM_CD_WT];



        temp_FEM_CD_WT=sortrows(temp_FEM_CD_WT,8);

        FEM_CD_WT{SecNo_FEM,LinerNo}=temp_FEM_CD_WT;


        % Fourier series expansion


        T0=2*pi;
        w0=2*pi/T0;

        x=zeros(length(dr_FEM_CD_WT),1);
        y=zeros(length(dr_FEM_CD_WT),1);
        u=zeros(length(dr_FEM_CD_WT),1);
        phi=zeros(length(dr_FEM_CD_WT),1);

        for n=0:FrOrd
        %     y=0;

            for jj=1:length(dr_FEM_CD_WT)
                x(jj,1)= dr_FEM_CD_WT(jj) * cos((n)*w0*theta_FEM_CD_WT(jj));
                y(jj,1)= dr_FEM_CD_WT(jj) * sin((n)*w0*theta_FEM_CD_WT(jj));
        %         u(jj,1)=x(jj,1)+y(jj,1)*1j;
            end

            A(n+1,1) = (2/T0) * trapz (theta_FEM_CD_WT,x);
            B(n+1,1) = (2/T0) * trapz (theta_FEM_CD_WT,y);
        %     Un(n+1,1) = 2000 * sqrt ((A(n+1)^2)+ (B(n+1)^2));
            U(n+1,1) = A(n+1)+1j*B(n+1);
            Uabs(n+1,1)=2000* abs(U(n+1,1));
            Phi(n+1,1)=angle(U(n+1,1))*180/pi;


        end

        FEM_CD_WT_Fourier{SecNo_FEM,LinerNo}=[[0:FrOrd]' Uabs Phi];

    %     plot3(datan(:,2),datan(:,3),datan(:,4));hold on
    %     plot(datan(:,2),datan(:,3),'g');hold on

        %polar(datan(:,5),datan(:,6),'y');hold on
    %     polar()
        clr=['b' 'g' 'r' 'm'];

    %     figure
        subplot(1,4,5-LinerNo)
        axis equal

    %     figure
    %     subplot(1,4,5-LinerNo)
%         polar(FEM_CD_WT{SecNo_FEM,LinerNo}(:,8),FEM_CD_WT{SecNo_FEM,LinerNo}(:,9)+scale*FEM_CD_WT{SecNo_FEM,LinerNo}(:,9),'m');hold on
    %     polar(FEM_OD_WT{SecNo_FEM,LinerNo}(:,4),FEM_OD_WT{SecNo_FEM,LinerNo}(:,5)+scale*FEM_OD_WT{SecNo_FEM,LinerNo}(:,6),'b');hold on
    %     axis equal
    %     
    %     polar(CMM_CD_WTn{SecNo_CMM,LinerNo}(:,4),CMM_CD_WTn{SecNo_CMM,LinerNo}(:,5),'--');hold on
    % %     text(datan(:,2),datan(:,3),num2str(datan(:,10)*180/pi))
    %     plot3(C(1),C(2),C(3),'*m')
    %     plot(C(1),C(2),'*m');hold on

        grid on
    %     
    end
end


% figure(8)
%     stem(nodes(:,1)-displacement(:,1))
% figure
% stem(nodes(:,1)-displacement(:,1))



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