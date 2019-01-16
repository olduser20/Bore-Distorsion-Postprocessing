tic

clc
clearvars

scale=500;


% Input: The section which user wants to study.
SecNo_CMM=11;


indx=1;
indy=2;

% There are 11 sections in each bore that has been measured by CMM.
% The ICE has 4 cylinders.

% CMM Open Deck with Torque Plate
CMM_OD_WT=cell(11,4);
% CMM Open Deck without Torque Plate
CMM_OD_WO=cell(11,4);
% CMM Close Deck with Torque Plate
CMM_CD_WT=cell(11,4);
% CMM Close Deck without Torque Plate
CMM_CD_WO=cell(11,4);

% FEM Open Deck
FEM_OD_WT=cell(41,4);
% FEM Close Deck
FEM_CD_WT=cell(41,4);

% Number of bores
n=4;

% Address of files which contain measured displacements.
% OD is short for "Open Deck".
% WT is short for "With Torque Plate"
% WO is short for "Without Torque Plate"

address_CMM_OD_WT='CMM_OD_WT\';
address_CMM_OD_WO='CMM_OD_WO\';
address_CMM_CD_WT='CMM_CD_WT\';
address_CMM_CD_WO='CMM_CD_WO\';

% Initializing an empty matrix which will contain maximum displacement in each case.
DIS=[];



for LinerNo=1:n
    
    displacement_names_CMM=['LINER0' int2str(LinerNo) 'SEC' strSecVec(SecNo_CMM*2-1:SecNo_CMM*2) '.txt'];
    displacement_names_FEM=['STEP1LINER0' int2str(LinerNo) 'SEC' strSecVec(SecNo_FEM*2-1:SecNo_FEM*2) '.txt'];
   
    
    displacement_CMM_CD_WT=importdata([address_CMM_CD_WT displacement_names_CMM(1,:)]);
    displacement_CMM_CD_WO=importdata([address_CMM_CD_WO displacement_names_CMM(1,:)]);
    displacement_CMM_OD_WT=importdata([address_CMM_OD_WT displacement_names_CMM(1,:)]);
    displacement_CMM_OD_WO=importdata([address_CMM_OD_WO displacement_names_CMM(1,:)]);
    
    
    nodes_FEM_OD_WT=importdata([address_FEM_OD_WT 'nodes.txt']);
    nodes_FEM_CD_WT=importdata([address_FEM_CD_WT 'nodes.txt']);
    displacement_FEM_OD_WT=importdata([address_FEM_OD_WT displacement_names_FEM(1,:)]);
    displacement_FEM_CD_WT=importdata([address_FEM_CD_WT displacement_names_FEM(1,:)]);
    
    temp_OD=node_extract(nodes_FEM_OD_WT,displacement_FEM_OD_WT);
    temp_CD=node_extract(nodes_FEM_CD_WT,displacement_FEM_CD_WT);
    
    
    N_CMM_CD_WT=length(displacement_CMM_CD_WT(:,1));
    N_CMM_CD_WO=length(displacement_CMM_CD_WO(:,1));
    N_CMM_OD_WT=length(displacement_CMM_OD_WT(:,1));
    N_CMM_OD_WO=length(displacement_CMM_OD_WO(:,1));
    
    N_FEM_OD_WT=length(displacement_FEM_OD_WT(:,1));
    N_FEM_CD_WT=length(displacement_FEM_CD_WT(:,1));
    
    


    C_CMM_CD_WT=mean(displacement_CMM_CD_WT(:,1:3));
    C_CMM_CD_WO=mean(displacement_CMM_CD_WO(:,1:3));
    C_CMM_OD_WT=mean(displacement_CMM_OD_WT(:,1:3));
    C_CMM_OD_WO=mean(displacement_CMM_OD_WO(:,1:3));
    
    C_FEM_OD_WT=mean(temp_OD(:,2:4));
    C_FEM_CD_WT=mean(temp_CD(:,2:4));
    
    
    R_CMM_CD_WT=norm(displacement_CMM_CD_WT(1,1:3)-C_CMM_CD_WT);
    R_CMM_CD_WO=norm(displacement_CMM_CD_WO(1,1:3)-C_CMM_CD_WO);
    R_CMM_OD_WT=norm(displacement_CMM_OD_WT(1,1:3)-C_CMM_OD_WT);
    R_CMM_OD_WO=norm(displacement_CMM_OD_WO(1,1:3)-C_CMM_OD_WO);
    
    R_FEM_OD_WT=norm(temp_OD(1,2:4)-C_FEM_OD_WT);
    R_FEM_CD_WT=norm(temp_CD(1,2:4)-C_FEM_CD_WT);
    
    

    [theta_CMM_CD_WT,r_CMM_CD_WT]=cart2pol(displacement_CMM_CD_WT(:,indx),displacement_CMM_CD_WT(:,indy));
    [theta_CMM_CD_WO,r_CMM_CD_WO]=cart2pol(displacement_CMM_CD_WO(:,indx),displacement_CMM_CD_WO(:,indy));
    [theta_CMM_OD_WT,r_CMM_OD_WT]=cart2pol(displacement_CMM_OD_WT(:,indx),displacement_CMM_OD_WT(:,indy));
    [theta_CMM_OD_WO,r_CMM_OD_WO]=cart2pol(displacement_CMM_OD_WO(:,indx),displacement_CMM_OD_WO(:,indy));
    
    [theta_FEM_OD_WT,r_FEM_OD_WT]=cart2pol(temp_OD(:,indx),temp_OD(:,indy));
    [theta_FEM_CD_WT,r_FEM_CD_WT]=cart2pol(temp_CD(:,indx),temp_CD(:,indy));
    
    

    
    dr_CMM_CD_WT=r_CMM_CD_WT-39.3;
    dr_CMM_CD_WO=r_CMM_CD_WO-39.3;
    dr_CMM_OD_WT=r_CMM_OD_WT-39.3;
    dr_CMM_OD_WO=r_CMM_OD_WO-39.3;
    
    temp_CMM_CD_WT=[displacement_CMM_CD_WT(:,1:3) theta_CMM_CD_WT r_CMM_CD_WT dr_CMM_CD_WT];
    temp_CMM_CD_WO=[displacement_CMM_CD_WO(:,1:3) theta_CMM_CD_WO r_CMM_CD_WO dr_CMM_CD_WO];
    temp_CMM_OD_WT=[displacement_CMM_OD_WT(:,1:3) theta_CMM_OD_WT r_CMM_OD_WT dr_CMM_OD_WT];
    temp_CMM_OD_WO=[displacement_CMM_OD_WO(:,1:3) theta_CMM_OD_WO r_CMM_OD_WO dr_CMM_OD_WO];
    
    
    temp_CMM_CD_WT=sortrows(temp_CMM_CD_WT,4);
    temp_CMM_CD_WTn=[temp_CMM_CD_WT;temp_CMM_CD_WT(1,:)];
    temp_CMM_CD_WO=sortrows(temp_CMM_CD_WO,4);
    temp_CMM_CD_WOn=[temp_CMM_CD_WO;temp_CMM_CD_WO(1,:)];
    temp_CMM_OD_WT=sortrows(temp_CMM_OD_WT,4);
    temp_CMM_OD_WTn=[temp_CMM_OD_WT;temp_CMM_OD_WT(1,:)];
    temp_CMM_OD_WO=sortrows(temp_CMM_OD_WO,4);
    temp_CMM_OD_WOn=[temp_CMM_OD_WO;temp_CMM_OD_WO(1,:)];
    
    
    
    CMM_CD_WT{SecNo_CMM,LinerNo}=temp_CMM_CD_WT;
    CMM_CD_WTn{SecNo_CMM,LinerNo}=temp_CMM_CD_WTn;
    CMM_CD_WO{SecNo_CMM,LinerNo}=temp_CMM_CD_WO;
    CMM_CD_WOn{SecNo_CMM,LinerNo}=temp_CMM_CD_WOn;
    CMM_OD_WT{SecNo_CMM,LinerNo}=temp_CMM_OD_WT;
    CMM_OD_WTn{SecNo_CMM,LinerNo}=temp_CMM_OD_WTn;
    CMM_OD_WO{SecNo_CMM,LinerNo}=temp_CMM_OD_WO;
    CMM_OD_WOn{SecNo_CMM,LinerNo}=temp_CMM_OD_WOn;
    
    

    clr=['b' 'g' 'r' 'm'];
    

    subplot(1,4,5-LinerNo)

    % Plotting deformed section for each case
    polar(CMM_CD_WT{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_CD_WT{SecNo_CMM,LinerNo}(:,6),'m');hold on
    polar(CMM_CD_WOn{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_CD_WOn{SecNo_CMM,LinerNo}(:,6),'c');hold on
    polar(CMM_OD_WT{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_OD_WT{SecNo_CMM,LinerNo}(:,6),'b');hold on
    polar(CMM_OD_WO{SecNo_CMM,LinerNo}(:,4),39.3+scale*CMM_OD_WO{SecNo_CMM,LinerNo}(:,6),'r')
    axis equal
    
    % Plotting the base circle
    polar(CMM_OD_WTn{SecNo_CMM,LinerNo}(:,4),CMM_OD_WTn{SecNo_CMM,LinerNo}(:,5),'--');hold on

    
    grid on


%     dr_CMM_CD_WT
%     dr_CMM_CD_WO
%     dr_CMM_OD_WT
%     dr_CMM_OD_WO
        DIS=[DIS;max(dr_CMM_CD_WO) max(dr_CMM_CD_WT) max(dr_CMM_OD_WO) max(dr_CMM_OD_WT)];

end

DIS


toc
