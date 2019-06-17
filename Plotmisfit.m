% About the program%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Nian Wang & Yang Shen, 06/16/2019
% This is an example about how to plot the misfit results of the  
% 3D full waveform location method used to invert a known 
% explosive source - the Non Proliferation Experiment at the 
% Nevada Test Site in year 1993.
% Refer to Wang et al. 2019 JGR (under minor revision) 
% for details 

% About the citation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you are using (part of) the code, please cite at least one % of the following papers:
% 1. N., Wang, Y., Shen, A., Flinders, W., Zhang (2016).
%    Accurate source location from waves scattered by surface 
%    topography. J. Geophys. Res., 121, 4538-4552.
% 2. N., Wang, Y., Shen, X. Y., Bao, A., Flinders (2019). 
%    Locating shallow seismic sources with waves scattered by 
%    surface topography: validation of the method at the Nevada 
%    Test Sites. J. Geophys. Res. (under minor revision).

% About the details of the program%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following example shows how to make the contour plots 
% (Figure 5b in Wang et al. (2019)) based on the misfit of P and % P coda waves in the 3D grid search region 


% The following is the detail of the program

% 1. Load the misfit array in the grid search region
     clear all;
     MGT=load('MisfitPcoda.mat')
     SumGT=MGT.SumGT;
     [n1 n2 n3]=size(SumGT)

   % Find the minimum and maximum value
     gmin=min(min(min(SumGT)));
     [i1,i2,i3]=ind2sub(size(SumGT),find(SumGT==gmin));
     gmax=max(max(max(SumGT)));
     [j1,j2,j3]=ind2sub(size(SumGT),find(SumGT==gmax));
 
     SumGT=SumGT/gmax;

   % Load the coordinates in the grid search region
     Cordx0=load('./gridsearchx.mat');
     Cordx=Cordx0.x1;
     Cordx1=Cordx/1000;
 
     Cordy0=load('./gridsearchy.mat');
     Cordy=Cordy0.y1;
     Cordy1=Cordy/1000;

     Cordz0=load('./gridsearchz.mat');
     Cordz=Cordz0.z1;
     Cordz1=Cordz/1000;

   % Get the inverted solution and true solution's coordinates 

     invind1=i1;
     invind2=i2;
     invind3=i3;
     Inv_loc=[Cordx1(invind1,invind2,invind3),Cordy1(invind1,invind2,invind3),Cordz1(invind1,invind2,invind3)]

     refind1=20;
     refind2=22;
     refind3=22;
     Ref_loc=[Cordx1(refind1,refind2,refind3),Cordy1(refind1,refind2,refind3),Cordz1(refind1,refind2,refind3)]

     dist=abs(Ref_loc-Inv_loc)
 
 
%2. Contour plot cut at the best solution in the y-z plane

     X=squeeze(Cordy1(i1,:,:));
     Y=squeeze(Cordz1(i1,:,:));
     AA=squeeze(SumGT(i1,:,:));
     v=[min(min(AA)):min(min(AA))/5:2*min(min(AA)) 2*min(min(AA)):max(max(AA))/6:max(max(AA))];

     figure;
     [C,h]=contourf(X,Y,AA,v); 
     set(gca,'Fontsize',40,'Linewidth',3);                                           
     % xlabel('Y-grid','Fontsize',40);
     % ylabel('Z-grid','Fontsize',40);
     colormap('jet');
     colorbar('linewidth',3);
     caxis([0.0 0.2]);
     hcb=colorbar;
     set(hcb,'YTick',[0 0.1 0.2],'Linewidth',4,'Fontsize',40);
     set(hcb,'YTick',[0 0.1 0.2],'Linewidth',4,'Fontsize',40);   
     scl_daspect=[1 1 1];
     daspect(scl_daspect)
     hold on
     scatter(Ref_loc(2),Ref_loc(3),700,'k','o','Linewidth',3)
     scatter(Inv_loc(2),Inv_loc(3),700,'k','v','Linewidth',3)
     ylim([-7.5 2.4]);
     % saveas(gcf,'GTmisfityz.fig','fig')


%3. Contour plot cut at the best solution in the x-z plane

     X=squeeze(Cordx1(:,i2,:));
     Y=squeeze(Cordz1(:,i2,:));
     AA=squeeze(SumGT(:,i2,:));
     v=[min(min(AA)):min(min(AA))/3:2*min(min(AA)) 2*min(min(AA)):max(max(AA))/10:max(max(AA))];
     figure;
     [C,h]=contourf(X,Y,AA,v);  
     set(gca,'Fontsize',40,'Linewidth',3);                                          
     % xlabel('X-grid','Fontsize',40);
     % ylabel('Z-grid','Fontsize',40);
     colormap('jet');
     colorbar('linewidth',3);
     caxis([0 0.2]);
     hcb=colorbar;
     set(hcb,'YTick',[0 0.1 0.2],'Linewidth',4,'Fontsize',40);
     set(hcb,'YTick',[0 0.1 0.2],'Linewidth',4,'Fontsize',40);
     scl_daspect=[1 1 1];
     daspect(scl_daspect)
     hold on
     scatter(Ref_loc(1),Ref_loc(3),700,'k','o','Linewidth',3)
     scatter(Inv_loc(1),Inv_loc(3),700,'k','v','Linewidth',3)
     ylim([-7.5 2.4]);
     % saveas(gcf,'GTmisfitxz.fig','fig')


%4. Contour plot cut at the best solution in the x-y plane  
 
     X=squeeze(Cordx1(:,:,i3));
     Y=squeeze(Cordy1(:,:,i3));
     AA=squeeze(SumGT(:,:,i3));   
     v=[min(min(AA)):min(min(AA))/4:1.05*min(min(AA)) 1.05*min(min(AA)):min(min(AA))/3.5:3*min(min(AA)) 3*min(min(AA)):max(max(AA))/5:max(max(AA))];
     figure;
     [C,h]=contourf(X,Y,AA,v);  
     set(gca,'Fontsize',40,'Linewidth',4);
     % xlabel('X-grid','Fontsize',40);
     % ylabel('Y-grid','Fontsize',40);
     colormap('jet');
     colorbar('linewidth',3);
     caxis([0 0.2]);
     hcb=colorbar;
     set(hcb,'YTick',[0 0.1 0.2],'Linewidth',4,'Fontsize',40);
     set(hcb,'YTick',[0 0.1 0.2],'Linewidth',4,'Fontsize',40);
     scl_daspect=[1 1 1];
     daspect(scl_daspect)
     hold on
     scatter(Ref_loc(1),Ref_loc(2),700,'k','o','Linewidth',3)
     scatter(Inv_loc(1),Inv_loc(2),700,'k','v','Linewidth',3)
     xlim([104.5 113.8]);
     set(gca,'XTick',[106 108 110 112])
     set(gca,'XTickLabel',{'106','108','110','112'})
     % saveas(gcf,'GTmisfitxy.fig','fig')

    
  
  
