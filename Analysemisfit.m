% About the program%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Nian Wang & Yang Shen, 06/16/2019
% This is an example about the 3D full waveform location method
% used to invert a known explosive source - the Non 
% Proliferation Experiment at the Nevada Test Site in year 1993.
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
% The following example shows how to locate the best solution 
% and estimate the uncertainty using the misfit of P and P coda 
% waves (Figure 5b in Wang et al. (2019)) in the 3D grid search % region 

% About the misfit array data
% At each grid point, the misfit array is the L2 norm misfit 
% between the data and synthetic seismogram. 
% The synthetic seismogram is calculated based on the 
% convolution between the strain Green’s tensor (SGT) and the 
% equivalent explosive source moment tensor.
% The user could use any open source code to calculate the SGT 
% and get the synthetic seismograms


% The following is the detail of the program

% 1. Load the misfit array
  clear all;
  MGT=load('MisfitPcoda.mat');
  SumGT=MGT.SumGT;
  [n1 n2 n3]=size(SumGT)

% 2. Find the best solution in the search region
    
   % Scale/normalize the misfit (the best location will not be changed)
   
   gmin=min(min(min(SumGT)));
   gmax=max(max(max(SumGT)));
   SumGT=SumGT/gmax;
   
   % Get the minimum and maximum values and their indexes in the search region

   gmin=min(min(min(SumGT)));
   [i1,i2,i3]=ind2sub(size(SumGT),find(SumGT==gmin));

   gmax=max(max(max(SumGT)));
   [j1,j2,j3]=ind2sub(size(SumGT),find(SumGT==gmax));

      
   % Load the coordinates of the search region

   Cordx0=load('./gridsearchx.mat');
   Cordx=Cordx0.x1;
   Cordx1=Cordx(1:n1,1:n2,1:n3)/1000; % change unit from m to km

   Cordy0=load('./gridsearchy.mat');
   Cordy=Cordy0.y1;
   Cordy1=Cordy(1:n1,1:n2,1:n3)/1000;

   Cordz0=load('./gridsearchz.mat');
   Cordz=Cordz0.z1;
   Cordz1=Cordz(1:n1,1:n2,1:n3)/1000;


   % Get the location of the best solution (with the minimum misfit)
   invind1=i1;
   invind2=i2;
   invind3=i3;  
   Inv_loc=[Cordx1(invind1,invind2,invind3),Cordy1(invind1,invind2,invind3),Cordz1(invind1,invind2,invind3)]

  % Provide the location of the true solution
   ref2ind1=20;
   ref2ind2=22;
   ref2ind3=22;
   Ref_loc=[Cordx1(ref2ind1,ref2ind2,ref2ind3),Cordy1(ref2ind1,ref2ind2,ref2ind3),Cordz1(ref2ind1,ref2ind2,ref2ind3)]

  % Calculate the distance between the true location and the best solution
   dist1=abs(Ref_loc-Inv_loc)


% 3. Analyse the uncertainty of the best solution
    
 
   fvalue=0.125;  % Fisher 1 sigma 0.68, N=33 data, and (3/(N-3))*F(3,N-3)=1.25*3/(N-3)=0.125 
                  % Please see Wang et al. (2016, 2019) for details 

   %To find grid points in the 1 sigma region   
   %To avoid a single point far away from the clustered points, 
     % add constrains on the searching length for the confidence level      
     % constrain in confx*0.25=2.5 km on each side of the best solution
     % check if the final uncertainty will be much smaller than the constraints
    confx=10;
    confy=10;
    confz=10;
    index=1;
    for i=1:n1
    for j=1:n2
    for k=1:n3
    % if SumGT(i,j,k)<=gmin*(1+fvalue) && SumGT(i,j,k)>=gmin*(1-fvalue)
    if (SumGT(i,j,k)<=gmin*(1+fvalue)) && (SumGT(i,j,k)>=gmin*(1-fvalue) && ...
      (i<=(i1+confx)) && (i>=(i1-confx)) && ...
      (j<=(i2+confy)) && (j>=(i2-confy)) && ...
      (k<=(i3+confz)) && (k>=(i3-confz)) )    
    AA(index,1)=i;
    AA(index,2)=j;
    AA(index,3)=k;
    index=index+1;
    end
    end
    end
    end

   [m1,m2]=size(AA); 
   for i=1:m1
   BB(i,1)=Cordx1(AA(i,1),AA(i,2),AA(i,3));
   BB(i,2)=Cordy1(AA(i,1),AA(i,2),AA(i,3));
   BB(i,3)=Cordz1(AA(i,1),AA(i,2),AA(i,3));
   end
   BB=double(BB);
 

   %Estimate the uncertainty
   solmin=[Cordx1(i1,i2,i3) Cordy1(i1,i2,i3) Cordz1(i1,i2,i3)]
   data1=BB;

    %Use the mean distance between the best solution and all the possible points
    reso(1,1)=mean(abs(data1(:,1)-solmin(1)));
    reso(1,2)=mean(abs(data1(:,2)-solmin(2)));
    reso(1,3)=mean(abs(data1(:,3)-solmin(3)));
    reso=double(reso);
    Uncertainty=reso % in km
    save Uncertainty68.dat reso -ascii

   % Other (reasonable) measurements could be also used to estimate the uncertainty  

   %% Measure the uncertainty along the best solution
    
    %Along the z-direction 
    %[n1 n2]=size(data1)
    %j=1;
    %for k=1:n1
    %if (data1(k,1) == solmin(1)) && (data1(k,2) == solmin(2)) 
    %  dataz(j,:)=data1(k,:);
    %  j=j+1;
    %end
    %end
    %tmp=(max(dataz(:,3))-min(dataz(:,3)))/2;
    %reso2(1,3)=tmp;
    %reso3(1,3)=mean(abs(dataz(:,3)-solmin(3)));

    %Along the x-direction 
    %[n1 n2]=size(data1)
    %j=1;
    %for k=1:n1
    %if (data1(k,2) == solmin(2)) 
    %  datax(j,:)=data1(k,:);
    %  j=j+1;
    %end
    %end
    %tmp=(max(datax(:,1))-min(datax(:,1)))/2;
    %reso2(1,1)=tmp;
    %reso3(1,1)=mean(abs(datax(:,1)-solmin(1)));

    %Along the y-direction
    %[n1 n2]=size(data1)
    %j=1;
    %for k=1:n1
    %if (data1(k,1) == solmin(1)) 
    %  datay(j,:)=data1(k,:);
    %  j=j+1;
    %end
    %end
    %tmp=(max(datay(:,2))-min(datay(:,2)))/2;
    %reso2(1,2)=tmp;
    %reso3(1,2)=mean(abs(datay(:,2)-solmin(2)));
 
    %reso2
    %reso3
 
   



   %%Analyse the upper and lower bounds of the possible solutions
     %testres(1,:)=[min(data1(:,1)-solmin(1)) max(data1(:,1)-solmin(1))];
     %testres(2,:)=[min(data1(:,2)-solmin(2)) max(data1(:,2)-solmin(2))];
     %testres(3,:)=[min(data1(:,3)-solmin(3)) max(data1(:,3)-solmin(3))];

     %testres

