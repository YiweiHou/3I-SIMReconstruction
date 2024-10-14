% Wirtten by Yiwei Hou, Xin Chen as part of the paper 
% "Triangle-beam interference structured illumination microscopy"
% Authority: Yiwei Hou, Xin Chen and Peng Xi
% College of fzuture technology, Peking University
% For any questions, please contact: houyiwei@stu.pku.edu.cn; xipeng@pku.edu.cn
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function varargout = Recon3ISIM(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Recon3ISIM_OpeningFcn, ...
                   'gui_OutputFcn',  @Recon3ISIM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
 
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
 
 
% --- Executes just before Recon3ISIM is made visible.
function Recon3ISIM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Recon3ISIM (see VARARGIN)
clear param
axes(handles.axes1); imshow(zeros(100,100),[])
axes(handles.axes2); imshow(zeros(100,100),[])
% Choose default command line output for Recon3ISIM
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
 
% UIWAIT makes Recon3ISIM wait for user response (see UIRESUME)
% uiwait(handles.figure1);
 
 
% --- Outputs from this function are returned to the command line.
function varargout = Recon3ISIM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Get default command line output from handles structure
varargout{1} = handles.output;
 
 
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param
%% Load raw SIM data
FlagLoad = 1;
[filename, pathname] = uigetfile({'*.tiff;*.tif'},'Select the raw data sequence','MultiSelect','off'); 
if isequal(filename,0)
    disp('Cancel data reading');
    return;
else 
    info = imfinfo(fullfile(pathname, filename));
    frames = numel(info);
  
    raw_data = zeros(info(1).Height, info(1).Width, frames, 'uint16');
     for k = 1:7
         raw_data(:,:,k) = im2uint16(imread(fullfile(pathname, filename), k));
     end
end
param.nz = frames/7; % Z-stack reconstruction or t-stack 7frames reconstruction
param.dz = 7;
Format = info.Format;
param.Format = Format;
param.filename = filename;
param.pathname = pathname;
FlagParameter = 0;                    
param.FlagParameter = FlagParameter;
param.FlagLoad = FlagLoad;
handles.param = param;
guidata(hObject,handles);
axes(handles.axes1);
I=0;
for i=1:7
    I=I+double(raw_data(:,:,i));
end
I=I/7;
param.I=I;
imshow(double(I),[]);
 
 
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param
global Mode
global f_stack_old 
addpath('.\functions')
tic
%% Step0: Pre-process
        param.NA = str2double(get(handles.NA,'String'));
        param.lambda = str2double(get(handles.wavelength,'String'));
        param.pixelsize = str2double(get(handles.pixelsize,'String'))*10^(-3);
        pathname=param.pathname;
        filename=param.filename;
        nz = param.nz;
        dz = param.dz;
% check the data size, if <256, force it to be 256
    [Np1,Np2] = size(param.I); 
    NPixel = max(Np1,Np2);
    if NPixel<256
       NPixel = 256;
    end
    
    parainheri=get(handles.Parainher,'value');  % whethre to overlap all time-lapse 
    % if overlap all images, generate the overlapped Recon3I image stack
    hb = waitbar(0,'doing global parameter estimation');
    if parainheri==1
        fprintf('doing global parameter estimation')
        raw_data=zeros(Np1,Np2,min(7*nz,70));
                     for k = 1:min(7*nz,70)
                         raw_data(:,:,k) = im2uint16(imread(fullfile(pathname, filename), k));
                         raw_data=double(raw_data);
                     end
             
        Hex_for_estimate=zeros(Np1,Np2,7);
        rounder=7;
        [~,~,d3]=size(raw_data);
        frame=d3/rounder;
        for i=1:rounder
            for j=1:frame
            Hex_for_estimate(:,:,i)=Hex_for_estimate(:,:,i)+raw_data(:,:,i+(j-1)*rounder);
            Hex_for_estimate(:,:,i)=Hex_for_estimate(:,:,i)/frame;
            ma=max(max(Hex_for_estimate(:,:,i)));
            mi=min(min(Hex_for_estimate(:,:,i)));
            Hex_for_estimate(:,:,i)=(2^16-1)*(Hex_for_estimate(:,:,i)-mi)/(ma-mi);
            end
        end
         clear raw_data 
         
        % do parameter estimation
        % Generate approximate OTF/PSF
   
        [OTFo, PSFe, Kotf, param] = Generate_OTF(Hex_for_estimate,param, NPixel);

        psf = fftshift(ifft2(ifftshift(OTFo)));
        Snoisy = deconvlucy(Hex_for_estimate,psf,5);
        %Snoisy = Hex_for_estimate;
        
        % parameter estimation
        s0_class{1} = single_orientation_class(OTFo,Kotf,Snoisy,PSFe);
        
        % frequency vector estimation
        kAo_L = zeros(3,2); kAmean_L = zeros(3,2); K_ini_L= zeros(1,3);
        kAo_H = zeros(3,2); kAmean_H = zeros(3,2); K_ini_H= zeros(1,3);
        flag = 1;
        s0_class{1} = KapproxEstimationF(s0_class{1},flag); % kAo contains three values along three directions
        s0_class{1}.kAmean = kmeanF(s0_class{1});  % kAmean contains three values
        kAo_L = s0_class{1}.kAo;
        kAmean_L = s0_class{1}.kAmean; 
        TE = kAmean_L;
        K_ini_L = sqrt(TE(:,1).^2+TE(:,2).^2);
        
        if max(K_ini_L)-min(K_ini_L)>1.5 
           flag = 0;         
           s0_class{1} = KapproxEstimationF(s0_class{1},flag);  
           s0_class{1}.kAmean = kmeanF(s0_class{1});
           kAo_H = s0_class{1}.kAo;
           kAmean_H = s0_class{1}.kAmean; 
           TE1 = kAmean_H;
           K_ini_H = sqrt(TE1(:,1).^2+TE1(:,2).^2);
           
           if max(K_ini_H)-min(K_ini_H)>1.5
                T_L = max(K_ini_L)-min(K_ini_L);
                T_H = max(K_ini_H)-min(K_ini_H);
                fprintf('The difference between K values is too great! The result may be not good!')
                if T_L > T_H
                    s0_class{1}.kAo = kAo_H;
                    s0_class{1}.kAmean = kAmean_H;
                else
                    s0_class{1}.kAo = kAo_L;
                    s0_class{1}.kAmean = kAmean_L;
                end
           end
        end
        % phase estimation
        s0_class{1}.phaseA = phasesF(s0_class{1});
        % modulation estimation
        [s0_class{1}.fDo,s0_class{1}.fDp,s0_class{1}.fDm,s0_class{1}.fDp1,s0_class{1}.fDm1,s0_class{1}.fDp2,s0_class{1}.fDm2] = PCMseparateF(s0_class{1});  
        s0_class{1}.modFac = ModulationFactorF(s0_class{1}); 
        fprintf('global parameter estimation complete')
             s=sprintf('doing global parameter estimation:%d',ceil(100));
             waitbar(1,hb,[s '%']);
             kA1 = s0_class{1}.kAmean(1,:);
             kB1 = s0_class{1}.kAmean(2,:);
             kC1 = s0_class{1}.kAmean(3,:);
             K1_dis=num2str(sqrt(kA1*kA1')); K1_dis=K1_dis(1:6);
             K2_dis=num2str(sqrt(kB1*kB1')); K2_dis=K2_dis(1:6);
             K3_dis=num2str(sqrt(kC1*kC1')); K3_dis=K3_dis(1:6);
             set(handles.k1,'String',K1_dis);
             set(handles.k2,'String',K2_dis);
             set(handles.k3,'String',K3_dis);
             if s0_class{1}.modFac(1)>1
                s0_class{1}.modFac(1)=0.01;
            end
             if s0_class{1}.modFac(2)>1
                s0_class{1}.modFac(2)=0.01;
             end
             if s0_class{1}.modFac(3)>1
                s0_class{1}.modFac(3)=0.01;
             end 
             angle1=num2str(rad2deg(atan(kA1(2)/kA1(1)))); angle1=angle1(1:6);
             angle2=num2str(rad2deg(atan(kB1(2)/kB1(1)))); angle2=angle2(1:6);
             angle3=num2str(rad2deg(atan(kC1(2)/kC1(1)))); angle3=angle3(1:6);
             set(handles.A1,'String',angle1);
             set(handles.A2,'String',angle2);
             set(handles.A3,'String',angle3);
    end
    close(hb);
    
    Rolling=get(handles.Rolling,'value');  % whethre to overlap all time-lapse 
    if Rolling==1
        nz=param.nz*7-6;
        param.nz=nz;
    end
    

 
    MRAreg=0.01;              %Modify according to requirement
    DataDir=param.pathname;
    mkdir([DataDir,'\','new3I_SIM']);
ha = waitbar(0,'Recon3I-SIM reconstruction processing');
filenum=1;
HexF=get(handles.HexF,'value');  % whethre to overlap all time-lapse 
    for k = 1:nz        
 % resotre the 7 raw Hex images
        Snoisy0 = zeros(Np1, Np2, 7);
        Snoisy = zeros(NPixel,NPixel,7);
        K_h = [Np1,Np2]; N_h = [NPixel, NPixel];
        L_h = ceil((N_h - K_h)/2); 
        v_h = colonvec(L_h+1,L_h+K_h);
        if Rolling==0
            for j = 1:7
                Snoisy0(:,:,j) = double(im2uint16(imread(fullfile(pathname, filename), (k-1)*dz+j)));
                Snoisy(v_h{:},j) = Snoisy0(:,:,j);
            end
        else
            for j = 1:7
                Snoisy0(:,:,j) =double(im2uint16(imread(fullfile(pathname, filename), (k-1)+j)));
                Snoisy(v_h{:},j) = Snoisy0(:,:,j);
            end
        end
 
       %% Generating WF
        WF1 = mean(Snoisy,3);
        WF1 = WF1(end/2+1-Np1/2:end/2+Np1/2,end/2+1-Np2/2:end/2+Np2/2);
        WF = imresize(WF1,[2*Np1,2*Np2]);
                
       %% Generate approximate OTF/PSF
        [OTFo, PSFe, Kotf, param] = Generate_OTF(Snoisy,param, NPixel);
        
       %% Preprocessing and parameter estimation
        % RL deconvlucy
        psf = fftshift(ifft2(ifftshift(OTFo)));
        [~,~,d3]=size(Snoisy);
        Snoisy = deconvlucy(Snoisy,psf,5);
        
        % parameter estimation
        s_class = single_orientation_class(OTFo,Kotf,Snoisy,PSFe);
        
        % frequency vector estimation
        if parainheri==0
        kAo_L = zeros(3,2); kAmean_L = zeros(3,2); K_ini_L= zeros(1,3);
        kAo_H = zeros(3,2); kAmean_H = zeros(3,2); K_ini_H= zeros(1,3);
        flag = 1;
        s_class = KapproxEstimationF(s_class,flag); % kAo contains three values along three directions
        s_class.kAmean = kmeanF(s_class);  % kAmean contains three values
        kAo_L = s_class.kAo;
        kAmean_L = s_class.kAmean; 
        TE = kAmean_L;
        K_ini_L = sqrt(TE(:,1).^2+TE(:,2).^2);
        
        if max(K_ini_L)-min(K_ini_L)>1.5 
           flag = 0;         
           s_class = KapproxEstimationF(s_class,flag);  
           s_class.kAmean = kmeanF(s_class);
           kAo_H = s_class.kAo;
           kAmean_H = s_class.kAmean; 
           TE1 = kAmean_H;
           K_ini_H = sqrt(TE1(:,1).^2+TE1(:,2).^2);
           
           if max(K_ini_H)-min(K_ini_H)>1.5
                T_L = max(K_ini_L)-min(K_ini_L);
                T_H = max(K_ini_H)-min(K_ini_H);
                fprintf('The difference between K values is too great! The result may be not good!')
                if T_L > T_H
                    s_class.kAo = kAo_H;
                    s_class.kAmean = kAmean_H;
                else
                    s_class.kAo = kAo_L;
                    s_class.kAmean = kAmean_L;
                end
           end
        end
 
        else
            s_class.kAo=s0_class{1}.kAo;
            s_class.kAmean=s0_class{1}.kAmean;
            [fDo_temp,fDp_temp,~,fDp_temp1,~,fDp_temp2,~] = PCMseparate_tempF(s_class);
            s_class.sepf(:,:,1) = fDo_temp;
            s_class.sepf(:,:,2) = fDp_temp;
            s_class.sepf(:,:,3) = fDp_temp1;
            s_class.sepf(:,:,4) = fDp_temp2;    
        end
        % phase estimation
        s_class.phaseA = phasesF(s_class);
        % modulation estimation
        [s_class.fDo,s_class.fDp,s_class.fDm,s_class.fDp1,s_class.fDm1,s_class.fDp2,s_class.fDm2] = PCMseparateF(s_class);  
        s_class.modFac = ModulationFactorF(s_class);
 
       %% Wiener filtering
        % object power parameters
        kA1 = s_class.kAmean(1,:);
        kB1 = s_class.kAmean(2,:);
        kC1 = s_class.kAmean(3,:);
        AA = [sqrt(kA1*kA1'),sqrt(kB1*kB1'),sqrt(kC1*kC1')];
        coe = max(AA)/s_class.Kotf;
        s_class.OBJpara = OBJpowerPara(s_class,OTFo,coe);
        % Wiener Filtering central frequency component 
        co = 1;
         s_class.fDof=s_class.fDo;
        % Wiener Filtering the noisy off-center frequency components    
        [s_class.FDpf,s_class.FDmf,s_class.NpDp,s_class.NpDm] = PCMfilteringF(s_class,co); 
       
       %% Reconstruction and Spectrum optimization
        param.attStrength = str2double(get(handles.att,'String'));
        param.a = 2.5;
        param.w1 = str2double(get(handles.Wienerpara,'String'));
        param.w2 = 0.1;
        
        param.w3 = param.w2;
        param.OtfProvider = Airy_giveOtf(param,param.a);
 
        kA = s_class.kAmean(1,:);
        kB = s_class.kAmean(2,:);
        kC = s_class.kAmean(3,:);
 
        param.Dir(1).px = kA(2); param.Dir(1).py = kA(1); 
        param.Dir(2).px = kB(2); param.Dir(2).py = kB(1); 
        param.Dir(3).px = kC(2); param.Dir(3).py = kC(1); 
        [w,h] = size(param.OtfProvider.otf);
        fftDirectlyCombined = zeros(2*h,2*w);
        shifted=zeros(2*h,2*w,7);
        shifted(:,:,1) = placeFreq(s_class.fDof);
        shifted(:,:,1) = applyOtf(shifted(:,:,1),param.OtfProvider,1,0,0,1,0);
 
        for I = 1:3
            shifted(:,:,2) = applyOtf(s_class.FDpf(:,:,I),param.OtfProvider,2,1*param.Dir(I).px,1*param.Dir(I).py,1,0);
            shifted(:,:,3) = applyOtf(s_class.FDmf(:,:,I),param.OtfProvider,2,-param.Dir(I).px,-param.Dir(I).py,1,0);
            fftDirectlyCombined = fftDirectlyCombined +sum(shifted,3);  
        end
 
        for i = 1:3
            param.Dir(i).phaOff = s_class.phaseA(i); 
            param.Dir(i).modul = s_class.modFac(i);
        end
        K0 = [sqrt(kA*kA') sqrt(kB*kB') sqrt(kC*kC')];
        K = max(ceil(K0));
        cutoff = floor(1*K)/param.sampleLateral+1.0;
        otfTS = zeros(2*h,2*w);
        otfTS = writeApoVector(otfTS,param.OtfProvider,cutoff);      % Ideal OTF
        MaskF = zeros(2*h,2*w);
        MaskF(otfTS~=0)=1;
        
        wFilter0=WienerFilterWiener_2D(param);
        Wk0=otfTS./(wFilter0.wDenom+param.w2^2);
        Wk03=otfTS./(wFilter0.wDenom+param.w3^2);

        directSIM=real(ifft2(fftshift((fftDirectlyCombined))));
        
        Wiener=real(ifft2(fftshift((fftDirectlyCombined.*Wk0.*MaskF))));
        SIM=Wiener;
        Wiener2=real(ifft2(fftshift((fftDirectlyCombined.*Wk03.*MaskF))));
        imwrite(uint16(Wiener/max(max(Wiener))*2^16),[param.pathname,'SIM',num2str(k),'WienerSIM.tif']);

% HiFi Spectrum optimization
HiFi=1;
if HiFi==1
        % Apply W1
        wFilter1=WienerFilterW1_2D(param);
        Wk1=otfTS./(wFilter1.wDenom+param.w1^2);
        fftInitialTS=fftDirectlyCombined.*Wk1;
        %figure; pcolor(Wk1); shading interp
        W1SIM=real(ifft2(fftshift(fftInitialTS)));
        %imwrite(uint16((2^16-1)*((W1SIM-min(min(W1SIM)))/(max(max(W1SIM))-min(min(W1SIM))))),[param.pathname,'SIM',num2str(k),'W1SIM.tif']);
        % Apply W2   
        ApoFWHM=0.5*(cutoff-1);
        ApoFWHM=min(0.5,round(ApoFWHM*100)/100);
        apo= apodize_gauss([2*h,2*w], struct('rad',ApoFWHM));  
        wFilter2=WienerFilterW2_2D(param);
        Wk2=apo./(wFilter2.wDenom+param.w2^2);     
        %figure; pcolor(Wk2); shading interp
        FW2 = fftInitialTS.*Wk2.*MaskF;
        W1W2SIM=real(ifft2(fftshift(FW2)));
        W1W2SIM=W1W2SIM/max(max(W1W2SIM))*(2^16-1);
        %imwrite(uint16(W1W2SIM),[param.pathname,'SIM',num2str(k),'W1W2SIM.tif']);
        SIM=W1W2SIM;
        SIM(SIM<0)=0;
%         Size1=2*param.Size1;
%         Size2=2*param.Size2;
%         TS=zeros(Size1,Size2);
%         Temp=Wiener;
%         Temp(Temp<0)=0;
%         Temp=255*Temp/max(max(Temp));
%         TS(1:Size1,1:Size2)=Temp(1:Size1,1:Size2);
%         TS=importImages2(TS);
%         SIM = TS(end/2+1-Np1:end/2+Np1,end/2+1-Np2:end/2+Np2);
end
if HexF==1
 centerbank=zeros(6,2);
 centerbank(1:3,:)=[kA1;kB1;kC1];
 centerbank(4:6,:)=-[kA1;kB1;kC1];
 FSIM_s=fftshift(fft2(SIM));
  sigma=20;
  [d1,d2]=size(FSIM_s);
  G_check=ones(d1,d2);
 x=1:d1;
 y=1:d2;
 x=x-fix(d1/2);
 y=y-fix(d2/2);
 [X,Y]=meshgrid(x,y);
 R=sqrt(X.^2+Y.^2);
 R=fix(R);
 I_curve=FRavg(FSIM_s);
 HexF=get(handles.HexF,'value');  % whethre to overlap all time-lapse 
        % Correct the SIM kink
for hex_dir=1:6    
     x_s=centerbank(hex_dir,2);
     y_s=centerbank(hex_dir,1);
     G=1/(2*pi*sigma^2)*exp(-(X-x_s).^2/(2*sigma^2)-(Y-y_s).^2/(2*sigma^2));
     G=G/max(max(G));
     G=1-G;
     G=G';
     mask1=G;
     mask1(mask1>0.95)=0;
     mask1(mask1>0)=1;
     Range=abs(FSIM_s).*mask1;
     [m,m1]=max(max(Range));
     [m,m2]=max(Range(:,m1));
     x_s=m1-fix(d1/2);
     y_s=m2-fix(d2/2);
     G=1/(2*pi*sigma^2)*exp(-(X-x_s).^2/(2*sigma^2)-(Y-y_s).^2/(2*sigma^2));
     G=G/max(max(G));
     G=1-G;
     Raw_int=max(max(abs(FSIM_s).*mask1));
     Sur_int=I_curve(fix(sqrt((y_s).^2+(x_s).^2)));
%      Sur_int=max(max(abs(FSIM_s).*mask2));
     c=Raw_int/Sur_int;
%      c=c+0.6;
     if c>1
     Gmin=0.6/(c-1);
     G=G+Gmin;
     G=G/max(max(G));
     G_check=G_check.*G;
     FSIM_s= FSIM_s.*G;
     %figure; pcolor(abs(G)); shading interp
     end
end
 %figure; pcolor(abs(G_check)); shading interp
 Wiener=real(ifft2(fftshift((FSIM_s))));
 Wiener2=Wiener;
 SIM=Wiener;
 end
 
    
% Inversion-based optimization
Deconvop=get(handles.Deconv,'value');  % whethre to overlap all time-lapse 
if Deconvop==1         
        SIM=MRA_deconv(SIM,param.lambda,param.pixelsize*10^3,param.NA,MRAreg,Mode);
        w1=size(SIM,2);
        w2=size(SIM,1);
        w1=size(SIM,2);
        w2=size(SIM,1);
        [X,Y]=meshgrid(linspace(1,w1,w1),linspace(1,w2,w2));
        X=X-fix(w1/2);
        Y=Y-fix(w2/2);
        R=sqrt(X.^2+Y.^2);
        Kf=0.5*sqrt((w1/2).^2+(w2/2).^2);
        GF=exp(-R.^2/(2*Kf^2));
        TempF=fftshift(fft2(SIM)).*GF;
        SIM=real(ifft2(fftshift(TempF)));        
end

        if k==1
           % if parainheri==0

            kA1 = s_class.kAmean(1,:);
            kB1 = s_class.kAmean(2,:);
            kC1 = s_class.kAmean(3,:);
                param.K0 = (sqrt(kA1*kA1')+sqrt(kB1*kB1')+sqrt(kC1*kC1'))/3;
                param.modul = (s_class.modFac(1)+s_class.modFac(2)+s_class.modFac(3))/3;
                param.coe = 1+param.K0/s_class.Kotf;

            fileID = fopen([param.pathname,'param.txt'], 'a');   
            fprintf(fileID, 'Frequency vector\n');
            fprintf(fileID, '%s\n', kA1);  
            fprintf(fileID, '%s\n', kB1);
            fprintf(fileID, '%s\n', kC1);
            fprintf(fileID, 'Modulation factor\n');
            fprintf(fileID, '%s\n', s_class.modFac);
            fprintf(fileID, 'Initial phase\n');
            fprintf(fileID, '%s\n', s_class.phaseA);
            fclose(fileID);
             K1_dis=num2str(sqrt(kA1*kA1')); %K1_dis=K1_dis(1:6);
             K2_dis=num2str(sqrt(kB1*kB1')); %K2_dis=K2_dis(1:6);
             K3_dis=num2str(sqrt(kC1*kC1')); %K3_dis=K3_dis(1:6);
             set(handles.k1,'String',K1_dis);
             set(handles.k2,'String',K2_dis);
             set(handles.k3,'String',K3_dis);
             if s_class.modFac(1)>1
                s_class.modFac(1)=1;
            end
             if s_class.modFac(2)>1
                s_class.modFac(2)=1;
             end
             if s_class.modFac(3)>1
                s_class.modFac(3)=1;
             end 

             if parainheri==0
             angle1=num2str(rad2deg(atan(kA1(2)/kA1(1)))); %angle1=angle1(1:6);
             angle2=num2str(rad2deg(atan(kB1(2)/kB1(1)))); %angle2=angle2(1:6);
             angle3=num2str(rad2deg(atan(kC1(2)/kC1(1)))); %angle3=angle3(1:6);
             set(handles.A1,'String',angle1);
             set(handles.A2,'String',angle2);
             set(handles.A3,'String',angle3);
             end
            %end
        end
        s=sprintf('Recon3I-SIM reconstruction processing:%d',ceil(k/nz*100));
        waitbar(k/nz,ha,[s '%']);
        %imwrite(uint16(SIM),[param.pathname,'SIM',num2str(k),'.tif']);
        clear s_class
        WF = (WF-min(WF(:)))/(max(WF(:))-min(WF(:)));
        SIM = (SIM-min(SIM(:)))/(max(SIM(:))-min(SIM(:)));

        avg1 = mean(WF(:));
        avg2= mean(SIM(:));
%         if avg2>avg1&&Deconvop~=1
%         SIM = SIM + avg1-avg2 + 0.05185 *log(avg1+0.2);
%         SIM(SIM<0)=0;
%         end 
        SIM = (SIM-min(SIM(:)))/(max(SIM(:))-min(SIM(:)));
        
           if k==1
             axes(handles.axes2);
            imshow(SIM,[])                
           end
        imwrite(uint16((2^16-1)*WF), [DataDir,'new3I_SIM','\','WF_',num2str(filenum),'.tif'], 'WriteMode','append')
        imwrite(uint16((2^16-1)*SIM), [DataDir,'new3I_SIM','\','SIM_',num2str(filenum),'.tif'], 'WriteMode','append')
        info = imfinfo(fullfile([DataDir,'new3I_SIM','\','SIM_',num2str(filenum),'.tif']));
        frames = numel(info);      
        if info(1).Height*info(1).Width*frames>1024*1024*1500
           filenum=filenum+1; 
        end
     end
          close(ha)
          toc
   
  
 
 
 
 
 
 
 
 
function att_Callback(hObject, eventdata, handles)
% hObject    handle to att (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of att as text
%        str2double(get(hObject,'String')) returns contents of att as a double
 
 
% --- Executes during object creation, after setting all properties.
function att_CreateFcn(hObject, eventdata, handles)
% hObject    handle to att (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function a_Callback(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of a as text
%        str2double(get(hObject,'String')) returns contents of a as a double
 
 
% --- Executes during object creation, after setting all properties.
function a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
% --- Executes on button press in Parainher.
function Parainher_Callback(hObject, eventdata, handles)
% hObject    handle to Parainher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hint: get(hObject,'Value') returns toggle state of Parainher
 
 
% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hint: get(hObject,'Value') returns toggle state of radiobutton2
 
 
% --- Executes on button press in Rolling.
function Rolling_Callback(hObject, eventdata, handles)
% hObject    handle to Rolling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hint: get(hObject,'Value') returns toggle state of Rolling
 
 
 
function wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of wavelength as text
%        str2double(get(hObject,'String')) returns contents of wavelength as a double
 
 
% --- Executes during object creation, after setting all properties.
function wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function pixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of pixelsize as text
%        str2double(get(hObject,'String')) returns contents of pixelsize as a double
 
 
% --- Executes during object creation, after setting all properties.
function pixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function NA_Callback(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of NA as text
%        str2double(get(hObject,'String')) returns contents of NA as a double
 
 
% --- Executes during object creation, after setting all properties.
function NA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function k1_Callback(hObject, eventdata, handles)
% hObject    handle to k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of k1 as text
%        str2double(get(hObject,'String')) returns contents of k1 as a double
 
 
% --- Executes during object creation, after setting all properties.
function k1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function k2_Callback(hObject, eventdata, handles)
% hObject    handle to k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of k2 as text
%        str2double(get(hObject,'String')) returns contents of k2 as a double
 
 
% --- Executes during object creation, after setting all properties.
function k2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function k3_Callback(hObject, eventdata, handles)
% hObject    handle to k3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of k3 as text
%        str2double(get(hObject,'String')) returns contents of k3 as a double
 
 
% --- Executes during object creation, after setting all properties.
function k3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function A1_Callback(hObject, eventdata, handles)
% hObject    handle to A1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of A1 as text
%        str2double(get(hObject,'String')) returns contents of A1 as a double
 
 
% --- Executes during object creation, after setting all properties.
function A1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function A2_Callback(hObject, eventdata, handles)
% hObject    handle to A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of A2 as text
%        str2double(get(hObject,'String')) returns contents of A2 as a double
 
 
% --- Executes during object creation, after setting all properties.
function A2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function A3_Callback(hObject, eventdata, handles)
% hObject    handle to A3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of A3 as text
%        str2double(get(hObject,'String')) returns contents of A3 as a double
 
 
% --- Executes during object creation, after setting all properties.
function A3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function M1_Callback(hObject, eventdata, handles)
% hObject    handle to M1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of M1 as text
%        str2double(get(hObject,'String')) returns contents of M1 as a double
 
 
% --- Executes during object creation, after setting all properties.
function M1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function M2_Callback(hObject, eventdata, handles)
% hObject    handle to M2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of M2 as text
%        str2double(get(hObject,'String')) returns contents of M2 as a double
 
 
% --- Executes during object creation, after setting all properties.
function M2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
 
function M3_Callback(hObject, eventdata, handles)
% hObject    handle to M3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of M3 as text
%        str2double(get(hObject,'String')) returns contents of M3 as a double
 
 
% --- Executes during object creation, after setting all properties.
function M3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
% --- Executes on selection change in pop1.
function pop1_Callback(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: contents = cellstr(get(hObject,'String')) returns pop1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop1
 
 
 
% --- Executes during object creation, after setting all properties.
function pop1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
% --- Executes on selection change in pop1.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: contents = cellstr(get(hObject,'String')) returns pop1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop1
 
 
% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
 
% --- Executes on button press in HexF.
function HexF_Callback(hObject, eventdata, handles)
% hObject    handle to HexF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hint: get(hObject,'Value') returns toggle state of HexF
 
 
 
function Wienerpara_Callback(hObject, eventdata, handles)
% hObject    handle to Wienerpara (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: get(hObject,'String') returns contents of Wienerpara as text
%        str2double(get(hObject,'String')) returns contents of Wienerpara as a double
 
 
% --- Executes during object creation, after setting all properties.
function Wienerpara_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wienerpara (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Deconv.
function Deconv_Callback(hObject, eventdata, handles)
% hObject    handle to Deconv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Deconv
