%% TIME RESOLVED ANISOTROPY METHOD (TRAM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%LAST UPDATED 160222 BY Thomas S van Zanten%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/thomas/Documents/MATLAB/added_codes/bfmatlab/'))
addpath('/Users/thomas/Documents/MATLAB/GitHub/TimeResolvedAnisotropy/')
var.pathname=uigetdir;
var.pathname=[var.pathname '/'];
cd(var.pathname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%DEFININING ALL THE VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var.scrsz = get(0,'ScreenSize');
var.opts=optimset('LargeScale', 'off','Display','off');
%life-/rotational-time indicators and proportional values of irf/data 
var.lt=1; var.rt=1; var.Nf1=1; var.Nf2=1; var.shift_ALL=0;
%lifetime variables
var.a=0.03; var.b=0.2; var.c=0.2; var.tau_lt1=0.5; var.tau_lt2=2; var.tau_lt3=6;
%rotational variables
var.r0=0.4; var.d=0.4; var.e=0.1; var.tau_r1=0.2; var.tau_r2=14; var.tau_r3=100;
%Gfactor and initial channel shift
var.Gf=1.04; var.shift_iPA=295; var.shift_iPE=260;
var.L=100; var.rep=12.5;%%%Interpolation value & RepRate in ns
%%%%%%%%OPEN DATA AND SEPARATE PA AND PE FILES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=bfopen([]);
data_PA=cell2mat(data{2,1}(:,1)); data_PA=double(data_PA);
data_PE=cell2mat(data{1,1}(:,1));  data_PE=double(data_PE);
%%%%%%%%OPEN IRF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irf=bfopen([]);
irf_PA=cell2mat(irf{2,1}(:,1)); irf_PA=double(irf_PA);
irf_PE=cell2mat(irf{1,1}(:,1)); irf_PE=double(irf_PE);
%%%%%%%%OPEN BACKGROUND%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BG=bfopen([]);
%BG_PA=cell2mat(BG{2,1}(:,1)); BG_PA=double(BG_PA);
%BG_PE=cell2mat(BG{1,1}(:,1)); BG_PE=double(BG_PE);
%data_PA=data_PA-BG_PA;
%data_PE=data_PE-BG_PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:(var.rep/length(irf_PA)):var.rep-(var.rep/length(irf_PA)); t=t';%time axis in ns
%create temporary file to be able to update continuously
Temp_PA=data_PA;
Temp_PE=data_PE;
%normalize such that irf goes down to 0
irf_PA=irf_PA-ceil(mean(irf_PA(find(t>5 & t<8))));
irf_PE=irf_PE-ceil(mean(irf_PE(find(t>5 & t<8))));
irf_PA(irf_PA<0)=0;%everything should be >0
irf_PE(irf_PE<0)=0;
%normalize the irf with the data file and create a temporary irf file
irf_PA=irf_PA./max(irf_PA).*max(data_PA); Temp_irf_PA=irf_PA;
irf_PE=irf_PE./max(irf_PE).*max(data_PA); Temp_irf_PE=irf_PE;

%calculation of the total intensity and total irf decay
data_total=data_PA+2.*var.Gf.*data_PE; irf_total=irf_PA+2.*var.Gf.*irf_PE;
%define the empty vectors that are to be filled during the fitting process
G0_total(1:length(t),1)=0; G0_PA(1:length(t),1)=0; G0_PE(1:length(t),1)=0;

%future fitting analysis is only performed on the regions containing real data in both channels
var.start=max(min(find(data_PE>0)),min(find(data_PA>0)));
var.end=min(max(find(data_PE(var.start:end)>0)+var.start-1),...
    max(find(data_PA(var.start:end)>0)+var.start-1))-50;
res_T=zeros(size((var.start:var.end)'));
res_PA=zeros(size((var.start:var.end)'));res_PE=zeros(size((var.start:var.end)'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE DISPLAYING UPDATED DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F1=figure('Name','Time Resolved Decays PA and PE','NumberTitle','off');
set(F1,'OuterPosition',[1 var.scrsz(4)/2 var.scrsz(3) var.scrsz(4)]);
subplot(3,3,[1,4]);
semilogy(t(var.start:var.end),data_total(var.start:var.end),'ko','XDataSource',...
    't(var.start:var.end)','YDataSource','data_total(var.start:var.end)')
hold on
semilogy(t(var.start:var.end),irf_total(var.start:var.end),'r-','XDataSource',...
    't(var.start:var.end)','YDataSource','irf_total(var.start:var.end)')
hold on
semilogy(t(var.start:var.end),G0_total(var.start:var.end),'g-','LineWidth',2,...
    'XDataSource','t(var.start:var.end)','YDataSource','G0_total(var.start:var.end)') 
title('Intensity decay')
axis([0.5 12 mean(data_total(var.start:var.start+10))/2 2*max(data_total)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,[2,5]);
semilogy(t(var.start:var.end),data_PA(var.start:var.end),'bo','XDataSource',...
    't(var.start:var.end)','YDataSource','data_PA(var.start:var.end)')
hold on
semilogy(t(var.start:var.end),irf_PA(var.start:var.end),'r-','XDataSource',...
    't(var.start:var.end)','YDataSource','irf_PA(var.start:var.end)')
hold on
semilogy(t(var.start:var.end),G0_PA(var.start:var.end),'k-','LineWidth',2,...
    'XDataSource','t(var.start:var.end)','YDataSource','G0_PA(var.start:var.end)')
title('PARALLEL (PA or VV)')
axis([0.5 12 mean(data_PA(var.start:var.start+10))/2 2*max(data_PA)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,[3,6]);
semilogy(t(var.start:var.end),data_PE(var.start:var.end),'go','XDataSource',...
    't(var.start:var.end)','YDataSource','data_PE(var.start:var.end)')
hold on
semilogy(t(var.start:var.end),irf_PE(var.start:var.end),'r-','XDataSource',...
    't(var.start:var.end)','YDataSource','irf_PE(var.start:var.end)')
hold on
semilogy(t(var.start:var.end),G0_PE(var.start:var.end),'k-','LineWidth',2,...
    'XDataSource','t(var.start:var.end)','YDataSource','G0_PE(var.start:var.end)')
title('PERPENDICULAR (PE or VH)')
axis([0.5 12 mean(data_PE(var.start:var.start+10))/2 2*max(data_PE)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%RESIDUALS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%res_T=(data_total(var.start:var.end)-G0_total(var.start:var.end))/data_total(var.start:var.end);
%res_PA=(data_PA(var.start:var.end)-G0_PA(var.start:var.end))/data_PA(var.start:var.end);
%res_PE=(data_PE(var.start:var.end)-G0_PE(var.start:var.end))/data_PE(var.start:var.end);
subplot(3,3,[7,8,9]);
plot(t(var.start:var.end),res_T,'k-','LineWidth',2,'XDataSource','t(var.start:var.end)',...
    'YDataSource','res_T')
hold on
plot(t(var.start:var.end),res_PA,'b-','LineWidth',2,'XDataSource','t(var.start:var.end)',...
    'YDataSource','res_PA')
hold on
plot(t(var.start:var.end),res_PE,'g-','LineWidth',2,'XDataSource','t(var.start:var.end)',...
    'YDataSource','res_PE')
title('RESIDUALS')
axis([0.5 12 -0.3 0.3])
linkdata on%ALLOWING FOR CONTINUOUS UPDATED CURVES WHILE CHANGING VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPDATING AND USING THE REQUIRED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var.r_inf=(sum(sum(data_PA))-sum(sum(data_PE)))/sum(sum(data_total));

uicontrol('Style','text','Position',[20 700 60 20],'String','No. tau_l');
han.lt=uicontrol('Style','edit','Position',[20 680 60 20],'String',num2str(var.lt));
uicontrol('Style','text','Position',[100 700 60 20],'String','No. tau_r');
han.rt=uicontrol('Style','edit','Position',[100 680 60 20],'String',num2str(var.rt));
uicontrol('Style','text','Position',[20 650 60 20],'String','Gfactor');
han.Gf=uicontrol('Style','edit','Position',[20 630 60 20],'String',num2str(var.Gf));
uicontrol('Style','text','Position',[100 650 60 20],'String','NormF1');
han.Nf1=uicontrol('Style','edit','Position',[100 630 60 20],'String',num2str(var.Nf1));
uicontrol('Style','text','Position',[100 600 60 20],'String','NormF2');
han.Nf2=uicontrol('Style','edit','Position',[100 580 60 20],'String',num2str(var.Nf2));
uicontrol('Style','text','Position',[20 600 60 20],'String','Ratio');
han.RNf=uicontrol('Style','edit','Position',[20 580 60 20],'String',num2str(var.Nf2./var.Nf1));

uicontrol('Style','text','Position',[30 540 120 20],'String','INTERPOLATION = 100');

uicontrol('Style','text','Position',[60 500 60 20],'String','shift PA');
han.shift_iPA=uicontrol('Style','edit','Position',[60 480 60 20],'String',num2str(var.shift_iPA));
uicontrol('Style','text','Position',[60 450 60 20],'String','shift PE');
han.shift_iPE=uicontrol('Style','edit','Position',[60 430 60 20],'String',num2str(var.shift_iPE));
uicontrol('Style','text','Position',[60 400 60 20],'String','shift ALL');
han.shift_ALL=uicontrol('Style','edit','Position',[60 380 60 20],'String',num2str(var.shift_ALL));

uicontrol('Style','text','Position',[var.scrsz(3)-120 700 60 20],'String','a');
han.a=uicontrol('Style','edit','Position',[var.scrsz(3)-120 680 60 20],'String',num2str(var.a));
uicontrol('Style','text','Position',[var.scrsz(3)-120 650 60 20],'String','b');
han.b=uicontrol('Style','edit','Position',[var.scrsz(3)-120 630 60 20],'String',num2str(var.b));
uicontrol('Style','text','Position',[var.scrsz(3)-120 600 60 20],'String','c');
han.c=uicontrol('Style','edit','Position',[var.scrsz(3)-120 580 60 20],'String',num2str(var.c));

uicontrol('Style','text','Position',[var.scrsz(3)-120 550 60 20],'String','tau_lt1');
han.tault1=uicontrol('Style','edit','Position',[var.scrsz(3)-120 530 60 20],'String',num2str(var.tau_lt1));
uicontrol('Style','text','Position',[var.scrsz(3)-120 500 60 20],'String','tau_lt2');
han.tault2=uicontrol('Style','edit','Position',[var.scrsz(3)-120 480 60 20],'String',num2str(var.tau_lt2));
uicontrol('Style','text','Position',[var.scrsz(3)-120 450 60 20],'String','tau_lt3');
han.tault3=uicontrol('Style','edit','Position',[var.scrsz(3)-120 430 60 20],'String',num2str(var.tau_lt3));

uicontrol('Style','text','Position',[var.scrsz(3)-120 400 60 20],'String','r_inf');
han.rinf=uicontrol('Style','edit','Position',[var.scrsz(3)-120 380 60 20],'String',num2str(var.r_inf));
uicontrol('Style','text','Position',[var.scrsz(3)-120 350 60 20],'String','r0');
han.r0=uicontrol('Style','edit','Position',[var.scrsz(3)-120 330 60 20],'String',num2str(var.r0));

uicontrol('Style','text','Position',[var.scrsz(3)-120 300 60 20],'String','d');
han.d=uicontrol('Style','edit','Position',[var.scrsz(3)-120 280 60 20],'String',num2str(var.d));
uicontrol('Style','text','Position',[var.scrsz(3)-120 250 60 20],'String','e');
han.e=uicontrol('Style','edit','Position',[var.scrsz(3)-120 230 60 20],'String',num2str(var.e));

uicontrol('Style','text','Position',[var.scrsz(3)-120 200 60 20],'String','tau1');
han.taur1=uicontrol('Style','edit','Position',[var.scrsz(3)-120 180 60 20],'String',num2str(var.tau_r1));
uicontrol('Style','text','Position',[var.scrsz(3)-120 150 60 20],'String','tau2');
han.taur2=uicontrol('Style','edit','Position',[var.scrsz(3)-120 130 60 20],'String',num2str(var.tau_r2));
uicontrol('Style','text','Position',[var.scrsz(3)-120 100 60 20],'String','tau3');
han.taur3=uicontrol('Style','edit','Position',[var.scrsz(3)-120 80 60 20],'String',num2str(var.tau_r3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FITTING, UPDATING, ACCEPTING AND SAVING INTENSITY DECAYS UPON BUTTON CLICK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
btn1 = uicontrol('Style', 'pushbutton', 'String', 'Refresh',...
'Position', [40 120 100 40],'Callback',...
'[data_total,data_PA,data_PE,irf_total,irf_PA,irf_PE,G0_total,G0_PA,G0_PE,han,var] = shiftNfit(t,Temp_PA,Temp_PE,Temp_irf_PA,Temp_irf_PE,han,var);res_T=(data_total(var.start:var.end)-G0_total(var.start:var.end))./data_total(var.start:var.end);res_PA=(data_PA(var.start:var.end)-G0_PA(var.start:var.end))./data_PA(var.start:var.end);res_PE=(data_PE(var.start:var.end)-G0_PE(var.start:var.end))./data_PE(var.start:var.end);');

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Accept',...
'Position', [40 80 100 40],'Callback',...
'[anis,anis_G0,tt,filename,var] = accept(t,data_total,irf_total,G0_total,data_PA,irf_PA,G0_PA,data_PE,irf_PE,G0_PE,var); clear btn1 btn2 F1 han Temp_irf_PA Temp_irf_PE Temp_PA Temp_PE; save([var.pathname filename])');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%