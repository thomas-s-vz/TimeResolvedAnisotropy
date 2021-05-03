%% The shiftNfit function performs the fitting on all the data after the IRF is shifted with respect to the data
%Last adaption Thomas S van Zanten 160222
function [data_total,data_PA,data_PE,irf_total,irf_PA,irf_PE,G0_total,G0_PA,G0_PE,han,var]=shiftNfit(t,Temp_PA,Temp_PE,Temp_irf_PA,Temp_irf_PE,han,var)
%% GET THE NUMBERS FROM THE FIGURE (TEMPORARY PRESENT IN HAN) %%%%%%%%%%%%%
var.shift_iPA=str2num(get(han.shift_iPA,'String')); var.shift_iPE=str2num(get(han.shift_iPE,'String')); var.shift_ALL=str2num(get(han.shift_ALL,'String')); 

var.Gf=str2num(get(han.Gf,'String')); var.lt=str2num(get(han.lt,'String')); var.rt=str2num(get(han.rt,'String'));
var.a=str2num(get(han.a,'String')); var.b=str2num(get(han.b,'String')); var.c=str2num(get(han.c,'String')); var.tau_lt1=str2num(get(han.tault1,'String')); var.tau_lt2=str2num(get(han.tault2,'String')); var.tau_lt3=str2num(get(han.tault3,'String'));
var.r0=str2num(get(han.r0,'String')); var.r_inf=str2num(get(han.rinf,'String')); var.Nf1=str2num(get(han.Nf1,'String')); var.Nf2=str2num(get(han.Nf2,'String')); 
var.d=str2num(get(han.d,'String')); var.e=str2num(get(han.e,'String')); var.tau_r1=str2num(get(han.taur1,'String')); var.tau_r2=str2num(get(han.taur2,'String')); var.tau_r3=str2num(get(han.taur3,'String'));
%% ASSIGN AND SHIFT THE IRF AND DATA FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irf_PA=Temp_irf_PA; irf_PE=Temp_irf_PE; data_PA=Temp_PA; %data_PE=Temp_PE; GETS DEFINED AT line28

ts = interp1(t,1:1/var.L:length(t))';%define new time scale based on interpolation value L
    IRF_PA = interp1(t,irf_PA,ts,'spline');%interpolate the irf to fit the new timescale
    IRF_PE = interp1(t,irf_PE,ts,'spline');
    IRF_PA = circshift(IRF_PA,var.shift_iPA);%shift the interpolated irf with the defined number with respect to the data file
    IRF_PE = circshift(IRF_PE,var.shift_iPE);
  
    for i=1:2000
        shift=-1000+i;
        IRF = circshift(IRF_PE,shift);
        shift_space(i)= sum((IRF_PA-IRF).^2);
    end
var.shift_PE=find(shift_space==min(shift_space))-1000;%find the different between the PA and PE starting points using maximum overlay values        
%var.shift_PE=round(find(IRF_PA==max(IRF_PA))-find(IRF_PE==max(IRF_PE)));%find the different between the PA and PE starting points using peak values

IRF_PE = circshift(IRF_PE,var.shift_PE+var.shift_ALL);%shift the irf of PE again to allign to the PA (+the additional forced shift)
    
    irf_PA = IRF_PA(1:var.L:end);%reduce the interpolated and shifted irf's to their original size
 	irf_PE = IRF_PE(1:var.L:end);
    irf_PA(irf_PA<0) = 0;%avoid the existence of values <0
    irf_PE(irf_PE<0) = 0;
    
    DATA_PE = interp1(t,Temp_PE,ts,'spline');%interpolate the PE data to fit the new timescale
    DATA_PE = circshift(DATA_PE,var.shift_PE+var.shift_ALL);%shift the data of PE to allign to the PA channel (+the additional forced shift)
    data_PE = DATA_PE(1:var.L:end);%reduce the interpolated and shifted PE data file to its original size
    data_PE(data_PE<0) = 0;%avoid the existence of values <0
    
data_total=data_PA+2.*var.Gf.*data_PE; irf_total=irf_PA+2.*var.Gf.*irf_PE;%calculation of the total intensity and total irf decay after all the shifting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA FITTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE LIFETIME COMPONENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if var.lt==1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_lt=[var.a var.tau_lt1];
    [P_lt]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,irf_total,0,0,data_total,0,0,0,var),par_lt,[],[],[],[],...
        [0 0],[1 10],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P_lt,var.lt,var.rt,t,irf_total,0,0,data_total,0,0,0,var);
        G0_total=G_T;

    %%%%AND SINGLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if var.rt==1        
s_par0=[var.Nf1 var.Nf2 var.r0 var.tau_r1]; 
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),s_par0,[],[],[],[],...
        [0 0 0 0],[10 10 0.6 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;

    %%%%AND DOUBLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        elseif var.rt==2          
d_par0=[var.Nf1 var.Nf2 var.r0 var.d var.tau_r1 var.tau_r2];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),d_par0,[],[],[],[],...
        [0 0 0 0 0 0],[10 10 0.6 1 5 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;
        
    %%%%AND TRIPLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        elseif var.rt==3          
d_par0=[var.Nf1 var.Nf2 var.r0 var.d var.tau_r1 var.tau_r2 var.e var.tau_r3];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),d_par0,[],[],[],[],...
        [0 0 0 0 0 0 0 0],[10 10 0.6 2 5 20 2 200],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;
        
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% DOUBLE LIFETIME COMPONENTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif var.lt==2

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_lt=[var.a var.b var.tau_lt1 var.tau_lt2];
    [P_lt]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,irf_total,0,0,data_total,0,0,0,var),par_lt,[],[],[],[],...
        [0 0 0 2],[1 1 5 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P_lt,var.lt,var.rt,t,irf_total,0,0,data_total,0,0,0,var);
        G0_total=G_T;
        
    %%%%AND SINGLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if var.rt==1
s_par0=[var.Nf1 var.Nf2 var.r0 var.tau_r1];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),s_par0,[],[],[],[],...
        [0 0 0 0],[10 10 0.6 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;

    %%%%AND DOUBLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        elseif var.rt==2
d_par0=[var.Nf1 var.Nf2 var.r0 var.d var.tau_r1 var.tau_r2];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),d_par0,[],[],[],[],...
        [0 0 0 0 0 0],[10 10 0.6 1 5 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;
        
    %%%%AND TRIPLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        elseif var.rt==3
d_par0=[var.Nf1 var.Nf2 var.r0 var.d var.tau_r1 var.tau_r2 var.e var.tau_r3];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),d_par0,[],[],[],[],...
        [0 0 0 0 0 0 0 0],[10 10 0.6 2 5 20 2 200],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;
        
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% TRIPLE LIFETIME COMPONENTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif var.lt==3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_lt=[var.a var.b var.c var.tau_lt1 var.tau_lt2 var.tau_lt3];
    [P_lt]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,irf_total,0,0,data_total,0,0,0,var),par_lt,[],[],[],[],...
        [0 0 0 0 0.5 2],[1 1 10 1 5 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P_lt,var.lt,var.rt,t,irf_total,0,0,data_total,0,0,0,var);
        G0_total=G_T;
        
    %%%%AND SINGLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if var.rt==1
s_par0=[var.Nf1 var.Nf2 var.r0 var.tau_r1];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),s_par0,[],[],[],[],...
        [0 0 0 0],[10 10 0.6 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;

    %%%%AND DOUBLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        elseif var.rt==2
d_par0=[var.Nf1 var.Nf2 var.r0 var.d var.tau_r1 var.tau_r2];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),d_par0,[],[],[],[],...
        [0 0 0 0 0 2],[10 10 0.6 1 5 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;
        
    %%%%AND TRIPLE ROTATIONAL COMPONENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        elseif var.rt==3
d_par0=[var.Nf1 var.Nf2 var.r0 var.d var.tau_r1 var.tau_r2 var.e var.tau_r3];
    [P]=fmincon(@(par) Gfit_min(par,var.lt,var.rt,t,0,irf_PA,irf_PE,0,data_PA,data_PE,P_lt,var),d_par0,[],[],[],[],...
        [0 0 0 0 0 0 0 0],[10 10 0.6 2 5 20 2 100],[],var.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,irf_PA,0,0,data_PA,0,P_lt,var);
        G0_PA=G_tPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [G_T,G_tPA,G_tPE]=Gfit_fit(P,var.lt,var.rt,t,0,0,irf_PE,0,0,data_PE,P_lt,var);
        G0_PE=G_tPE;
        
        end
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESETTING ALL THE VALUES IN THE MAIN FIGURE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var.r_inf=(sum(sum(data_PA))-var.Gf.*sum(sum(data_PE)))/sum(sum(data_total));
if var.lt==1
    set(han.a,'String',num2str(P_lt(1)))
    set(han.b,'String',0)     
    set(han.c,'String',0)
    set(han.tault1,'String',num2str(P_lt(2)))    
    set(han.tault2,'String',0)  
	set(han.tault3,'String',0)
var.a=P_lt(1); var.b=0; var.c=0; var.tau_lt1=P_lt(2); var.tau_lt2=0; var.tau_lt3=0;
elseif var.lt==2
    set(han.a,'String',num2str(P_lt(1)))
    set(han.b,'String',num2str(P_lt(2)))     
    set(han.c,'String',0)
    set(han.tault1,'String',num2str(P_lt(3)))    
    set(han.tault2,'String',num2str(P_lt(4)))  
	set(han.tault3,'String',0)
var.a=P_lt(1); var.b=P_lt(2); var.c=0; var.tau_lt1=P_lt(3); var.tau_lt2=P_lt(4); var.tau_lt3=0;
elseif var.lt==3
    set(han.a,'String',num2str(P_lt(1)))
    set(han.b,'String',num2str(P_lt(2)))     
    set(han.c,'String',num2str(P_lt(3)))
    set(han.tault1,'String',num2str(P_lt(4)))    
    set(han.tault2,'String',num2str(P_lt(5)))  
	set(han.tault3,'String',num2str(P_lt(6)))
var.a=P_lt(1); var.b=P_lt(2); var.c=P_lt(3); var.tau_lt1=P_lt(4); var.tau_lt2=P_lt(5); var.tau_lt3=P_lt(6);
end
    
if var.rt==1
    set(han.Nf1,'String',num2str(P(1)))
    set(han.Nf2,'String',num2str(P(2)))
    set(han.r0,'String',num2str(P(3)))
	set(han.d,'String',0)
    set(han.e,'String',0)
    set(han.taur1,'String',num2str(P(4)))
    set(han.taur2,'String',0)
    set(han.taur3,'String',0)    
    set(han.rinf,'String',num2str(var.r_inf))
	set(han.RNf,'String',num2str(P(2)/P(1)))
var.Nf1=P(1); var.Nf2=P(2); var.r0=P(3); var.d=0; var.e=0; var.tau_r1=P(4); var.tau_r2=0; var.tau_r3=0;
elseif var.rt==2
    set(han.Nf1,'String',num2str(P(1)))
    set(han.Nf2,'String',num2str(P(2)))
    set(han.r0,'String',num2str(P(3)))
	set(han.d,'String',num2str(P(4)))
    set(han.e,'String',0)    
    set(han.taur1,'String',num2str(P(5)))
    set(han.taur2,'String',num2str(P(6)))
    set(han.taur3,'String',0)    
    set(han.rinf,'String',num2str(var.r_inf))
	set(han.RNf,'String',num2str(P(2)/P(1)))
var.Nf1=P(1); var.Nf2=P(2); var.r0=P(3); var.d=P(4); var.e=0; var.tau_r1=P(5); var.tau_r2=P(6); var.tau_r3=0;
elseif var.rt==3
    set(han.Nf1,'String',num2str(P(1)))
    set(han.Nf2,'String',num2str(P(2)))
    set(han.r0,'String',num2str(P(3)))
	set(han.d,'String',num2str(P(4)))
    set(han.e,'String',num2str(P(7)))
    set(han.taur1,'String',num2str(P(5)))
    set(han.taur2,'String',num2str(P(6)))
    set(han.taur3,'String',num2str(P(8)))
    set(han.rinf,'String',num2str(var.r_inf))
	set(han.RNf,'String',num2str(P(2)/P(1)))
var.Nf1=P(1); var.Nf2=P(2); var.r0=P(3); var.d=P(4); var.e=P(7); var.tau_r1=P(5); var.tau_r2=P(6); var.tau_r3=P(8);
end

end