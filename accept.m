%% The accept function generates the respective fit figures and anisotropy figures to be saved as .png
%Last adaption Thomas S van Zanten 160222
function [anis,anis_G0,tt,filename,var]=accept(t,data_total,irf_total,G0_total,data_PA,irf_PA,G0_PA,data_PE,irf_PE,G0_PE,var)
%%%%%check lines 29 and 59-60!!!!!! for usage with 256 channels as well
close all
%% RECALCULATE THE PERCENTAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var.b=var.b/(var.b+(1-var.b)+var.c);
var.c=var.c/(var.b+(1-var.b)+var.c);
var.d=var.d/(var.d+(1-var.d)+var.e);
var.e=var.e/(var.d+(1-var.d)+var.e);
var.r0=var.r0/(var.d+(1-var.d)+var.e);
%% FILENAME CHOICE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B{1,1} = 'Test';
B{1,2} = 'FileName';
promptB = B(1,2); defB = B(1,1); AddOpts.Resize = 'on'; AddOpts.WindowStyle = 'normal';
AddOpts.Interpreter = 'tex'; ParametersB(1) = inputdlg(promptB,'Please give values', 1, defB, AddOpts);
filename_s=ParametersB{1}; filename=[filename_s '_ANA.mat'];
%% FINAL FIGURE PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
F2=figure('Name','TOTAL INTENSITY DECAY (LIFETIME)','NumberTitle','off');
set(F2,'OuterPosition',[1 var.scrsz(4)/2 var.scrsz(3)/2 var.scrsz(4)]);
subplot(3,2,[1,2,3,4]);
semilogy(t(var.start:var.end),data_total(var.start:var.end),'ko')
hold on
semilogy(t(var.start:var.end),irf_total(var.start:var.end),'r-')
hold on
semilogy(t(var.start:var.end),G0_total(var.start:var.end),'g-','LineWidth',2) 
title('INTENSITY DECAY')
axis([0.5 12 mean(data_total(var.start:var.start+10))/2 2*max(data_total)])
subplot(3,2,[5,6]);
plot(t(var.start:var.end),data_total(var.start:var.end)-G0_total(var.start:var.end),'k-','LineWidth',2)
title('RESIDUALS')
%axis([0.5 12 min(data_total(70:950)-G0_total(70:950)) max(data_total(70:950)-G0_total(70:950))])
set(findall(F2,'-property','FontSize'),'FontSize',18)
set(get(F2,'Children'),'LineWidth',2,'FontSize',20,'fontWeight','bold')
    print ('-r300', '-dpng', [var.pathname filename_s 'T.png'])

F3=figure('Name','Time Resolved Decays PA and PE','NumberTitle','off');
set(F3,'OuterPosition',[1 var.scrsz(4)/2 var.scrsz(3) var.scrsz(4)]);
subplot(3,2,[1,3]);
semilogy(t(var.start:var.end),data_PA(var.start:var.end),'bo')
hold on
semilogy(t(var.start:var.end),irf_PA(var.start:var.end),'r-')
hold on
semilogy(t(var.start:var.end),G0_PA(var.start:var.end),'k-','LineWidth',2)
title('PARALLEL (PA or VV)')
axis([0.5 12 mean(data_PA(var.start:var.start+10))/2 2*max(data_PA)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,[2,4]);
semilogy(t(var.start:var.end),data_PE(var.start:var.end),'go')
hold on
semilogy(t(var.start:var.end),irf_PE(var.start:var.end),'r-')
hold on
semilogy(t(var.start:var.end),G0_PE(var.start:var.end),'k-','LineWidth',2)
title('PERPENDICULAR (PE or VH)')
axis([0.5 12 mean(data_PE(var.start:var.start+10))/2 2*max(data_PE)])

subplot(3,2,[5,6]);
plot(t(var.start:var.end),data_PA(var.start:var.end)-G0_PA(var.start:var.end),'b-','LineWidth',2)
hold on
plot(t(var.start:var.end),data_PE(var.start:var.end)-G0_PE(var.start:var.end),'g-','LineWidth',2)
title('RESIDUALS')
%axis([0.5 12 min([min(data_PA(70:950)-G0_PA(70:950)) min(data_PE(70:950)-G0_PE(70:950))])...
%    max([max(data_PA(70:950)-G0_PA(70:950)) max(data_PE(70:950)-G0_PE(70:950))]) ])
set(findall(F3,'-property','FontSize'),'FontSize',18)
set(get(F3,'Children'),'LineWidth',2,'FontSize',20,'fontWeight','bold')

    print ('-r300', '-dpng', [var.pathname filename_s 'PAnPE.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%ANISOTROPY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
anis=(data_PA-(data_PE.*var.Gf))./(data_PA+2*(data_PE.*var.Gf)+eps);
anis_G0=(G0_PA-(G0_PE.*var.Gf))./(G0_PA+2*(G0_PE.*var.Gf)+eps);
A_start=find(irf_PA==max(irf_PA));
tt=t-t(A_start); L_tt=length(tt);
for ti=1:L_tt
    if tt(ti)<0
        tt(ti+L_tt)=tt(ti)+12.5;
        anis(ti+L_tt)=anis(ti);
        anis_G0(ti+L_tt)=anis_G0(ti);
    end
end
A_end=min(find(tt>10));

F4=figure('Name','ANISOTROPY DECAY','NumberTitle','off');
set(F4,'OuterPosition',[1 var.scrsz(4)/2 var.scrsz(3)/2 var.scrsz(4)]);%%@ SECOND VALUE FOR PC=1 FOR MAC=var.scrsz(4)/2
subplot(3,2,[1,2,3,4]);
plot(tt(A_start:A_end),anis(A_start:A_end),'ro')
hold on
plot(tt(A_start:A_end),anis_G0(A_start:A_end),'b-','LineWidth',2) 
title('ANISOTROPY DECAY')
axis([0 10 0 0.571])
subplot(3,2,[5,6]);
plot(tt(A_start:A_end),anis(A_start:A_end)-anis_G0(A_start:A_end),'r-','LineWidth',2)
title('RESIDUALS')
axis([0 10 min(anis(A_start:A_end)-anis_G0(A_start:A_end)) max(anis(A_start:A_end)-anis_G0(A_start:A_end))])
set(findall(F4,'-property','FontSize'),'FontSize',18)
set(get(F4,'Children'),'LineWidth',2,'FontSize',20,'fontWeight','bold')
    print ('-r300', '-dpng', [var.pathname filename_s 'A.png'])
end