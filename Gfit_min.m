%% The Gfit_min function performs the fitting on all the data and outputs chi-squared value for minimization purposes
%Last adaption Thomas S van Zanten 160222
function [chi2] = Gfit_min(par0,lt,rt,t,IRF_T,IRF_PA,IRF_PE,T,PA,PE,par_lt,var)
m=0;%%%BACKGROUND IS ASSUMED TO BE CLOSE TO 0 PER CHANNEL
ll=length(t);%%%No of channels
rep=var.rep; ind_start=var.start; ind_end=var.end;
t1=[t(1:ll);t(1:ll)+rep;t(1:ll)+2*rep];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lt==1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE EXPONENTIAL LIFETIME FUNCTION WITH OFFSET%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(T)>0
    G=@(t,a,tau_lt) 3.*a.*exp(-t/tau_lt);
    G0=conv(IRF_T, G(t1,par0(1),par0(2)))+m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%taking into account photon overflow due to long lifetime vs repetition rate
    G0_T(1:ll)=G0(1:ll)+G0(ll+1:2*ll)+G0(2*ll+1:3*ll); G0_T=G0_T';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2= sum(sum( ((G0_T(ind_start:ind_end)-T(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%% 
    elseif rt==1
a=par_lt(1); tau_lt=par_lt(2);
%%%%%%%%%%%%%%SINGLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,tau1) (Nf1.*a).*exp(-t/tau_lt).*(1+2.*r0.*exp(-t/tau1)); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4)))+m/3;
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA';    
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%SINGLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,tau1) (Nf2.*a).*exp(-t/tau_lt).*(1-r0.*exp(-t/tau1));
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4)))+m/3;    
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';    
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOUBLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%%
    elseif rt==2
a=par_lt(1); tau_lt=par_lt(2);
%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,d,tau1,tau2) (Nf1.*a).*exp(-t/tau_lt).*(1+2.*(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)))); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4),par0(5),par0(6)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA'; 
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,d,tau1,tau2) (Nf2.*a).*exp(-t/tau_lt).*(1-(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)))); 
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4),par0(5),par0(6)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRIPLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%%
    elseif rt==3
a=par_lt(1); tau_lt=par_lt(2);
%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,d,e,tau1,tau2,tau3) (Nf1.*a).*exp(-t/tau_lt).*(1+2.*(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)+e.*exp(-t/tau3)))); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4),par0(7),par0(5),par0(6),par0(8)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA'; 
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,d,e,tau1,tau2,tau3) (Nf2.*a).*exp(-t/tau_lt).*(1-(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)+e.*exp(-t/tau3)))); 
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4),par0(7),par0(5),par0(6),par0(8)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif lt==2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOUBLE EXPONENTIAL LIFETIME FUNCTION WITH OFFSET%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(T)>0
    G=@(t,a,b,tau_lt1,tau_lt2) 3.*a.*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2));
    G0=conv(IRF_T, G(t1,par0(1),par0(2),par0(3),par0(4)))+m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%taking into account photon overflow due to long lifetime vs repetition rate
    G0_T(1:ll)=G0(1:ll)+G0(ll+1:2*ll)+G0(2*ll+1:3*ll); G0_T=G0_T';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2= sum(sum( ((G0_T(ind_start:ind_end)-T(ind_start:ind_end)).^2 )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%%
    elseif rt==1
a=par_lt(1); b=par_lt(2); tau_lt1=par_lt(3); tau_lt2=par_lt(4);
%%%%%%%%%%%%%%SINGLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,tau1) (Nf1.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)).*(1+2.*r0.*exp(-t/tau1)); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4)))+m/3;
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA';  
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%SINGLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,tau1) (Nf2.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)).*(1-r0.*exp(-t/tau1));
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4)))+m/3;    
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';   
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOUBLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%%
    elseif rt==2 && sum(PA)>0
a=par_lt(1); b=par_lt(2); tau_lt1=par_lt(3); tau_lt2=par_lt(4);
%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,d,tau1,tau2) (Nf1.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)).*(1+2.*r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2))); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4),par0(5),par0(6)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA';   
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,d,tau1,tau2) (Nf2.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)).*(1-r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2))); 
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4),par0(5),par0(6)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';   
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRIPLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%%
    elseif rt==3 && sum(PA)>0
a=par_lt(1); b=par_lt(2); tau_lt1=par_lt(3); tau_lt2=par_lt(4);
%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,d,e,tau1,tau2,tau3) (Nf1.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)).*(1+2.*r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)+e.*exp(-t/tau3))); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4),par0(7),par0(5),par0(6),par0(8)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA';   
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,d,e,tau1,tau2,tau3) (Nf2.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)).*(1-r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)+e.*exp(-t/tau3))); 
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4),par0(7),par0(5),par0(6),par0(8)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';   
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif lt==3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRIPLE EXPONENTIAL LIFETIME FUNCTION WITH OFFSET%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(T)>0
    G=@(t,a,b,c,tau_lt1,tau_lt2,tau_lt3) 3.*a.*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)+(c).*exp(-t/tau_lt3));
    G0=conv(IRF_T, G(t1,par0(1),par0(2),par0(3),par0(4),par0(5),par0(6)))+m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%taking into account photon overflow due to long lifetime vs repetition rate
    G0_T(1:ll)=G0(1:ll)+G0(ll+1:2*ll)+G0(2*ll+1:3*ll); G0_T=G0_T';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2= sum(sum( ((G0_T(ind_start:ind_end)-T(ind_start:ind_end)).^2 )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%% 
    elseif rt==1
a=par_lt(1); b=par_lt(2); c=par_lt(3); tau_lt1=par_lt(4); tau_lt2=par_lt(5); tau_lt3=par_lt(6);
%%%%%%%%%%%%%%SINGLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,tau1) (Nf1.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)+(c).*exp(-t/tau_lt3)).*(1+2.*r0.*exp(-t/tau1)); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4)))+m/3;
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA';    
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%SINGLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,tau1) (Nf2.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)+(c).*exp(-t/tau_lt3)).*(1-r0.*exp(-t/tau1));
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4)))+m/3;    
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';    
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOUBLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%%
    elseif rt==2
a=par_lt(1); b=par_lt(2); c=par_lt(3); tau_lt1=par_lt(4); tau_lt2=par_lt(5); tau_lt3=par_lt(6);
%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,d,tau1,tau2) (Nf1.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)+(c).*exp(-t/tau_lt3)).*(1+2.*(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)))); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4),par0(5),par0(6)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA'; 
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,d,tau1,tau2) (Nf2.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)+(c).*exp(-t/tau_lt3)).*(1-(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)))); 
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4),par0(5),par0(6)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRIPLE EXPONENTIAL ROTATIONAL FUNCTION WITH OVERFLOW%%%%%%%%%%%%%%%%%%%%
    elseif rt==3
a=par_lt(1); b=par_lt(2); c=par_lt(3); tau_lt1=par_lt(4); tau_lt2=par_lt(5); tau_lt3=par_lt(6);
%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PA%%%%%%%%%%%%%%%%%
    G_PA=@(t,Nf1,r0,d,e,tau1,tau2,tau3) (Nf1.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)+(c).*exp(-t/tau_lt3)).*(1+2.*(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)+e.*exp(-t/tau3)))); 
    G0a=conv(IRF_PA, G_PA(t1,par0(1),par0(3),par0(4),par0(7),par0(5),par0(6),par0(8)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PA(1:ll)=G0a(1:ll)+G0a(ll+1:2*ll)+G0a(2*ll+1:3*ll); G0_PA=G0_PA'; 
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PA= sum(sum( ((G0_PA(ind_start:ind_end)-PA(ind_start:ind_end)).^2 )));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%DOUBLE EXPONENTIAL FUNCTION WITH OVERFLOW PE%%%%%%%%%%%%%%%%%
    G_PE=@(t,Nf2,r0,d,e,tau1,tau2,tau3) (Nf2.*a).*(b.*exp(-t/tau_lt1)+(1-b).*exp(-t/tau_lt2)+(c).*exp(-t/tau_lt3)).*(1-(r0.*(d.*exp(-t/tau1)+(1-d).*exp(-t/tau2)+e.*exp(-t/tau3)))); 
    G0e=conv(IRF_PE, G_PE(t1,par0(2),par0(3),par0(4),par0(7),par0(5),par0(6),par0(8)))+m/3;  
%taking into account photon overflow due to long lifetime vs repetition rate
	G0_PE(1:ll)=G0e(1:ll)+G0e(ll+1:2*ll)+G0e(2*ll+1:3*ll); G0_PE=G0_PE';
%%%%%%%%%%%%%%MINIMALIZATION PROBLEM (BETWEEN DATA AND FIT)%%%%%%%%%%%%%%%%
chi2_PE= sum(sum( ((G0_PE(ind_start:ind_end)-PE(ind_start:ind_end)).^2 )));

chi2=chi2_PA+chi2_PE;
    end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end