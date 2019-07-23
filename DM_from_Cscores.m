load('../../../../ownCloud2/ForDropbox/euclid/skt_agarch11/skt_agarch11_-0.5_5_T2000_H100_II4_PCP_NEW_MULTI_THRES_MCLE_2.mat');
model = 'skt_gas'; %'skt_agarch11';
data_name = 'IBM_T2000_crisis_up';

addpath(genpath('include/'));

%% time constant evaluation
C_score_post = results0.C_score ;

C_score_Ca = results_a.C_score_C ;
C_score_PCa = results_a.C_score_PC; 

C_score_Cb = results_b.C_score_C ;
C_score_PCb = results_b.C_score_PC;

C_score_Ca2 = results_a2.C_score_C ;
C_score_PCa2 = results_a2.C_score_PC; 

C_score_Cb2 = results_b2.C_score_C ;
C_score_PCb2 = results_b2.C_score_PC;

C_score_Cc = results_c.C_score_C ;
C_score_PCc = results_c.C_score_PC;

C_score_Cd = results_d.C_score_C ;
C_score_PCd = results_d.C_score_PC; 
% 
% C_score_Ce = results_e.C_score_C ;
% C_score_PCe = results_e.C_score_PC;

C_score_C_z = results_z.C_score_C ;
C_score_PC_z = results_z.C_score_PC;



C_score_Cm_a = results_m_a.C_score_C ;
C_score_PCm_a = results_m_a.C_score_PC; 

C_score_Cm_a2 = results_m_a2.C_score_C ;
C_score_PCm_a2 = results_m_a2.C_score_PC; 

C_score_Cm_b2 = results_m_b2.C_score_C ;
C_score_PCm_b2 = results_m_b2.C_score_PC; 

C_score_Cm_b = results_m_b.C_score_C ;
C_score_PCm_b = results_m_b.C_score_PC; 

C_score_Cm_c = results_m_c.C_score_C ;
C_score_PCm_c = results_m_c.C_score_PC; 

C_score_Cm_d = results_m_d.C_score_C ;
C_score_PCm_d = results_m_d.C_score_PC; 

% C_score_Cm_e = results_m_e.C_score_C ;
% C_score_PCm_e = results_m_e.C_score_PC; 


%% time varying evaluation
Cv_score_post = results0.Cv_score;

Cv_score_C_z = results_z.Cv_score_C ;
Cv_score_PC_z = results_z.Cv_score_PC; 

Cv_score_Ca = results_a.Cv_score_C ;
Cv_score_PCa = results_a.Cv_score_PC; 

Cv_score_Cb = results_b.Cv_score_C ;
Cv_score_PCb = results_b.Cv_score_PC;

Cv_score_Ca2 = results_a2.Cv_score_C ;
Cv_score_PCa2 = results_a2.Cv_score_PC; 

Cv_score_Cb2 = results_b2.Cv_score_C ;
Cv_score_PCb2 = results_b2.Cv_score_PC;

Cv_score_Cc = results_c.Cv_score_C ;
Cv_score_PCc = results_c.Cv_score_PC;

Cv_score_Cd = results_d.Cv_score_C ;
Cv_score_PCd = results_d.Cv_score_PC; 

% Cv_score_Ce = results_e.Cv_score_C ;
% Cv_score_PCe = results_e.Cv_score_PC;


Cv_score_Cm_a = results_m_a.Cv_score_C ;
Cv_score_PCm_a = results_m_a.Cv_score_PC; 

Cv_score_Cm_a2 = results_m_a2.Cv_score_C ;
Cv_score_PCm_a2 = results_m_a2.Cv_score_PC; 

Cv_score_Cm_b = results_m_b.Cv_score_C;
Cv_score_PCm_b = results_m_b.Cv_score_PC; 

Cv_score_Cm_b2 = results_m_b2.Cv_score_C ;
Cv_score_PCm_b2 = results_m_b2.Cv_score_PC; 

Cv_score_Cm_c = results_m_c.Cv_score_C ;
Cv_score_PCm_c = results_m_c.Cv_score_PC; 

Cv_score_Cm_d = results_m_d.Cv_score_C ;
Cv_score_PCm_d = results_m_d.Cv_score_PC; 

% Cv_score_Cm_e = results_m_e.Cv_score_C ;
% Cv_score_PCm_e = results_m_e.Cv_score_PC; 

%%


dens_post = predictive_dens_skt_gas(y(T:(T+H)), results0.fT, results0.draw);     

dens_C_z = predictive_dens_skt_gas(y(T:(T+H)), results_z.fT_C, results_z.draw_C);     
dens_PC_z = predictive_dens_skt_gas(y(T:(T+H)), results_z.fT_PC, results_z.draw_PC);     

dens_Ca = predictive_dens_skt_gas(y(T:(T+H)), results_a.fT_C, results_a.draw_C);     
dens_PCa = predictive_dens_skt_gas(y(T:(T+H)), results_a.fT_PC, results_a.draw_PC);     

dens_Cb = predictive_dens_skt_gas(y(T:(T+H)), results_b.fT_C, results_b.draw_C);     
dens_PCb = predictive_dens_skt_gas(y(T:(T+H)), results_b.fT_PC, results_b.draw_PC); 

dens_Cc = predictive_dens_skt_gas(y(T:(T+H)), results_c.fT_C, results_c.draw_C);     
dens_PCc = predictive_dens_skt_gas(y(T:(T+H)), results_c.fT_PC, results_c.draw_PC); 

dens_Cd = predictive_dens_skt_gas(y(T:(T+H)), results_d.fT_C, results_d.draw_C);     
dens_PCd = predictive_dens_skt_gas(y(T:(T+H)), results_d.fT_PC, results_d.draw_PC); 

% dens_Ce = predictive_dens_skt_gas(y(T:(T+H)), results_e.fT_C, results_e.draw_C);     
% dens_PCe = predictive_dens_skt_gas(y(T:(T+H)), results_e.fT_PC, results_e.draw_PC); 

dens_Cm_a = predictive_dens_skt_gas(y(T:(T+H)), results_m_a.fT_C, results_m_a.draw_C);     
dens_PCm_a = predictive_dens_skt_gas(y(T:(T+H)), results_m_a.fT_PC, results_m_a.draw_PC);     

dens_Cm_a2 = predictive_dens_skt_gas(y(T:(T+H)), results_m_a2.fT_C, results_m_a2.draw_C);     
dens_PCm_a2 = predictive_dens_skt_gas(y(T:(T+H)), results_m_a2.fT_PC, results_m_a2.draw_PC);  

dens_Cm_b = predictive_dens_skt_gas(y(T:(T+H)), results_m_b.fT_C, results_m_b.draw_C);     
dens_PCm_b = predictive_dens_skt_gas(y(T:(T+H)), results_m_b.fT_PC, results_m_b.draw_PC); 

dens_Cm_b2 = predictive_dens_skt_gas(y(T:(T+H)), results_m_b2.fT_C, results_m_b2.draw_C);     
dens_PCm_b2 = predictive_dens_skt_gas(y(T:(T+H)), results_m_b2.fT_PC, results_m_b2.draw_PC); 

dens_Cm_c = predictive_dens_skt_gas(y(T:(T+H)), results_m_c.fT_C, results_m_c.draw_C);     
dens_PCm_c = predictive_dens_skt_gas(y(T:(T+H)), results_m_c.fT_PC, results_m_c.draw_PC); 

dens_Cm_d = predictive_dens_skt_gas(y(T:(T+H)), results_m_d.fT_C, results_m_d.draw_C);     
dens_PCm_d = predictive_dens_skt_gas(y(T:(T+H)), results_m_d.fT_PC, results_m_d.draw_PC); 



%% DM stats
aaa = who('-regexp','C_score\w*$');  
bbb = who('-regexp','Cv_score\w*$');  
ccc = who('-regexp','dens\w*$');  
met = length(aaa);

tmp = (met-1)/2;
ind = [met,repmat([1,(tmp+1)],1,tmp) + reshape(repmat((0:1:tmp-1),2,1),1,2*tmp) ];
% ind = [met,1,5,2,6,3,7,4,8]
aaa = aaa(ind);
bbb = bbb(ind);
ccc = ccc(ind);

DM = NaN(met,met,3);
DMv = NaN(met,met,3);

BF = NaN(met,met);
 
cBF = NaN(met,met,3);
cBFv = NaN(met,met,3);

% NWSE =  NaN(met,met,3);

for ii = 2:met
    for jj = 1:(ii-1)
        d = eval(char(aaa{ii})) - eval(char(aaa{jj}));
        for dd = 1:3
            nwse = sqrt(NeweyWest(d(dd,:)));
%             NWSE(ii,jj,dd) = nwse;                        
            DM(ii,jj,dd) = mean(d(dd,:))/nwse;            
        end
    end
end


for ii = 2:met
    for jj = 1:(ii-1)
        d = eval(char(bbb{ii})) - eval(char(bbb{jj}));
        for dd = 1:3
            nwse = sqrt(NeweyWest(d(dd,:)));
            DMv(ii,jj,dd) = mean(d(dd,:))/nwse;            
        end
    end
end




for ii = 2:met
    for jj = 1:(ii-1)
        d1 = eval(char(ccc{ii}));
        d2 = eval(char(ccc{jj}));            
        BF(ii,jj) = exp(sum(log(d1)) - sum(log(d2)));            
    end
end
 
BFone = (BF > 1) ;





for ii = 2:met
    for jj = 1:(ii-1)
        d1 = eval(char(aaa{ii}));
        d2 = eval(char(aaa{jj}));
        for dd = 1:3
            cBF(ii,jj,dd) = exp(sum(d1(dd)) - sum(d2(dd)));            
        end
    end
end


for ii = 2:met
    for jj = 1:(ii-1)
        d1 = eval(char(bbb{ii}));
        d2 = eval(char(bbb{jj}));
        for dd = 1:3
            cBFv(ii,jj,dd) = exp(sum(d1(dd)) - sum(d2(dd)));            
        end
    end
end


name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_DM_DMv_3.mat'];
save(name,'DM','DMv','aaa','bbb');

name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_BF_BFv_3.mat'];
save(name,'BF','BFone','cBF','cBFv','aaa','bbb','ccc');
