#!/usr/bin/env python3
"""
Specific-etiology detection by cancer type, eSS versus COSMIC v3.6.

For each etiology group, compares the fraction of samples in which the eSS panel
and the matched COSMIC signatures are detected, per cancer type. Also produces the
treatment-slice pie (Figure 2), the detection McNemar bars (Figure 3), and the
eSS-COSMIC decomposition panel (Figure 5).
"""
import os, glob
import numpy as np, pandas as pd
from math import floor, log10
from statsmodels.stats.contingency_tables import mcnemar
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE=os.environ.get("ESS_BASE", "../data/pancan_eSS_assignment")
OUT=os.path.join(BASE,"eSS_cancer_Etiology")
META=os.path.join(BASE,"HD Signatures WGS Datasets v2.0 - final_sample_summary_v2.tsv")
ESSF=os.path.join(BASE,"Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
COSMIC_DIR=os.path.join(BASE,"output_assignment_OGSBS5_cosmic_3_6")
COHORTS=["TCGA","PCAWG","Mutographs"]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})

ESS_ETIO_SETS={
 'Tobacco/PAH (eSS16-21,45)':['eSS16','eSS17','eSS18','eSS19','eSS20','eSS21','eSS45'],
 'DMBA/DBPDE (PAH, T>A) (eSS32,34)':['eSS32','eSS34'],
 'UV/solar (eSS8,26,27)':['eSS8','eSS26','eSS27'],
 'Aflatoxin (eSS10,11)':['eSS10','eSS11'],
 'Aristolochic acid (eSS33)':['eSS33'],
 'Colibactin (eSS36)':['eSS36'],
 'Acrylamide/glycidamide (eSS39)':['eSS39'],
 'Bromochloroacetic acid (eSS37)':['eSS37'],
 'Aromatic amine/diet (eSS9,13,22,44)':['eSS9','eSS13','eSS22','eSS44'],
 'Alkylating/indust. (eSS24,25,28,29,40-42,46,47)':['eSS24','eSS25','eSS28','eSS29','eSS40','eSS41','eSS42','eSS46','eSS47'],
 'Other chemical (eSS2,12,15)':['eSS2','eSS12','eSS15'],
 'Treatment (eSS1,14,23,30,31,35,43,48,49)':['eSS1','eSS14','eSS23','eSS30','eSS31','eSS35','eSS43','eSS48','eSS49'],
}
COS_ETIO_SETS={
 'Tobacco (SBS4,92,100)':['SBS4','SBS92','SBS100'],'Tobacco-chewing (SBS29)':['SBS29'],
 'UV/solar (SBS7,38)':['SBS7a','SBS7b','SBS7c','SBS7d','SBS38'],'Aflatoxin (SBS24)':['SBS24'],
 'Aristolochic acid (SBS22)':['SBS22a','SBS22b','SBS22c'],'Colibactin (SBS88)':['SBS88'],
 'Duocarmycin (SBS90)':['SBS90'],'Temozolomide (SBS11)':['SBS11'],'Platinum (SBS31,35)':['SBS31','SBS35'],
 'Chemo (SBS25,86,99)':['SBS25','SBS86','SBS99'],'Haloalkane (SBS42)':['SBS42'],
 'APOBEC (SBS2,13)':['SBS2','SBS13'],'MMR (SBS6,15,20,21,26,44)':['SBS6','SBS15','SBS20','SBS21','SBS26','SBS44'],
 'POLE/POLD (SBS10,14)':['SBS10a','SBS10b','SBS10c','SBS10d','SBS14'],'AID (SBS9,84,85)':['SBS9','SBS84','SBS85'],
 'Clock-like (SBS1,5,40)':['SBS1','SBS5','SBS40a','SBS40b','SBS40c'],
}
# treatment membership for fig2 pie typing
TRT_ESS=set(ESS_ETIO_SETS['Treatment (eSS1,14,23,30,31,35,43,48,49)'])
ENV_ESS=set().union(*[set(v) for k,v in ESS_ETIO_SETS.items() if k!='Treatment (eSS1,14,23,30,31,35,43,48,49)'])
FLAT_ESS={'eSS3','eSS4','eSS5','eSS6','eSS7'}

meta=pd.read_csv(META,sep='\t',dtype=str); meta=meta[meta.Cohort.isin(COHORTS)].copy()
ess=pd.read_csv(ESSF,sep='\t').set_index('Samples'); ess=ess[[c for c in ess.columns if c.startswith('eSS')]]
cos=pd.concat([pd.read_csv(f,sep='\t').set_index('Samples') for f in
    glob.glob(os.path.join(COSMIC_DIR,'*','Assignment_Solution','Activities','Assignment_Solution_Activities.txt'))])
cos=cos[~cos.index.duplicated()]
s=[x for x in meta.Sample_ID if x in ess.index and x in cos.index]
meta=meta.set_index('Sample_ID').loc[s]; ess=ess.loc[s]; cos=cos.loc[s]
er=ess.div(ess.sum(1).replace(0,np.nan),axis=0).fillna(0)
cr=cos.div(cos.sum(1).replace(0,np.nan),axis=0).fillna(0)
N=len(meta); print("Patients:",N)
ct=meta['Cancer_Type_HD_Grouping']

def fmtp(p):
    if p>=0.05: return 'ns'
    if p<1e-300: return r'$p<10^{-300}$'
    e=int(floor(log10(p))); m=p/10**e
    return rf'$p={m:.1f}\times10^{{{e}}}$'

# ============ FIG2v3 pie with visible treatment slice ============
cd=cr.idxmax(1); ed=er.idxmax(1)
def cosmic_flat(sig):
    return sig in COS_ETIO_SETS['Clock-like (SBS1,5,40)'] or sig in ('SBS39','SBS93','SBS94','SBS95','SBS12','SBS16','SBS17a','SBS17b','SBS19','SBS23','SBS28','SBS33','SBS34','SBS37','SBS41')
flat=pd.Series([cosmic_flat(x) for x in cd],index=cr.index)
etype=pd.Series(['treatment' if x in TRT_ESS else 'flat' if x in FLAT_ESS else 'environmental' for x in ed],index=er.index)
cnts={'Environmental exposure':int((flat&(etype=='environmental')).sum()),
      'Treatment':int((flat&(etype=='treatment')).sum()),
      'Background / clock-like':int((flat&(etype=='flat')).sum())}
fig,ax=plt.subplots(figsize=(7.6,5.8))
vals=list(cnts.values())
pielabs=[f'{k}\n{v:,} ({100*v/flat.sum():.1f}%)' if k!='Treatment' else '' for k,v in cnts.items()]
wedges,_=ax.pie(vals,labels=pielabs,colors=['#c0392b','#e67e22','#95a5a6'],startangle=90,
                labeldistance=1.13,wedgeprops=dict(edgecolor='white'))
wt=wedges[1]; ang=np.deg2rad((wt.theta1+wt.theta2)/2)
ax.annotate(f'Treatment\n{cnts["Treatment"]:,} ({100*cnts["Treatment"]/flat.sum():.1f}%)',
            xy=(np.cos(ang),np.sin(ang)),xytext=(1.5,0.05),fontsize=9,va='center',
            arrowprops=dict(arrowstyle='-',color='#e67e22',lw=1.2))
ax.set_title(f'eSS attribution of {int(flat.sum()):,} COSMIC clock-like patients',fontweight='bold')
plt.tight_layout()
for e in('png','pdf'): fig.savefig(os.path.join(OUT,f'fig2v3_rescue_pie.{e}'),bbox_inches='tight')
plt.close()

# ============ FIG3 positive controls — McNemar on detection, staggered brackets, two thresholds ============
def ev_nv(x):
    e=x.str.contains('current|ex-|ever|reformed|smoker|yes|^1$',case=False,regex=True,na=False)&\
      ~x.str.contains('never|non-?smoker|^no$|^0$',case=False,regex=True,na=False)
    n=x.str.contains('never|non-?smoker|lifelong',case=False,regex=True,na=False)
    return e,n
ever,never=ev_nv(meta['Smoking_Status'])
lung=ct.str.contains('Lung',case=False,na=False); skin=ct.str.contains('Skin',case=False,na=False)
TOB_E=ESS_ETIO_SETS['Tobacco/PAH (eSS16-21,45)']; TOB_S=COS_ETIO_SETS['Tobacco (SBS4,92,100)']
UV_E=ESS_ETIO_SETS['UV/solar (eSS8,26,27)']; UV_S=COS_ETIO_SETS['UV/solar (SBS7,38)']
def mcp(e,c,mask):
    a=int((e&c&mask).sum()); b=int((e&~c&mask).sum()); cc=int((~e&c&mask).sum()); d=int((~e&~c&mask).sum())
    return mcnemar([[a,b],[cc,d]],exact=True).pvalue
def build(THR):
    e_tob=er[TOB_E].sum(1)>=THR; c_tob=cr[[c for c in TOB_S if c in cr]].sum(1)>=THR
    e_uv=er[UV_E].sum(1)>=THR;   c_uv=cr[[c for c in UV_S if c in cr]].sum(1)>=THR
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    def panel(ax,title,rows):
        y=np.arange(len(rows)); h=0.36
        for yi,(lab_,e,c,mask) in enumerate(rows):
            ep,en=int((e&mask).sum()),int(mask.sum()); cp,cn=int((c&mask).sum()),int(mask.sum())
            ev=100*ep/en if en else 0; cv=100*cp/cn if cn else 0
            ax.barh(yi+h/2,ev,h,color='#2980b9'); ax.barh(yi-h/2,cv,h,color='#95a5a6')
            ltxt=f'{ev:.0f}% ({ep}/{en})'; rtxt=f'{cv:.0f}% ({cp}/{cn})'
            ax.text(ev+1,yi+h/2,ltxt,va='center',fontsize=8)
            ax.text(cv+1,yi-h/2,rtxt,va='center',fontsize=8,color='#555')
            # bracket placed past the longer of the two value-labels for this row (no overlap)
            bx=max(ev,cv)+max(len(ltxt),len(rtxt))*1.65+6
            ax.plot([bx,bx+1.5,bx+1.5,bx],[yi-h/2,yi-h/2,yi+h/2,yi+h/2],color='#444',lw=0.8)
            p=mcp(e,c,mask); ax.text(bx+3,yi,fmtp(p),va='center',fontsize=8,color='#b30000' if p<0.05 else '#777')
        ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.set_xlim(0,148)
        ax.set_xlabel(f'% patients with etiology detected (≥{int(THR*100)}%)'); ax.set_title(title,fontweight='bold',pad=8)
        ax.legend(handles=[Patch(facecolor='#2980b9',label='eSS'),Patch(facecolor='#95a5a6',label='COSMIC')],frameon=False,loc='upper right',fontsize=8)
    panel(axes[0],'Tobacco/PAH (eSS16-21,45 vs SBS4/92/100)',[
     ('Lung ever-smoker',e_tob,c_tob,lung&ever),('Lung never-smoker',e_tob,c_tob,lung&never),
     ('All ever-smoker',e_tob,c_tob,ever),('All never-smoker',e_tob,c_tob,never)])
    panel(axes[1],'UV/solar (eSS8,26,27 vs SBS7)',[('Skin',e_uv,c_uv,skin),('Non-skin',e_uv,c_uv,~skin)])
    fig.suptitle(f'Positive controls — eSS vs COSMIC detection (≥{int(THR*100)}%)',fontweight='bold',y=1.0)
    plt.tight_layout(rect=[0,0.02,1,0.97])
    tag='5pct' if THR==0.05 else '10pct'
    for e in('png','pdf'): fig.savefig(os.path.join(OUT,f'fig3v4_positive_controls_{tag}.{e}'),bbox_inches='tight')
    plt.close()
build(0.05); build(0.10)

# ============ FIG5v5 multi-label co-detection — relabelled; full + no-endogenous ============
edet=pd.DataFrame({g:(er[m].sum(1)>=0.05) for g,m in ESS_ETIO_SETS.items()})
cdet=pd.DataFrame({g:(cr[[c for c in m if c in cr]].sum(1)>=0.05) for g,m in COS_ETIO_SETS.items()})
ENDO={'APOBEC (SBS2,13)','MMR (SBS6,15,20,21,26,44)','POLE/POLD (SBS10,14)','AID (SBS9,84,85)'}
erow=list(ESS_ETIO_SETS)
def draw_codetect(ccol,fname,note):
    cnt=np.zeros((len(erow),len(ccol))); jac=np.zeros_like(cnt)
    for i,e in enumerate(erow):
        for j,c in enumerate(ccol):
            both=int((edet[e]&cdet[c]).sum()); uni=int((edet[e]|cdet[c]).sum())
            cnt[i,j]=both; jac[i,j]=both/uni if uni else 0
    fig,ax=plt.subplots(figsize=(max(10,len(ccol)*0.85),7.5))
    im=ax.imshow(jac,cmap='YlOrRd',aspect='auto',vmin=0)
    ax.set_xticks(range(len(ccol))); ax.set_xticklabels([f'{c}\n[n={int(cdet[c].sum()):,}]' for c in ccol],rotation=45,ha='right',fontsize=7)
    ax.set_yticks(range(len(erow))); ax.set_yticklabels([f'{e}  [n={int(edet[e].sum()):,}]' for e in erow],fontsize=7.5)
    for i in range(len(erow)):
        for j in range(len(ccol)):
            if cnt[i,j]>0: ax.text(j,i,int(cnt[i,j]),ha='center',va='center',fontsize=6.5,color='white' if jac[i,j]>0.25 else '#333')
    ax.set_xlabel('COSMIC etiology (signatures)  —  n = patients detected ≥5%'); ax.set_ylabel('eSS etiology (signatures)')
    ax.set_title('eSS vs COSMIC etiology co-detection (≥5%)  —  text = patients, colour = Jaccard overlap',fontweight='bold',fontsize=11)
    plt.colorbar(im,label='Jaccard (overlap / union)',shrink=.6)
    plt.tight_layout()
    for e in('png','pdf'): fig.savefig(os.path.join(OUT,f'{fname}.{e}'),bbox_inches='tight')
    plt.close()
    pd.DataFrame(cnt,index=erow,columns=ccol).astype(int).to_csv(os.path.join(OUT,f'{fname}_counts.csv'))
draw_codetect(list(COS_ETIO_SETS),'fig5v5_etiology_codetection','')
draw_codetect([c for c in COS_ETIO_SETS if c not in ENDO],'fig5v5_etiology_codetection_no_endogenous','')

# ============ FIG5v6 agreement (labels already carry eSS ids) ============
shared=[('Tobacco/PAH (eSS16-21,45)','Tobacco (SBS4,92,100)'),('UV/solar (eSS8,26,27)','UV/solar (SBS7,38)'),
        ('Aflatoxin (eSS10,11)','Aflatoxin (SBS24)'),('Aristolochic acid (eSS33)','Aristolochic acid (SBS22)'),
        ('Colibactin (eSS36)','Colibactin (SBS88)')]
fig,ax=plt.subplots(figsize=(9.5,5)); both=[];eo=[];co=[];labels=[]
for eg,cg in shared:
    e=edet[eg]; c=cdet[cg]; labels.append(f'{eg.split(" (")[0]}\n{eg.split("(")[1][:-1]} vs {cg.split("(")[1][:-1]}')
    both.append(int((e&c).sum())); eo.append(int((e&~c).sum())); co.append(int((~e&c).sum()))
y=np.arange(len(shared))
ax.barh(y,both,color='#27ae60',label='both + (agree)')
ax.barh(y,eo,left=both,color='#2980b9',label='eSS only +')
ax.barh(y,co,left=[b+e for b,e in zip(both,eo)],color='#95a5a6',label='COSMIC only +')
for yi in range(len(shared)):
    tot=both[yi]+eo[yi]+co[yi]
    rec=100*both[yi]/(both[yi]+co[yi]) if (both[yi]+co[yi]) else 0
    ppv=100*both[yi]/(both[yi]+eo[yi]) if (both[yi]+eo[yi]) else 0
    ax.text(tot+8,yi,f'agree={both[yi]}  recall={rec:.0f}%  PPV={ppv:.0f}%',va='center',fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels(labels,fontsize=8); ax.invert_yaxis()
ax.set_xlabel('patients detected (≥5%)'); ax.set_xlim(0,max(b+e+c for b,e,c in zip(both,eo,co))*1.6)
ax.set_title('Per-etiology detection agreement: eSS vs COSMIC (≥5%)',fontweight='bold',fontsize=11)
ax.legend(frameon=False,loc='lower right',fontsize=8)
plt.tight_layout()
for e in('png','pdf'): fig.savefig(os.path.join(OUT,f'fig5v6_etiology_agreement.{e}'),bbox_inches='tight')
plt.close()

# ============ FIG6 NEW: specific-etiology detection BY CANCER TYPE, eSS vs COSMIC ============
ETIO_TISSUE={'Colibactin':('Colibactin (eSS36)','Colibactin (SBS88)'),
             'Aristolochic acid':('Aristolochic acid (eSS33)','Aristolochic acid (SBS22)'),
             'Aflatoxin':('Aflatoxin (eSS10,11)','Aflatoxin (SBS24)'),
             'Tobacco/PAH':('Tobacco/PAH (eSS16-21,45)','Tobacco (SBS4,92,100)'),
             'UV/solar':('UV/solar (eSS8,26,27)','UV/solar (SBS7,38)')}
# FIXED cancer-type order across all panels (by frequency) so shared y-axis labels are consistent
vc=ct.value_counts(); types=vc[vc>=100].index.tolist()
yy=np.arange(len(types)); h=0.38
fig,axes=plt.subplots(1,len(ETIO_TISSUE),figsize=(4.2*len(ETIO_TISSUE),6.5),sharey=True)
for ax,(name,(eg,cg)) in zip(axes,ETIO_TISSUE.items()):
    e=edet[eg]; c=cdet[cg]
    evals=[int((e&(ct==t)).sum()) for t in types]; cvals=[int((c&(ct==t)).sum()) for t in types]
    ax.barh(yy+h/2,evals,h,color='#2980b9',label='eSS')
    ax.barh(yy-h/2,cvals,h,color='#95a5a6',label='COSMIC')
    ax.set_title(name,fontweight='bold',fontsize=10); ax.set_xlabel('patients ≥5%')
axes[0].set_yticks(yy); axes[0].set_yticklabels(types,fontsize=8); axes[0].invert_yaxis()
axes[0].legend(frameon=False,loc='lower right',fontsize=8)
fig.suptitle('Specific-etiology detection by cancer type — eSS vs COSMIC (≥5%)',fontweight='bold',y=1.02,fontsize=12)
plt.tight_layout()
for e in('png','pdf'): fig.savefig(os.path.join(OUT,f'fig6_specific_etiology_by_cancer.{e}'),bbox_inches='tight')
plt.close()

# ============ FIG1 (regenerated) dominant etiology class — single grey background, simple title ============
GTYPE={'Tobacco (SBS4,92,100)':'environmental','Tobacco-chewing (SBS29)':'environmental','UV/solar (SBS7,38)':'environmental',
 'Aflatoxin (SBS24)':'environmental','Aristolochic acid (SBS22)':'environmental','Colibactin (SBS88)':'environmental','Haloalkane (SBS42)':'environmental',
 'Duocarmycin (SBS90)':'treatment','Temozolomide (SBS11)':'treatment','Platinum (SBS31,35)':'treatment','Chemo (SBS25,86,99)':'treatment',
 'APOBEC (SBS2,13)':'endogenous','MMR (SBS6,15,20,21,26,44)':'endogenous','POLE/POLD (SBS10,14)':'endogenous','AID (SBS9,84,85)':'endogenous',
 'Clock-like (SBS1,5,40)':'flat'}
sig2type={sig:GTYPE[g] for g,sl in COS_ETIO_SETS.items() for sig in sl}
cos_dtype=pd.Series([sig2type.get(x,'flat') for x in cd],index=cr.index)
ess_dtype=etype  # environmental/treatment/flat
cats=['environmental','treatment','endogenous','flat']
clab={'environmental':'Environmental','treatment':'Treatment','endogenous':'Endogenous (COSMIC only)','flat':'Clock-like / background'}
ccol2={'environmental':'#c0392b','treatment':'#e67e22','endogenous':'#2c7fb8','flat':'#95a5a6'}
fig,ax=plt.subplots(figsize=(7,5))
for i,(name,ser) in enumerate([('eSS',ess_dtype),('COSMIC',cos_dtype)]):
    bottom=0
    for cclass in cats:
        pct=100*(ser==cclass).mean(); n=int((ser==cclass).sum())
        if pct==0: continue
        ax.bar(i,pct,0.55,bottom=bottom,color=ccol2[cclass],edgecolor='white')
        if pct>3: ax.text(i,bottom+pct/2,f'{pct:.1f}%\n(n={n:,})',ha='center',va='center',fontsize=8,
                          color='white' if cclass!='flat' else '#333')
        bottom+=pct
ax.set_xticks([0,1]); ax.set_xticklabels(['eSS','COSMIC'],fontsize=12)
ax.set_ylabel('% of patients'); ax.set_ylim(0,100)
ax.set_title(f'Dominant signature etiology (n={N:,})',fontweight='bold')
ax.legend(handles=[Patch(facecolor=ccol2[c],label=clab[c]) for c in cats],loc='center left',bbox_to_anchor=(1.01,0.5),frameon=False)
plt.tight_layout()
for e in('png','pdf'): fig.savefig(os.path.join(OUT,f'fig1v2_dominant_class_pooled.{e}'),bbox_inches='tight')
plt.close()

# ============ FIG7 PPV + recall per etiology (eSS vs COSMIC reference, >=5%) ============
shared7=[('Tobacco/PAH','Tobacco/PAH (eSS16-21,45)','Tobacco (SBS4,92,100)'),
         ('UV/solar','UV/solar (eSS8,26,27)','UV/solar (SBS7,38)'),
         ('Aflatoxin','Aflatoxin (eSS10,11)','Aflatoxin (SBS24)'),
         ('Aristolochic acid','Aristolochic acid (eSS33)','Aristolochic acid (SBS22)'),
         ('Colibactin','Colibactin (eSS36)','Colibactin (SBS88)')]
names=[]; ppv=[]; rec=[]
for nm,eg,cg in shared7:
    e=edet[eg]; c=cdet[cg]; tp=int((e&c).sum()); fp=int((e&~c).sum()); fn=int((~e&c).sum())
    names.append(nm); ppv.append(100*tp/(tp+fp) if (tp+fp) else 0); rec.append(100*tp/(tp+fn) if (tp+fn) else 0)
y=np.arange(len(names)); h=0.38
fig,ax=plt.subplots(figsize=(8.5,4.8))
ax.barh(y+h/2,ppv,h,color='#8e44ad',label='PPV (precision)')
ax.barh(y-h/2,rec,h,color='#16a085',label='Recall (sensitivity)')
for yi in range(len(names)):
    ax.text(ppv[yi]+1,yi+h/2,f'{ppv[yi]:.0f}%',va='center',fontsize=8)
    ax.text(rec[yi]+1,yi-h/2,f'{rec[yi]:.0f}%',va='center',fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(names); ax.invert_yaxis(); ax.set_xlim(0,105)
ax.set_xlabel('% (COSMIC detection as reference)')
ax.set_title('eSS precision (PPV) and recall vs COSMIC, per etiology (≥5%)',fontweight='bold',fontsize=11)
ax.legend(frameon=False,loc='lower right')
plt.tight_layout()
for e in('png','pdf'): fig.savefig(os.path.join(OUT,f'fig7_PPV_recall.{e}'),bbox_inches='tight')
plt.close()
print("v5 figures done (incl fig1 regen + fig7 PPV).")
