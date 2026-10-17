"""Display selected results only: no score construction, fitting or hypothesis tests."""
from figure_theme import ROOT, rd, q, figure, ax, text, panel, clean, PRIMARY, NODE, INK, NS, SUB
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from palette import SUPPORT, BURDEN, TRANSITION
from PIL import Image
import json, sys

plt.rcParams.update({'font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'hatch.linewidth':.5})
FINAL='--final' in sys.argv
def head(f,l,title,x,y):
    panel(f,x,y,l);text(f,x+7,y,title,9,'bold')
def save(f,name,supp=False):
    d=ROOT/('supplementary_figures' if supp else 'main_figures');base=d/name
    f.savefig(str(base)+'.pdf',dpi=600);f.savefig(str(base)+'.svg',dpi=600);f.savefig(str(base)+'.png',dpi=200)
    if FINAL:
        f.savefig(str(base)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(str(base)+'.tiff') as im:im.convert('RGB').save(str(base)+'.tiff',compression='tiff_lzw',dpi=(600,600))
    plt.close(f)
def paired(a,p,n,title,qval,delta=None):
    p=np.asarray(p,float);n=np.asarray(n,float)
    for u,v in zip(p,n):a.plot([0,1],[u,v],color='#BFC3C7',alpha=.7,lw=.45,zorder=1)
    for i,(v,c) in enumerate([(p,PRIMARY),(n,NODE)]):
        a.scatter(np.full(len(v),i),v,s=5,color=c,alpha=.7,zorder=2)
        b=a.boxplot(v,positions=[i],widths=.43,showfliers=False,patch_artist=True,
          medianprops={'color':INK,'linewidth':.9},boxprops={'facecolor':'white','edgecolor':c},
          whiskerprops={'color':c,'linewidth':.6},capprops={'color':c,'linewidth':.6})
        for box in b['boxes']:box.set_alpha(.65)
    a.set_xticks([0,1],['Primary','Node']);a.set_xlim(-.5,1.5);clean(a)
    st='q = '+q(qval)
    if delta is not None:st=rf'$\Delta$ = {delta:+.2f}; '+st
    a.set_title(title+'\n'+st,fontsize=8,pad=4)

def fig1():
    f=figure(210);d=rd('bulk')
    head(f,'a','Paired differential expression',4,5);a=ax(f,17,17,62,69)
    for mask,c in [(d['adj.P.Val']>=.05,'#CFD3D6'),((d['adj.P.Val']<.05)&(d.logFC<0),PRIMARY),((d['adj.P.Val']<.05)&(d.logFC>0),NODE)]:
        a.scatter(d.loc[mask,'logFC'],-np.log10(d.loc[mask,'P.Value']),s=3,c=c,linewidths=0)
    a.set_xlabel('Log2 fold change (Node - Primary)');a.set_ylabel('-log10 P');clean(a)
    a.text(.03,.98,'12,502 genes; 36 pairs\n562 primary-higher\n328 node-higher',transform=a.transAxes,va='top',fontsize=8)
    head(f,'b','Paired expression differences',88,5)
    hm=pd.read_csv(ROOT/'source_data/heatmap_paired_deltas.csv',index_col=0)
    a=ax(f,110,17,56,69);lim=np.max(np.abs(hm.values))
    im=a.imshow(hm,aspect='auto',cmap='RdBu',vmin=-lim,vmax=lim,interpolation='nearest')
    a.set_yticks(range(len(hm)),hm.index,fontsize=8);a.set_xticks([0,17,35],['1','18','36']);a.set_xlabel('Matched patient (display index)');a.tick_params(length=0)
    cb=f.colorbar(im,cax=ax(f,113,100,49,2.5),orientation='horizontal');cb.ax.tick_params(labelsize=8,length=2,pad=1)
    cb.set_ticks([-lim,0,lim]);cb.set_ticklabels([f'{-lim:.2f}','0',f'{lim:.2f}'])
    cb.set_label('Node - Primary (log2 expression)',fontsize=8,labelpad=1)
    head(f,'c','Candidate-set enrichment',4,117)
    o=rd('fig1_ora');a=ax(f,45,133,34,52)
    labels=['ECM organization [R]','Collagen metabolism [G]','Cell-substrate\nadhesion [G]','ECM-receptor [K]','B-cell activation [G]','B-cell receptor [G]','Adaptive immunity [G]','DNA replication [K]']
    z=a.scatter(100*o.gene_ratio_numeric,range(8),s=20+o.gene_count*1.2,c=-np.log10(o.BH_q),cmap='Blues',edgecolor=INK,lw=.3)
    a.set_yticks(range(8),labels,fontsize=8);a.set_ylim(7.6,-.7);a.set_xlim(0,21);a.set_xlabel('Query overlap (%)');a.axhline(3.5,color=NS,lw=.5,ls=':');clean(a)
    a.text(.99,1.025,'Primary-higher',transform=a.transAxes,ha='right',va='bottom',fontsize=8,color=PRIMARY)
    a.text(.99,.48,'Node-higher',transform=a.transAxes,ha='right',va='top',fontsize=8,color=NODE)
    handles=[a.scatter([],[],s=20+v*1.2,color='#BDD7E7',edgecolor=INK,lw=.3,label=str(v)) for v in [10,40,70]]
    f.legend(handles=handles,title='Gene count',loc='lower left',bbox_to_anchor=(.01,.0),ncol=3,frameon=False,fontsize=8,title_fontsize=8,handletextpad=.3,columnspacing=.7)
    cb=f.colorbar(z,cax=ax(f,48,202,31,2),orientation='horizontal');cb.ax.tick_params(labelsize=8,length=2,pad=1);cb.ax.set_title('-log10 BH q',fontsize=8,pad=2)
    head(f,'d','Ranked enrichment',88,117);g=rd('fig1_gsea');a=ax(f,127,133,28,52)
    names=['B-cell receptor [G]','E2F targets [H]','G2M checkpoint [H]','ECM degradation [R]','EMT [H]','Adipogenesis [H]']
    a.axvline(0,color=NS,lw=.5)
    for i,row in g.iterrows():a.scatter(row.NES,i,s=25,color=NODE if row.NES>0 else PRIMARY)
    a.set_yticks(range(6),names,fontsize=8);a.set_ylim(5.6,-.6);a.set_xlim(-4.3,2.8);a.set_xticks([-4,0,2]);a.set_xlabel('NES');clean(a)
    text(f,164,128,'BH q',8,'bold',ha='center')
    for i,row in g.iterrows():text(f,164,133+(i+.6)/6.2*52,q(row.padj),8,ha='center')
    save(f,'Fig1_final')

def fig2():
    f=figure(157);pairs=rd('pairs');est=rd('estimate_samples').set_index('score');stats=rd('estimate').set_index('score')
    head(f,'a','Paired ESTIMATE scores',4,5)
    for i,(key,label) in enumerate([('ImmuneScore','Immune'),('StromalScore','Stromal')]):
        a=ax(f,16+i*39,23,29,43);z=stats.loc[key];paired(a,est.loc[key,pairs.primary],est.loc[key,pairs.node],label,z.BH_q,z.median_delta)
        if i==0:a.set_ylabel('ESTIMATE score')
    head(f,'b','MCP-counter differences',89,5);m=rd('mcp');a=ax(f,120,18,30,54)
    a.set_xlim(-.65,1.30)
    for i,z in m.iterrows():
        if pd.notna(z.median_delta):
            a.plot([0,z.median_delta],[i,i],color='#BCC2C6',lw=.5)
            a.scatter(z.median_delta,i,s=23,facecolor=INK if z.BH_q<.05 else 'white',edgecolor=INK)
        else:a.text(0,i,'NA',fontsize=8,va='center')
    a.axvline(0,color=NS,lw=.5);a.set_yticks(range(10),['T cells','Cytotoxic','B lineage','NK cells','Monocytic','Myeloid DCs','Neutrophils','Endothelial','Fibroblasts','CD8 T cells'],fontsize=8);a.set_ylim(9.5,-.5);a.set_xlabel('Median paired Δ',fontsize=8);clean(a)
    text(f,164,13,'BH q',8,'bold',ha='center')
    for i,z in m.iterrows():text(f,164,18+(i+.5)/10*54,q(z.BH_q),8,ha='center')
    text(f,94,83,'Filled: q < 0.05; open: q ≥ 0.05',8)
    head(f,'c','Nonmalignant support for nine context-compatible genes',4,92)
    genes=rd('class_b').gene_symbol.tolist();n=rd('nonmalignant');sub=n[n.gene_symbol.isin(genes)].copy()
    assert sub.N_status.value_counts().to_dict()=={'N0':52,'not_evaluable':28,'N+':10}
    pops=['B','CAFs','CD4+ T','CD8+ T','DCs','Macrophage','NK','TECs','plasma B']
    assert len(pops)==9 and 'Mast cells' not in pops
    sub=sub[sub.population.isin(pops)].copy()
    assert sub.N_status.value_counts().to_dict()=={'N0':52,'not_evaluable':19,'N+':10}
    a=ax(f,35,101,127,36)
    for i,g in enumerate(genes):
        for j,p in enumerate(pops):
            state=sub[(sub.gene_symbol==g)&(sub.population==p)].iloc[0].N_status
            a.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=SUPPORT if state=='N+' else ('white' if state=='N0' else '#DFE2E5'),edgecolor='#A9AFB3',lw=.35,hatch='////' if state=='not_evaluable' else None))
    a.set_xlim(-.5,len(pops)-.5);a.set_ylim(len(genes)-.5,-.5);a.set_yticks(range(len(genes)),genes,fontsize=8)
    a.set_xticks(range(len(pops)),['B','CAFs','CD4 T','CD8 T','DCs','Mac.','NK','TECs','Plasma B'],fontsize=8,rotation=35,ha='right');a.tick_params(length=0)
    handles=[Patch(facecolor=SUPPORT,label='Supported (10)'),Patch(facecolor='white',edgecolor=NS,label='Evaluable, not supported (52)'),Patch(facecolor='#DFE2E5',edgecolor=NS,hatch='////',label='Not evaluable (19)')]
    f.legend(handles=handles,loc='lower center',bbox_to_anchor=(.52,.015),ncol=3,frameon=False,fontsize=8,columnspacing=1)
    sub.to_csv(ROOT/'source_data/Fig2_support_states.csv',index=False);save(f,'Fig2_final')

def coverage_supplement():
    genes=rd('class_b').gene_symbol.tolist();n=rd('nonmalignant')
    pops=n.population.drop_duplicates().tolist();sub=n[n.gene_symbol.isin(genes)].copy()
    assert sub.N_status.value_counts().to_dict()=={'N0':52,'not_evaluable':28,'N+':10}
    f=figure(76);head(f,'a','Full nonmalignant support coverage, including ineligible Mast cells',4,5)
    a=ax(f,35,15,127,39)
    for i,g in enumerate(genes):
        for j,p in enumerate(pops):
            state=sub[(sub.gene_symbol==g)&(sub.population==p)].iloc[0].N_status
            a.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=SUPPORT if state=='N+' else ('white' if state=='N0' else '#DFE2E5'),edgecolor='#A9AFB3',lw=.35,hatch='////' if state=='not_evaluable' else None))
    a.set_xlim(-.5,9.5);a.set_ylim(8.5,-.5);a.set_yticks(range(9),genes,fontsize=8)
    a.set_xticks(range(10),['B','CAFs','CD4 T','CD8 T','DCs','Mac.','Mast','NK','TECs','Plasma B'],fontsize=8,rotation=35,ha='right');a.tick_params(length=0)
    f.legend(handles=[Patch(facecolor=SUPPORT,label='Supported (10)'),Patch(facecolor='white',edgecolor=NS,label='Evaluable, not supported (52)'),Patch(facecolor='#DFE2E5',edgecolor=NS,hatch='////',label='Not evaluable (28)')],loc='lower center',bbox_to_anchor=(.52,.01),ncol=3,frameon=False,fontsize=8)
    save(f,'FigS3_full_coverage',True)
    sub.to_csv(ROOT/'source_data/FigS3_full_coverage.csv',index=False)

def alluvial(a,p):
    levels=list(SUB);flows=pd.crosstab(p.Primary_PAM50,p.Node_PAM50).reindex(index=levels,columns=levels,fill_value=0)
    left=flows.sum(axis=1);right=flows.sum(axis=0);lb={};rb={};pos=0
    for s in levels:lb[s]=pos;pos+=left[s]+1
    pos=0
    for s in levels:rb[s]=pos;pos+=right[s]+1
    lo=lb.copy();ro=rb.copy()
    for s in levels:
        for t in levels:
            nn=int(flows.loc[s,t])
            if not nn:continue
            xx=np.linspace(.2,.8,80);u=(xx-.2)/.6;yy=lo[s]+(ro[t]-lo[s])*(3*u*u-2*u*u*u)
            a.add_patch(Polygon(np.column_stack([np.r_[xx,xx[::-1]],41-np.r_[yy,(yy+nn)[::-1]]]),facecolor=SUB[s],alpha=.55,edgecolor='white',lw=.3));lo[s]+=nn;ro[t]+=nn
    for counts,base,x in [(left,lb,.17),(right,rb,.8)]:
        for s,sh in zip(levels,['LumA','LumB','HER2-E','Basal','Normal-like']):
            nn=int(counts[s])
            if nn:
                a.add_patch(Rectangle((x,41-base[s]-nn),.03,nn,color=SUB[s]));a.text(x-.015 if x<.5 else x+.045,41-base[s]-nn/2,f'{sh} ({nn})',ha='right' if x<.5 else 'left',va='center',fontsize=8)
    a.set_xlim(0,1);a.set_ylim(0,42);a.axis('off')

def fig3():
    f=figure(234);p=rd('pam50');r=rd('risk').set_index('signature');v=rd('readout_patients');s=rd('readouts').set_index('score')
    head(f,'a','Research PAM50 (43/50 features)',4,5);a=ax(f,8,16,155,44);alluvial(a,p)
    text(f,36,13,'Primary',8,ha='center');text(f,135,13,'Node',8,ha='center');z=r.loc['PAM50']
    text(f,87,64,f'28 concordant; 8/36 discordant ({z.discordance_pct:.1f}%); '+rf'$\kappa$ = {z.kappa:.3f}',8,ha='center')
    for l,y,title,keys,names in [('b',73,'Programs (mean gene z-score)',['ER_EARLY','ER_LATE','E2F','G2M'],['ER early','ER late','E2F','G2M']),('c',124,'Transcripts (processed log2 expression)',['ESR1','PGR','ERBB2','MKI67'],['ESR1','PGR','ERBB2','MKI67'])]:
        head(f,l,title,4,y)
        for i,(key,label) in enumerate(zip(keys,names)):
            paired(ax(f,13+i*41, y+15,30,27),v[key+'_Primary'],v[key+'_Node'],label,s.loc[key,'BH_q'])
    patient=rd('risk_patients');records=[]
    for ix,(l,key) in enumerate(zip('def',['ROR-S-like','GGI-like','RS-like'])):
        x=4+ix*57;head(f,l,key,x,176);z=r.loc[key]
        levels=['Low','Intermediate','High'];mat=pd.crosstab(patient[key+'_Primary_category'],patient[key+'_Node_category']).reindex(index=levels,columns=levels,fill_value=0)
        assert mat.to_numpy().sum()==36 and np.trace(mat)==int(z.concordant_N)
        assert np.triu(mat.to_numpy(),1).sum()==int(z.upward_shift_N)
        a=ax(f,x+14,185,36,27)
        for i in range(3):
            for j in range(3):
                col=TRANSITION['same' if i==j else ('up' if j>i else 'down')]
                a.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=col,edgecolor='white',lw=.8))
                a.text(j,i,str(mat.iloc[i,j]),ha='center',va='center',fontsize=10,weight='bold')
                records.append({'signature':key,'primary_stratum':['lower','middle','upper'][i],'node_stratum':['lower','middle','upper'][j],'patient_count':int(mat.iloc[i,j])})
        a.set_xlim(-.5,2.5);a.set_ylim(2.5,-.5);a.set_xticks(range(3),['Lower','Middle','Upper'],fontsize=8);a.set_yticks(range(3),['Lower','Middle','Upper'],fontsize=8);a.tick_params(length=0)
        a.set_xlabel('Node stratum',fontsize=8,labelpad=1)
        if ix==0:a.set_ylabel('Primary stratum',fontsize=8,labelpad=1)
        text(f,x+29,182,f'Different: {int(z.discordant_N)}/36 ({z.discordance_pct:.1f}%)',8,ha='center')
        text(f,x+29,223,rf'Continuous Δ = {z.median_delta:.3f}; q = '+q(z.BH_q),8,ha='center')
    f.legend(handles=[Patch(facecolor=TRANSITION[k],edgecolor=NS,label=v) for k,v in [('same','Unchanged stratum'),('up','Upward stratum'),('down','Downward stratum')]],loc='lower center',bbox_to_anchor=(.52,.002),ncol=3,frameon=False,fontsize=8)
    pd.DataFrame(records).to_csv(ROOT/'source_data/Fig3_transition_counts.csv',index=False);save(f,'Fig3_final')

ORDER=['BCR','E2F','G2M','ECM','Adipogenesis']
def burden_boxes(f,y):
    v=rd('burden_patients')
    for i,key in enumerate(ORDER):
        a=ax(f,14+i*32,y,25,35);arr=[v.loc[v.nodal_group==g,key].dropna().values for g in ['N0','LOW','HIGH']]
        assert list(map(len,arr))==[983,600,307]
        bx=a.boxplot(arr,showfliers=False,patch_artist=True,widths=.65,medianprops={'color':INK,'linewidth':.8},whiskerprops={'linewidth':.6},capprops={'linewidth':.6})
        for p,c in zip(bx['boxes'],BURDEN):p.set_facecolor(c);p.set_linewidth(.5)
        a.set_xticks([1,2,3],['N0','1–3','≥4']);a.set_title(key,fontsize=8);clean(a)
        if i==0:a.set_ylabel('Mean gene z-score',fontsize=8,labelpad=1)
def fig4():
    f=figure(125);head(f,'a','Program distributions by positive-node burden',4,5);burden_boxes(f,18)
    text(f,87,62,'METABRIC primary tumors: N0, n = 983; 1–3, n = 600; ≥4, n = 307',8,ha='center')
    head(f,'b','Unadjusted association with continuous positive-node count',4,72)
    b=rd('burden').set_index('score').loc[ORDER];a=ax(f,37,83,74,30)
    a.axvline(0,color=NS,lw=.6)
    for i,(_,z) in enumerate(b.iterrows()):a.scatter(z.rho,i,s=24,facecolor=INK if z.trend_BH_q<.05 else 'white',edgecolor=INK)
    a.set_yticks(range(5),ORDER);a.set_ylim(4.6,-.6);a.set_xlim(-.14,.14);a.set_xticks([-.1,0,.1]);a.set_xlabel('Spearman ρ');clean(a)
    text(f,126,78,'ρ',8,'bold',ha='center');text(f,151,78,'BH q',8,'bold',ha='center')
    for i,(_,z) in enumerate(b.iterrows()):
        yy=83+(i+.6)/5.2*30;text(f,126,yy,f'{z.rho:+.3f}',8,ha='center');text(f,151,yy,q(z.trend_BH_q),8,ha='center')
    save(f,'Fig4_final')

def supplements():
    f=figure(171);pairs=rd('pairs');raw=rd('estimate_samples').set_index('score');s=rd('estimate').set_index('score')
    head(f,'a','Additional ESTIMATE outputs',4,5)
    for i,(key,label) in enumerate([('ESTIMATEScore','ESTIMATE score'),('TumorPurity','Estimated purity')]):
        paired(ax(f,23+i*80,21,50,43),raw.loc[key,pairs.primary],raw.loc[key,pairs.node],label,s.loc[key,'BH_q'])
    head(f,'b','Integrated candidate classification',4,81);c=rd('class_counts');c=c[c.bulk_direction_group=='ALL'];a=ax(f,16,93,55,46)
    a.bar(c.integrated_class,c.N,color=['#4F849F','#93B2B7','#C0C6CA','#E4E7E9'],edgecolor=INK,lw=.4);a.set_ylim(0,570);a.set_ylabel('Candidates');clean(a)
    for i,z in enumerate(c.itertuples()):a.text(i,z.N+10,str(z.N),ha='center',fontsize=8)
    head(f,'c','Nonmalignant evaluability',88,81);n=rd('nonmalignant');tab=pd.crosstab(n.population,n.N_status).reindex(columns=['N+','N0','not_evaluable'],fill_value=0);a=ax(f,117,92,47,47);left=np.zeros(len(tab))
    for key,col,h in [('N+',NODE,None),('N0','white',None),('not_evaluable','#DFE2E5','////')]:
        a.barh(range(len(tab)),tab[key],left=left,color=col,edgecolor=NS,lw=.3,hatch=h);left+=tab[key]
    a.set_yticks(range(len(tab)),tab.index,fontsize=8);a.invert_yaxis();a.set_xlabel('Candidate-population records');clean(a)
    text(f,9,152,'A: malignant-compatible; B: context-compatible; C: unresolved; D: malignant non-evaluable.',8)
    text(f,9,159,'Blue: supported; white: evaluable, not supported; hatched: not evaluable.',8)
    save(f,'FigS3_additional_panels',True)
    f=figure(94);head(f,'a','Exploratory ordinal nodal-burden associations',4,5);o=rd('ordinal').set_index('score').loc[ORDER];a=ax(f,39,18,73,49)
    a.axvline(1,color=NS,lw=.5)
    for i,(_,z) in enumerate(o.iterrows()):
        a.plot([z.CI95_low,z.CI95_high],[i,i],color=INK,lw=.7);a.scatter(z.OR_per_SD,i,s=25,color=NODE)
    a.set_yticks(range(5),ORDER);a.set_ylim(4.6,-.6);a.set_xlabel('Odds ratio per score SD (95% CI)');clean(a)
    text(f,142,12,'BH q',8,'bold',ha='center')
    for i,(_,z) in enumerate(o.iterrows()):text(f,142,18+(i+.6)/5.2*49,q(z.BH_q),8,ha='center')
    text(f,9,81,'Univariable exploratory models; n = 1,890. Proportional-odds assumption not checked.',8)
    save(f,'FigS6_ordinal_context',True)

if __name__=='__main__':
    for fn in [fig1,fig2,fig3,fig4,coverage_supplement,supplements]:fn()
