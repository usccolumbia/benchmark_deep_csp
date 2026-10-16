import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R='results/new_t180/'; F='P623_Deep_CSP_benchmark/figures/'
t=pd.read_csv(R+'t180_main_table.csv').sort_values('top1',ascending=False)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(10.4,3.8)); x=np.arange(len(t)); w=.27
cols=['#0072B2','#E69F00','#009E73']
for i,(c,l) in enumerate([('top1_pct','StructureMatcher'),('sg_pct','Space group'),('cons_pct','Consensus')]):
    b=ax.bar(x+(i-1)*w,t[c],w,label=l,color=cols[i],edgecolor='white',linewidth=.4)
    for xi,v in zip(x+(i-1)*w,t[c]): ax.text(xi,v+.8,f'{v:.0f}',ha='center',va='bottom',fontsize=6.5)
ax.set_xticks(x); ax.set_xticklabels(t.algorithm,rotation=30,ha='right'); ax.set_ylabel('Top-1 success rate (%)'); ax.set_ylim(0,80); ax.legend(frameon=False,ncol=3,loc='upper right')
fig.tight_layout(); fig.savefig(F+'T180_new_benchmark_metric.png',dpi=220); plt.close(fig)
for m,title,fn in [('sm1','StructureMatcher','tcsp_comparison_is_sm_match_new.png'),('sg1','Space Group','tcsp_comparison_is_sg_match_new.png'),('cons1','Consensus','tcsp_comparison_is_consensus_match_new.png')]:
    d=pd.read_csv(R+f't180_decomp_{m}.csv'); d['tot']=d.both+d.algo_only; d=d.sort_values('tot',ascending=False).reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(10.6,4.8)); bottom=np.zeros(len(d)); x=np.arange(len(d))
    for c,l,col in [('both','Both methods succeed','#43A047'),('algo_only','Compared CSP method only','#1E88E5'),('tcsp_only','TCSPv3 only','#FFB300'),('neither','Both methods fail','#F44336')]:
        ax.bar(x,d[c],0.5,bottom=bottom,color=col,label=l); bottom+=d[c].values
    ax.set_xticks(x); ax.set_xticklabels(d.algorithm,rotation=45,ha='right'); ax.set_ylim(0,185); ax.set_yticks(range(0,176,25))
    ax.set_ylabel('Count of targets'); ax.set_xlabel('Compared CSP Method'); ax.set_title(f'{title} Agreement with TCSPv3 (top-1)')
    ax.grid(axis='y',linestyle='--',alpha=.5); ax.set_axisbelow(True)
    for sp in ('top','right'): ax.spines[sp].set_visible(True)
    h,l=ax.get_legend_handles_labels(); ax.legend(h,l,title='Comparison Outcome',loc='center left',bbox_to_anchor=(1.01,.5),frameon=True,fontsize=8)
    fig.tight_layout(); fig.savefig(F+fn,dpi=220); plt.close(fig)
print('ok')
