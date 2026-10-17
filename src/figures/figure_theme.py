"""Shared physical-coordinate theme; selected data are read without inference."""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from palette import PRIMARY, NODE, INK, NS, SUB

ROOT = Path(os.environ['BCPN_FIGURE_ROOT']).resolve()
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
 'xtick.labelsize':8,'ytick.labelsize':8,'axes.linewidth':.5,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'path','text.color':INK,'axes.labelcolor':INK,
 'xtick.color':INK,'ytick.color':INK,'savefig.facecolor':'white'})
def rd(key): return pd.read_csv(ROOT/'source_data'/f'{key}.csv')
def q(x):
    if pd.isna(x): return 'NA'
    if x >= .001: return f'{x:.4f}' if x < .01 else f'{x:.3f}'
    e=int(np.floor(np.log10(x)))
    return rf'${x/10**e:.1f}\times10^{{{e}}}$'
def figure(h):
    f=plt.figure(figsize=(174/25.4,h/25.4)); f.mmheight=h; return f
def ax(f,x,y,w,h): return f.add_axes([x/174,1-(y+h)/f.mmheight,w/174,h/f.mmheight])
def text(f,x,y,s,size=8,weight='normal',color=INK,ha='left',va='center',**kw):
    return f.text(x/174,1-y/f.mmheight,s,fontsize=size,weight=weight,color=color,ha=ha,va=va,**kw)
def panel(f,x,y,s): text(f,x,y,s,11,'bold')
def clean(a):
    a.spines[['top','right']].set_visible(False); a.tick_params(width=.5,length=2.5)
