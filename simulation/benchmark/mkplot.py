import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, warnings; warnings.filterwarnings("ignore")

d=pd.read_csv(os.environ.get("RESULTS","../results/bench_results.csv")); d=d[d.method!="method"]
for c in ["snr","TPR","FPR","MCC","edges"]: d[c]=pd.to_numeric(d[c])
NREP=int(d.groupby(["method","snr"]).size().min())
d=d.groupby(["method","snr"]).head(NREP)          # 保持平衡

ORDER=["EvoTopo-GRN","GENIE3","dynGENIE3_RF","dynGENIE3_ET",
       "SINDY_STLSQ","SINDY_SR3","SINDY_Lasso","nonlinear_ODE"]
COL=dict(zip(ORDER,["#2a78d6","#eb6834","#1baf7a","#eda100",
                    "#e87ba4","#008300","#4a3aa7","#e34948"]))
SURF,INK,INK2,GRID="#fcfcfb","#0b0b0b","#52514e","#e4e3df"
snrs=sorted(d.snr.unique())

fig,axes=plt.subplots(1,3,figsize=(15.5,4.9),facecolor=SURF,constrained_layout=True)
for ax,met,ylab in zip(axes,["MCC","TPR","FPR"],
        ["MCC","TPR (true positive rate)","FPR (false positive rate)"]):
    ax.set_facecolor(SURF)
    if met=="MCC": ax.axhline(0,color="#c9c8c3",lw=1.2,zorder=1)
    for m in ORDER:
        s=d[d.method==m].groupby("snr")[met].mean().reindex(snrs)
        ax.plot(snrs,s.values,color=COL[m],lw=2.0,marker="o",ms=4.5,
                mfc=COL[m],mec=SURF,mew=0.9,zorder=3)
    ax.set_xlabel("SNR",color=INK2,fontsize=11)
    ax.set_title(ylab,color=INK,fontsize=12.5,loc="left",pad=8)
    ax.set_xticks(snrs)
    ax.grid(color=GRID,lw=.7); ax.set_axisbelow(True)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2,labelsize=9.5,length=3,color=GRID)

# 直接标注可区分的几条
for ax,met in zip(axes,["MCC","TPR","FPR"]):
    for m in (["EvoTopo-GRN","GENIE3"] if met=="MCC" else
              ["SINDY_STLSQ","EvoTopo-GRN","GENIE3"] if met=="TPR" else
              ["SINDY_STLSQ","EvoTopo-GRN"]):
        v=d[d.method==m].groupby("snr")[met].mean().reindex(snrs)
        ax.annotate(m,(snrs[-1],v.values[-1]),textcoords="offset points",
                    xytext=(6,0),ha="left",va="center",color=COL[m],fontsize=9.5)
    ax.set_xlim(min(snrs)-0.3, max(snrs)+2.6)

h=[plt.Line2D([],[],color=COL[m],lw=2.4,marker="o",ms=5,mec=SURF,label=m) for m in ORDER]
fig.legend(handles=h,frameon=False,fontsize=10,labelcolor=INK2,ncol=8,
           loc="upper center",bbox_to_anchor=(0.5,1.10))
fig.suptitle(f"20 genes × 100 samples — SNR 1–20, {NREP} replicates per cell",
             color=INK,fontsize=14.5,x=0.006,ha="left",y=1.17)
fig.savefig(os.environ.get("FIG","../results/bench_snr_sweep.png"),dpi=140,facecolor=SURF,bbox_inches="tight")

# 导出汇总表
out=[]
for met in ["MCC","TPR","FPR"]:
    g=d.groupby(["method","snr"])[met].agg(["mean","std"]).reset_index()
    g["metric"]=met; out.append(g)
tab=pd.concat(out).pivot_table(index=["metric","method"],columns="snr",values="mean")
tab=tab.reindex([(m2,m) for m2 in ["MCC","TPR","FPR"] for m in ORDER])
tab.round(4).to_csv(os.environ.get("SUMMARY","../results/bench_summary.csv"))
print(f"NREP={NREP}"); print(tab.round(3).to_string())
