"""EvoTopo-GRN 基准：8 方法 x SNR 区间 x 重复次数（用 SNR_LO / SNR_HI / NREP 环境变量指定）
   设计：外层循环重复次数、内层循环 snr —— 每完成一轮就追加写盘，
   任何时刻 results.csv 都是一份"每个格子重复数相同"的平衡表。"""
import numpy as np, pandas as pd, warnings, time, os, sys; warnings.filterwarnings("ignore")
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(_here, "..", "simulation"), os.path.join(_here, "..", "methods")]
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sim_power_EI import sim_power_EI
from evotopo_grn import power_equation_fit, detrend, select_edges as select

SRC=os.environ.get("MTL_SRC", "../Compare_with_pySINDY_dynGenie3_nonlinearODE.py")  # 用环境变量 MTL_SRC 指定路径
raw=open(SRC).read().replace("X,X_obs,t,true_adj = sim20(snr)","")
ns={"__name__":"cmp"}; exec(compile(raw,SRC,"exec"),ns)
run_pySINDY, run_nonliner_ODE = ns["run_pySINDY"], ns["run_nonliner_ODE"]

# 把 run_dynGenie3 改造成可以切换 RF / ET
i0=raw.index("def run_dynGenie3(X,t):"); i1=raw.index("def run_nonliner_ODE(X,t):")
body=raw[i0:i1].replace("def run_dynGenie3(X,t):", "def run_dynGenie3_tm(X,t,tree_method='RF'):")
body=body.replace("dynGENIE3(TS_data, time_points)", "dynGENIE3(TS_data, time_points, tree_method=tree_method)")
ns2=dict(ns); exec(compile(body,"<dyn>","exec"),ns2); run_dynGenie3_tm=ns2["run_dynGenie3_tm"]

P, N, CUT, NTREES = 20, 100, 0.1, 1000
o=~np.eye(P,dtype=bool)
OUT=os.environ.get("RESULTS","../results/bench_results.csv")

def genie3_ss(X, est=RandomForestRegressor, seed=0):
    """GENIE3 稳态模式：x_j ~ x_-j，不用时间轴"""
    p=X.shape[1]; V=np.zeros((p,p))
    for j in range(p):
        idx=[k for k in range(p) if k!=j]
        m=est(n_estimators=NTREES,max_features="sqrt",random_state=seed,n_jobs=-1).fit(X[:,idx],X[:,j])
        fi=m.feature_importances_; s=fi.sum()
        V[j,idx]= fi/s if s>0 else fi
    return V

def metrics(E,A):
    E=np.asarray(E,int)[o]; T=A[o]
    TP=int(((E==1)&(T==1)).sum()); FP=int(((E==1)&(T==0)).sum())
    TN=int(((E==0)&(T==0)).sum()); FN=int(((E==0)&(T==1)).sum())
    den=np.sqrt(float(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    return dict(TP=TP,FP=FP,TN=TN,FN=FN,
                TPR=TP/(TP+FN) if TP+FN else 0.0,
                FPR=FP/(FP+TN) if FP+TN else 0.0,
                MCC=(TP*TN-FP*FN)/den if den else 0.0, edges=int(E.sum()))

def one(snr, seed):
    d=sim_power_EI(p=P,n=N,snr=float(snr),seed=seed)
    Y,EI,A = d["Y"].values, d["EI"], d["adj_main"].values
    fit=power_equation_fit(d["expr"].values, n_grid=20)
    det=detrend(fit)
    s0,s1,s2 = run_pySINDY(Y,EI)
    preds={
        "EvoTopo-GRN":   np.vstack([select(det,c) for c in range(P)]),
        "GENIE3":        np.where(genie3_ss(Y,RandomForestRegressor)>=CUT,1,0),
        "dynGENIE3_RF":  np.where(run_dynGenie3_tm(Y,EI,"RF")>=CUT,1,0),
        "dynGENIE3_ET":  np.where(run_dynGenie3_tm(Y,EI,"ET")>=CUT,1,0),
        "SINDY_STLSQ":   s0, "SINDY_SR3": s1, "SINDY_Lasso": s2,
        "nonlinear_ODE": np.where(run_nonliner_ODE(Y,EI)>=CUT,1,0),
    }
    return [dict(method=m, snr=snr, seed=seed, true_edges=int(A[o].sum()), **metrics(E,A))
            for m,E in preds.items()]

done=set()
if os.path.exists(OUT):
    prev=pd.read_csv(OUT); done={(r.snr,r.seed) for r in prev.itertuples()}
    print(f"续跑：已完成 {len(done)} 个 (snr,seed) 组合",flush=True)

SNR_LO=int(os.environ.get("SNR_LO","1")); SNR_HI=int(os.environ.get("SNR_HI","10"))
NREP=int(os.environ.get("NREP","50"))
print(f"目标: snr {SNR_LO}..{SNR_HI} x {NREP} 次",flush=True)

for rep in range(NREP):
    seed=1000+rep
    t_rep=time.time()
    for snr in range(SNR_LO,SNR_HI+1):
        if (snr,seed) in done: continue
        t0=time.time(); rows=one(snr,seed)
        pd.DataFrame(rows).to_csv(OUT, mode="a", header=not os.path.exists(OUT), index=False)
        print(f"rep {rep+1}/{NREP}  snr={snr:<3} ({time.time()-t0:.0f}s)",flush=True)
    print(f"=== 第 {rep+1}/{NREP} 轮完成，用时 {(time.time()-t_rep)/60:.1f} 分钟 ===",flush=True)
print("ALL DONE",flush=True)
