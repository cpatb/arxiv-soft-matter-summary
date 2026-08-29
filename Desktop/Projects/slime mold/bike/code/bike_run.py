"""
bike_run.py -- grow the ideal low-stress (protected) bike network for Chicago and
score it against the REAL protected-lane network to find the highest-value gaps.
"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import networkx as nx

import physarum as P
import data_bike as db
P.dc = db                      # reuse CandidateGraph with the bike instance
import data_chicago            # noqa (kept for parity)

# Write figures/data into the repo (matches run_experiments.py convention) so
# outputs land where README points; falls back to CWD-relative if __file__ absent.
_HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(_HERE, "..", "figures"); DAT = os.path.join(_HERE, "..", "data")
os.makedirs(FIG, exist_ok=True); os.makedirs(DAT, exist_ok=True)
LAKE="#bcd6e8"; LAND="#f4f2ec"; C_PROT="#1d6f42"; C_GAP="#c1121f"; C_IDEAL="#7a4eab"
plt.rcParams.update({"font.size":9,"figure.dpi":140})
SPACING=1.1; GAMMA=1.2


def lonlat(cg,n):
    x,y=cg.xy[n]; return (db.LON0+x/db._KM_PER_DEG_LON, db.LAT0+y/db._KM_PER_DEG_LAT)

def basemap(ax):
    poly=np.array(db.LAND_POLYGON); ax.set_facecolor(LAKE)
    ax.fill(poly[:,0],poly[:,1],color=LAND,zorder=0,ec="#9bb7c9",lw=0.8)
    ax.set_xlim(-87.97,-87.5); ax.set_ylim(41.62,42.10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect(1/np.cos(np.radians(41.85)))

def draw_protected(ax,lw=2.4,color=C_PROT,alpha=0.9,label=True):
    for name,(pts,mi) in db.PROTECTED_CORRIDORS.items():
        xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
        ax.plot(xs,ys,color=color,lw=lw,alpha=alpha,solid_capstyle="round",zorder=3)

def dist_to_protected(mx,my):
    best=1e9
    for name,(pts,mi) in db.PROTECTED_CORRIDORS.items():
        for a,b in zip(pts[:-1],pts[1:]):
            ax_,ay=db.project(*a); bx,by=db.project(*b)
            best=min(best,P._point_to_segment_km(mx,my,float(ax_),float(ay),float(bx),float(by)))
    return best


def main():
    cg=P.CandidateGraph(spacing_km=SPACING, weight_coeffs=(0.35,0.30,0.35))
    print("nodes",cg.N,"edges",cg.E,"terminals",len(cg.terminals))
    res=P.physarum_solve(cg,gamma=GAMMA,steps=350,batch=10,I0_scale=6.0,seed=4,record_every=0)
    G=repair(cg,P.physarum_network(cg,res['D'],frac=0.12)[0])
    m=P.network_metrics(G,cg)
    print("ideal net: edges=%d TL=%.0fkm FT=%.2f MD=%.1f"%(m['n_edges'],m['TL'],m['FT'],m['MD']))

    # classify edges: aligned with existing protected corridor vs gap
    THR=0.45
    for u,v,d in G.edges(data=True):
        mx=0.5*(cg.xy[u,0]+cg.xy[v,0]); my=0.5*(cg.xy[u,1]+cg.xy[v,1])
        d['gap']= dist_to_protected(mx,my) > THR
    aligned_km=sum(d['length'] for _,_,d in G.edges(data=True) if not d['gap'])
    gap_km=sum(d['length'] for _,_,d in G.edges(data=True) if d['gap'])
    print("aligned %.1f km, gap %.1f km"%(aligned_km,gap_km))

    # rank gap corridors: connected components of gap edges
    Gg=nx.Graph()
    for u,v,d in G.edges(data=True):
        if d['gap']: Gg.add_edge(u,v,**d)
    # equity/demand lookup by nearest demand node
    dem=[(k,db.project(*a['lonlat']),db.demand_weight(a),a['underserved']) for k,a in db.DEMAND_NODES.items()]
    def near_demand(nid):
        x,y=cg.xy[nid]; best=None;bd=1e9
        for k,(dx,dy),w,u in dem:
            dd=np.hypot(x-float(dx),y-float(dy))
            if dd<bd: bd=dd;best=(k,w,u)
        return best,bd
    comps=[c for c in nx.connected_components(Gg) if sum(Gg.degree(n) for n in c)>=2]
    rows=[]
    for c in nx.connected_components(Gg):
        sub=Gg.subgraph(c)
        L=sum(d['length'] for *_,d in sub.edges(data=True))
        if L<1.2: continue
        flux=sum(d['D']*d['length'] for *_,d in sub.edges(data=True))
        names={}; eq=[]
        for n in c:
            (k,w,u),bd=near_demand(n)
            if bd<1.6:
                names[k]=names.get(k,0)+1; eq.append(u)
        equity=float(np.mean(eq)) if eq else 0.3
        label=", ".join(sorted(names,key=names.get,reverse=True)[:2]) or "corridor"
        score=flux*(0.6+0.8*equity)
        rows.append(dict(label=label,length_km=round(L,1),flux=round(flux,1),
                         equity=round(equity,2),score=round(score,1),
                         nodes=list(c)))
    rows.sort(key=lambda r:r['score'],reverse=True)
    top=rows[:10]
    print("top gap corridors:")
    for i,r in enumerate(top,1):
        print(f"  {i}. {r['label']:30s} {r['length_km']:5.1f}km equity={r['equity']:.2f} score={r['score']:.1f}")

    # ---- figures ----
    figA(cg); figB(cg,G,res); figC(cg,G); figD(cg,G,top)

    out=dict(facility_miles=db.FACILITY_MILES,
             ideal_metrics={k:(round(v,3) if isinstance(v,float) else v) for k,v in m.items()},
             aligned_km=round(aligned_km,1),gap_km=round(gap_km,1),
             threshold_km=THR, gamma=GAMMA, spacing_km=SPACING,
             top_gaps=[{k:r[k] for k in ['label','length_km','flux','equity','score']} for r in top])
    json.dump(out,open(os.path.join(DAT,"bike_results.json"),"w"),indent=2)
    print("saved results")
    return top


def figA(cg):
    fig,(ax,axb)=plt.subplots(1,2,figsize=(12,6.6),gridspec_kw={'width_ratios':[1.5,1]})
    basemap(ax); draw_protected(ax)
    for k,a in db.DEMAND_NODES.items():
        lo,la=a['lonlat']
        col={'p':'#1f5fa8','u':'#e08a00'}['u' if a['underserved']>=0.7 else 'p']
        ax.scatter([lo],[la],s=14+50*db.demand_weight(a),c=col,ec='white',lw=0.5,zorder=5)
    ax.set_title("Chicago's existing PROTECTED bike lanes (green)\n"
                 "real CDOT data: only 68.7 protected mi citywide",fontsize=9)
    ax.legend(handles=[Line2D([0],[0],color=C_PROT,lw=3,label='protected bike lane (real)'),
                       Line2D([0],[0],marker='o',color='w',markerfacecolor='#e08a00',markersize=8,label='high-equity demand'),
                       Line2D([0],[0],marker='o',color='w',markerfacecolor='#1f5fa8',markersize=8,label='other demand')],
              loc='lower left',fontsize=7)
    # bar chart of facility mileage
    types=list(db.FACILITY_MILES.keys()); vals=[db.FACILITY_MILES[t] for t in types]
    cols=['#1d6f42' if 'Protected' in t else ('#7fb069' if 'Greenway' in t else '#c9a227') for t in types]
    axb.barh(range(len(types)),vals,color=cols)
    axb.set_yticks(range(len(types))); axb.set_yticklabels([t.replace(' Bike Lane','') for t in types],fontsize=8)
    axb.invert_yaxis(); axb.set_xlabel("centerline miles")
    for i,v in enumerate(vals): axb.text(v+1,i,f"{v:.0f}",va='center',fontsize=8)
    axb.set_title("Facility mileage by type\n(protected = 15% of total)",fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"figB1_existing_protected.png")); plt.close(fig)

def figB(cg,G,res):
    fig,ax=plt.subplots(figsize=(7,7.4)); basemap(ax)
    segs=[[lonlat(cg,u),lonlat(cg,v)] for u,v in G.edges()]
    ws=np.array([d['D'] for *_,d in G.edges(data=True)])
    lc=LineCollection(segs,colors=C_IDEAL,linewidths=0.6+2.6*ws/ws.max(),alpha=0.9,zorder=3)
    ax.add_collection(lc)
    for k,a in db.DEMAND_NODES.items():
        lo,la=a['lonlat']; ax.scatter([lo],[la],s=10,c='k',zorder=5)
    ax.set_title("Ideal connected low-stress bike network\n(slime-mold model, demand+equity weighted)",fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"figB2_ideal_network.png")); plt.close(fig)

def figC(cg,G):
    fig,ax=plt.subplots(figsize=(7,7.4)); basemap(ax)
    draw_protected(ax,lw=1.4,color="#9c9c9c",alpha=0.8)
    seg_al=[[lonlat(cg,u),lonlat(cg,v)] for u,v,d in G.edges(data=True) if not d['gap']]
    seg_gap=[[lonlat(cg,u),lonlat(cg,v)] for u,v,d in G.edges(data=True) if d['gap']]
    ax.add_collection(LineCollection(seg_al,colors=C_PROT,linewidths=2.0,alpha=0.85,zorder=3))
    ax.add_collection(LineCollection(seg_gap,colors=C_GAP,linewidths=2.4,alpha=0.9,zorder=4))
    for k,a in db.DEMAND_NODES.items():
        lo,la=a['lonlat']
        if a['underserved']>=0.8: ax.scatter([lo],[la],s=30,c='#e08a00',ec='k',lw=0.4,zorder=6)
    ax.set_title("Gap analysis: ideal network vs. reality\n"
                 "green = already protected   red = PRIORITY protected-lane gap",fontsize=9)
    ax.legend(handles=[Line2D([0],[0],color=C_PROT,lw=3,label='ideal corridor already protected'),
                       Line2D([0],[0],color=C_GAP,lw=3,label='priority gap (build/upgrade)'),
                       Line2D([0],[0],color='#9c9c9c',lw=2,label='existing protected (real)'),
                       Line2D([0],[0],marker='o',color='w',markerfacecolor='#e08a00',markersize=8,label='high-equity area')],
              loc='lower left',fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"figB3_gap_overlay.png")); plt.close(fig)

def figD(cg,G,top):
    fig,ax=plt.subplots(figsize=(8,5)); 
    labels=[f"{i+1}. {r['label']}" for i,r in enumerate(top)][::-1]
    scores=[r['score'] for r in top][::-1]
    eq=[r['equity'] for r in top][::-1]
    cols=[plt.cm.YlOrRd(0.3+0.6*e) for e in eq]
    ax.barh(range(len(top)),scores,color=cols)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(labels,fontsize=8)
    ax.set_xlabel("priority score  (flux-weighted length x equity)")
    ax.set_title("Top-10 priority protected-lane gap corridors\n(color = equity/crash severity)",fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"figB4_priority_ranking.png")); plt.close(fig)

def repair(cg,G):
    if G.number_of_edges()==0: return G
    main=max(nx.connected_components(G),key=len); Gc=cg.edge_graph(weight="cost")
    for d in cg.demand_terms:
        if d in main: continue
        dist,path=nx.single_source_dijkstra(Gc,d,weight="cost")
        tg=[(dist[n],n) for n in main if n in dist]
        if not tg: continue
        _,t=min(tg); nodes=path[t]
        for a,b in zip(nodes[:-1],nodes[1:]):
            ek=P._edge_lookup(cg,a,b)
            G.add_edge(a,b,length=float(cg.length[ek]),cost=float(cg.cost[ek]),D=1e-3,gap=True)
        main=max(nx.connected_components(G),key=len)
    return G

if __name__=="__main__":
    main()
