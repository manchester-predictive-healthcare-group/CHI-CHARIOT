# Python version 3.14.3
import functools
import networkx as nx # Version 3.6.1
import pandas as pd # Version 2.3.3
import numpy as np # Version 2.3.5   
import matplotlib.pyplot as plt # Version 3.10.8
from fire import Fire # Version 0.7.1

def format_edge_labels(graph):
    edge_labels = {}
    for k in graph.edges:
        effect_str = "?"
        or_str = ""
        if 'effect' in graph.edges[*k]:
            effect_str = f"{graph.edges[*k]['effect']:.4g}"
        if 'is_ratio' in graph.edges[*k] and graph.edges[*k]['is_ratio']:
            or_str = "*"
        edge_labels[k] = f"{effect_str}{or_str}"
    return edge_labels
    
def is_ratio(effect):
    effect = effect.upper().strip()
    return effect=="OR" or effect=="HR"

def dag_solver(
        dag_edges_path = "chariot_methods_edges.csv",
        evidence_path = "chariot_methods_evidence.csv",
        output_path = "output.csv",
        dag_node_pos_path = None,
        plot_intermediate_graphs = True,
        print_intermediate_information = True,
        cstyle = "arc3,rad=-0.45"
        ):
    """ Expects three arguments: `dag_edges_path`, `evidence_path`, and `output_path`, all of which are CSV files.
    `dag_edges_path` must contain two columns, one names "source" and another names "target", indicating the source and target nodes of each edge
    `evidence_path` must contain at least five columns, "source and "target", as above, as well as "effect_type", "effect", and "measure".
        "effect_type" must take one of two values: "direct" or "total", showing that the evidence is either a total effect or as a direct effect
        "effect" is the value of the evidence, which is expected all to be in consistent units.
        "measure" is a indicator to show either that the effect is a linear effect ("" or "linear"), or that it's an Odds Ratio ("OR") or Hazard Ratio ("HR").
    `output_path` is the filename which the estimated direct effects will be saved.
    A third, optional path might be provided through the `dag_node_pos_path` variable, with three columns "node", "y", "x", indicating the node and its x/y coordinates for plotting.
    The variables `plot_intermediate_graphs` and `print_intermediate_information` may be set to `False` for the program to not produce intermediate outputs.
    
    This is only a simple example implementation of the algorithm, and it expects the unknown effects to be direct effects between the total effect's source and target nodes.
    If there is a path with only a single unknown effect between the total effect's source and target node, but it isn't solvable, this implementation must be modified.
    We assume an exponential link function between the current running effect and any odds-ratio edge.
    """
    gdf = pd.read_csv(dag_edges_path)
    
    G = nx.DiGraph(gdf[["source","target"]].to_numpy()[:,:].tolist())

    if dag_node_pos_path is not None:
        ndf = pd.read_csv(dag_node_pos_path).set_index("node")
        pos = {}
        for n in ndf.index:
            pos[n] = tuple(map(float,ndf.loc[n].values))
    else:
        pos = nx.drawing.layout.spring_layout(G)

    nx.draw_networkx(G, pos, connectionstyle=cstyle)
    x, xx = plt.xlim()
    y, yy = plt.ylim()
    full_xlim = x-1,xx+1
    full_ylim = y-1,yy+1
    plt.xlim(*full_xlim)
    plt.ylim(*full_ylim)
    plt.savefig("dag_0_0raw.png")
    if False:
        plt.show()
    plt.close()

    edf = pd.read_csv(evidence_path)
    
    direct_evidence = edf[edf["effect_type"].str.lower().str.strip()=="direct"]
    for i in direct_evidence.index:
        source, target, effect, measure = direct_evidence.loc[i,["source","target","effect","measure"]]
        if print_intermediate_information: print(source, target, effect)
        G.edges[source,target]["effect"] = float(effect)
        measure = "" if pd.isna(measure) else measure
        G.edges[source,target]["is_ratio"] = measure.strip().upper()=="OR" or measure.strip().upper()=="HR"

    if plot_intermediate_graphs:
        nx.draw_networkx(G,pos, connectionstyle=cstyle)
        nx.draw_networkx_edge_labels(G,pos,format_edge_labels(G), connectionstyle=cstyle)
        plt.xlim(*full_xlim)
        plt.ylim(*full_ylim)
        plt.savefig("dag_0_1dir.png")
        plt.close()
    
    #%%
    indirect_evidence = edf[edf["effect_type"].str.lower().str.strip()=="total"]
        
    #%%
    j = 0
    k = 0
    any_solvable = True
    while any_solvable:
        j += 1
        any_solvable = False
        for i in indirect_evidence.index:
            source, target, effect, measure = indirect_evidence.loc[i,["source","target","effect","measure"]]
            measure = "" if pd.isna(measure) else measure
            if print_intermediate_information: print(source, target, effect)
            paths = list(nx.all_simple_paths(G, source, target))
            n_missing = [sum("effect" not in G.edges[s,t] for s,t in zip(p[:-1],p[1:])) for p in paths]
            if len(paths)>=1 and sum([n==0 for n in n_missing])==len(paths)-1: # Might be better [n==1 for n in n_missing]==1 ?
                if print_intermediate_information: print("Solvable")
                k += 1
                indirect_effects = []
                for p, n_miss in zip(paths,n_missing):
                    if n_miss>0:
                        continue
                    running_effect = 1
                    for s,t in zip(p[:-1],p[1:]):
                        if "effect" in G.edges[s,t]:
                            if G.edges[s,t]["is_ratio"]:
                                running_effect = G.edges[s,t]["effect"]**(running_effect)
                            else:
                                running_effect *= G.edges[s,t]["effect"]
                    if print_intermediate_information: 
                        print(f"Path {p}")
                        print(f"Direct effect {running_effect}")
                    indirect_effects.append(running_effect)

                is_ratio = measure.strip().upper()=="OR" or measure.strip().upper()=="HR"
                ieff = (np.prod if is_ratio else np.sum)(indirect_effects)
                deff = float(effect/ieff if is_ratio else effect-ieff)
    
                if plot_intermediate_graphs:
                    plt.close()
                    sgraph = G.subgraph(list(functools.reduce(set.union,paths,set())))
                    nx.draw_networkx(sgraph,pos, connectionstyle=cstyle)
                    nx.draw_networkx_edge_labels(sgraph, pos, format_edge_labels(sgraph), connectionstyle=cstyle)
                    plt.title(f"Solvable path from {source} to {target} with TE={effect:.4g} and $\\sum$ IE={ieff:.4g}{'*' if is_ratio else ''}\n DE={deff:.4g}{'*' if is_ratio else ''}")
                    plt.ylim(*full_ylim)
                    plt.xlim(*full_xlim)
                    plt.savefig(f"dag_{j}_{k}_{i}.png")
                    plt.close()
    
                for p, n_miss in zip(paths,n_missing):
                    if n_miss==0:
                        continue
                    if len(p)>2:
                        raise NotImplementedError("Haven't implemented how to solve when we know the total effect but the missing effect is not a direct effect")
    
                    G.edges[source,target]["is_ratio"] = is_ratio
                    if is_ratio:
                        G.edges[source,target]["effect"] = float(effect/np.prod(indirect_effects))
                    else:
                        G.edges[source,target]["effect"] = float(effect-np.sum(indirect_effects))
                    G.edges[source,target]["is_ratio"] = is_ratio
                    if print_intermediate_information: print(f"Direct effect of {source} in {target} is {G.edges[source,target]['effect']}")

                any_solvable = True

    #%%
    if plot_intermediate_graphs:
        nx.draw_networkx(G,pos, connectionstyle=cstyle)
        nx.draw_networkx_edge_labels(G,pos, format_edge_labels(G), connectionstyle=cstyle)
        plt.xlim(*full_xlim)
        plt.ylim(*full_ylim)
        plt.savefig("dag_0_2end.png")
        plt.close()

    d = {
        "source": [],
        "target": [],
        "effect": [],
        "type": [],
    }
    for s,t in G.edges:
        d["source"].append(s)
        d["target"].append(t)
        d["effect"].append(G.edges[s,t]["effect"] if "effect" in G.edges[s,t] else np.nan)
        d["type"].append("HR" if "is_ratio" in G.edges[s,t] else "")
    
    pd.DataFrame(d).to_csv(output_path, index=False)

if __name__=="__main__":
    Fire(dag_solver)
