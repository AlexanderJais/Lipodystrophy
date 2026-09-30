"""STRING protein-interaction network of male DEGs with no leptin effect
(persistent down) versus those normalized by leptin (normalized down).

Edges: STRING v12, Mus musculus, combined score >= 0.7 (high confidence).
Outputs results/network/.
"""
from pathlib import Path

import networkx as nx
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "network"
OUT.mkdir(parents=True, exist_ok=True)
API = "https://string-db.org/api/json"
SCORE = 700
ARTEFACTS = ("Gm", "mt-", "Mir", "Lars2", "Tcaf2")


def lfc(name):
    df = pd.read_csv(ROOT / "data" / f"DEG{name}.csv").drop_duplicates("SYMBOL")
    return df.set_index("SYMBOL").log2FoldChange


def partitions():
    s, l = lfc("male_LDsaline_WT"), lfc("male_LDleptin_WT")
    s = s[~s.index.str.startswith(ARTEFACTS)]
    pers = s.index.intersection(l.index)
    norm = s.index.difference(l.index)
    return {"persistent_down": list(pers[s[pers] < 0]),
            "normalized_down": list(norm[s[norm] < 0])}


def string(endpoint, genes):
    r = requests.post(f"{API}/{endpoint}", data={
        "identifiers": "\r".join(genes), "species": 10090,
        "required_score": SCORE, "caller_identity": "lipodystrophy_analysis"},
        timeout=120)
    r.raise_for_status()
    return r.json()


if __name__ == "__main__":
    s, l = lfc("male_LDsaline_WT"), lfc("male_LDleptin_WT")
    summary = []
    for name, genes in partitions().items():
        edges = pd.DataFrame(string("network", genes))
        edges = edges[["preferredName_A", "preferredName_B", "score"]].drop_duplicates()
        edges.to_csv(OUT / f"{name}_edges.csv", index=False)
        enr = string("ppi_enrichment", genes)[0]
        g = nx.from_pandas_edgelist(edges, "preferredName_A", "preferredName_B")
        g.add_nodes_from(genes)
        comps = sorted(nx.connected_components(g), key=len, reverse=True)
        main = g.subgraph(comps[0])
        nodes = pd.DataFrame({
            "gene": list(g.nodes),
            "degree": [g.degree(n) for n in g.nodes],
            "betweenness": pd.Series(nx.betweenness_centrality(g)).reindex(list(g.nodes)).values,
            "in_main_component": [n in comps[0] for n in g.nodes],
            "lfc_saline_vs_wt": [s.get(n) for n in g.nodes],
            "lfc_leptin_vs_wt": [l.get(n) for n in g.nodes]})
        nodes.sort_values("degree", ascending=False).to_csv(OUT / f"{name}_nodes.csv", index=False)
        summary.append({"set": name, "genes": len(genes),
                        "edges": enr["number_of_edges"],
                        "expected_edges": enr["expected_number_of_edges"],
                        "ppi_enrichment_p": enr["p_value"],
                        "main_component": len(comps[0])})
        print(name, summary[-1])
        print("  hubs:", ", ".join(nodes.sort_values("degree", ascending=False)
                                   .head(15).apply(lambda r: f"{r.gene}({r.degree})", axis=1)))
    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
