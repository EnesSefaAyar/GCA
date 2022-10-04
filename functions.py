import functions as fun
import pandas as pd
import numpy as np
import networkx as nx
from copy import deepcopy
from scipy import stats

"""Read data into memory"""
time_series = pd.read_csv(r"data\GSE84095.csv", low_memory=False)
reference = pd.read_csv(r"data\Custom_Reference.csv", low_memory=False)
expressions = pd.read_csv(r"data\CCLE_expression.csv", low_memory=False)
interaction_types = pd.read_csv(r"data\BIOGRID_times_series.csv", low_memory=False)
CRISPR = pd.read_csv(r"data\CRISPR_gene_effect.csv", low_memory=False, index_col="DepMap_ID")
keggs = pd.read_csv(r'data\KEGG_Pathway.csv', index_col="Gene", low_memory=False)

# Parameters
essential_thr = -0.75
r_thr = 0.3
p_thr = 0.05
base_exp = 0.58
min_data_point = 300
number_of_iter = 100
step_size = 350


def create_reference():
    Network = nx.DiGraph()
    edges = set([])
    for i, row in reference.iterrows():
        if row[0] in time_series.columns and row[1] in time_series.columns and row[1] in expressions.columns and row[0] in expressions.columns:
            forward_edge = (row[0], row[1], float(row[2]))
            reverse_edge = (row[1], row[0], float(row[2]))
            edges.add(forward_edge)
            edges.add(reverse_edge)
    Network.add_weighted_edges_from(edges)  # Weights are confidence values.
    # Remove isolated nodes and sel-edges from Graph
    Network.remove_nodes_from(list(nx.isolates(Network)))
    Network.remove_edges_from(list(nx.selfloop_edges(Network)))
    # Initialize attributes as None
    nx.set_edge_attributes(Network, None, "co_expression")
    nx.set_edge_attributes(Network, None, "interaction_type")
    nx.set_node_attributes(Network, None, "essential")
    nx.set_node_attributes(Network, None, "KEGG_Pathways")
    nx.set_node_attributes(Network, 0, "scope")  # Marks desired nodes for cytoscape (binary)
    # Set essential genes
    data = CRISPR.to_numpy()
    columns = np.array(CRISPR.columns)
    c_array = np.median(data, axis=0)  # find essential in all samples.
    nx.set_node_attributes(Network, dict.fromkeys(columns[c_array <= essential_thr], True), "essential")
    # Set biogrid edge attributes (BIOGRID: Direct Interaction, Physical Association, Association.)
    A = interaction_types["Alt IDs Interactor A"]
    B = interaction_types["Alt IDs Interactor B"]
    T = interaction_types["Interaction Types"]
    BIOGRID_edges = []
    for k in range(len(A)):
        if (A[k], B[k]) in Network.edges:
            if (A[k], B[k]) not in BIOGRID_edges:
                BIOGRID_edges.append((A[k], B[k]))
                BIOGRID_edges.append((B[k], A[k]))
                nx.set_edge_attributes(Network, {(A[k], B[k]): T[k]}, "interaction_type")
                nx.set_edge_attributes(Network, {(B[k], A[k]): T[k]}, "interaction_type")
    # Calculate and set co_expressions
    fun.coexpression(Network, r_thr, p_thr, base_exp, min_data_point)
    print("Reference Interactome:")
    print("# of nodes:", len(Network.nodes()), "\n# of edges:", int(len(Network.edges())/2))
    return Network


def construct(Gr, cell_line):
    G1 = deepcopy(Gr)
    G1.name = cell_line
    return G1


def add_attributes(Ga):
    if Ga.name in list(time_series["Gene"]):
        i_exp = list(time_series["Gene"] == Ga.name).index(True)
        for vertice in Ga.nodes:
            # expression
            if vertice in time_series.columns:
                nx.set_node_attributes(Ga, {vertice: time_series.loc[i_exp].at[vertice]}, "expression")
            else:
                nx.set_node_attributes(Ga, {vertice: False}, "expression")
            #  Kegg Pathways
            if vertice in list(keggs.index):
                nx.set_node_attributes(Ga, {vertice: str(list(keggs.loc[vertice])[0]).split(",")[:-1]}, "KEGG_Pathways")
            else:
                nx.set_node_attributes(Ga, {vertice: None}, "KEGG_Pathways")


def coexpression(Ref, r_th=0.6, p_th=0.01, base_expression=0.58, min_data_points=10):
    counter = 0
    for edge in list(Ref.edges):
        x = np.array(expressions[edge[0]])
        y = np.array(expressions[edge[1]])
        index = np.array([False] * len(x))
        index += (x > base_expression) & (y > base_expression)
        x_temp = deepcopy(x)
        x_temp = x_temp[index]
        if len(x_temp) > min_data_points:
            r, p = stats.pearsonr(x, y)
            if np.abs(r) > r_th and p < p_th:
                m, n = np.polyfit(x, y, 1)
                nx.set_edge_attributes(Ref, {edge: (r, m, n)}, "co_expression")
                counter += 1


"""
def mutation_curator(primary, secondary):
    networks = [deepcopy(primary), deepcopy(secondary)]
    # Shrink the maf data
    index = np.array([False] * len(mutations))
    for cell_line in [networks[0].name, networks[1].name]:
        index += list(mutations["DepMap_ID"] == cell_line)
    BRCA_mutations = mutations[index]
    
    # Different mutations to be added to primary (Mutations in secondary but not in primary)
    mut_genes_p = []
    mut_type_p = []
    mut_exp_p = []
    for row in BRCA_mutations.iterrows():
        if row[1]["Hugo_Symbol"] in list(networks[0].nodes) and row[1]["Hugo_Symbol"] in list(networks[1].nodes):
            if networks[1].nodes[row[1]["Hugo_Symbol"]]["mutated"] is not None:
                if networks[0].nodes[row[1]["Hugo_Symbol"]]["mutated"] != networks[1].nodes[row[1]["Hugo_Symbol"]]["mutated"]:
                    if networks[0].nodes[row[1]["Hugo_Symbol"]]["mutated"] is not None:
                        mut_genes_p.append(row[1]["Hugo_Symbol"])
                        mutation_list = []
                        mut_exp_p.append(networks[1].nodes[row[1]["Hugo_Symbol"]]["expression"])
                        for g_change, value in networks[1].nodes[row[1]["Hugo_Symbol"]]["mutated"].items():
                            if g_change not in list(networks[0].nodes[row[1]["Hugo_Symbol"]]["mutated"].keys()):
                                mutation_list.append(g_change)
                        mut_type_p.append(mutation_list)
                    elif networks[0].nodes[row[1]["Hugo_Symbol"]]["mutated"] is None:
                        mut_genes_p.append(row[1]["Hugo_Symbol"])
                        mutation_list = []
                        mut_exp_p.append(networks[1].nodes[row[1]["Hugo_Symbol"]]["expression"])
                        for g_change, value in networks[1].nodes[row[1]["Hugo_Symbol"]]["mutated"].items():
                            mutation_list.append(g_change)
                        mut_type_p.append(mutation_list)
    return mut_genes_p, mut_type_p, mut_exp_p
"""
