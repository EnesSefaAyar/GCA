from copy import deepcopy
import networkx as nx
import numpy as np
import functions as fun
import time
from matplotlib import pyplot as plt

# Parameters
essential_thr = -0.75
r_thr = 0.3
p_thr = 0.05
base_exp = 0.58
min_data_point = 300
number_of_iter = 250
step_size = 350

# Main starts here
start_time = time.time()
CL1 = "GSM2227192"  # 2h cell line name
CL2 = "GSM2227201"  # 48h cell line name
static_networks = []  # Stores initial and final states only
Graph = fun.create_reference()
start = fun.construct(Graph, CL1)
fun.add_attributes(start)
static_networks.append(start)
end = fun.construct(Graph, CL2)
fun.add_attributes(end)
static_networks.append(end)
dynamic_networks = [static_networks[0], deepcopy(static_networks[0])]  # Store networks through iterations

# Starter nodes (effected ones), changed nodes determined in the first iteration to be used in next iterations.
deltas = []  # stores altered nodes in first iteration
neighbors = set([])  # stores neighbors of altered nodes
for node in dynamic_networks[-2].nodes:  # First iteration
    total = 0  # temporarily keeps predicted expression level for each node
    weights = 0  # temporarily keeps total weight for each node
    c = 0  # temporarily keeps number of neighbors affecting the node
    for neighbor in nx.all_neighbors(dynamic_networks[-2], node):  # iterate over neighbors of a node
        if (neighbor, node) in dynamic_networks[-2].edges:  # determine incoming edges from neighbors
            if dynamic_networks[-2].edges[(neighbor, node)]["co_expression"] is not None:  # checks if edge has co_expression
                if dynamic_networks[-2].edges[(neighbor, node)]["interaction_type"] is not None:  # checks if edge has been annotated
                    r, m, n = dynamic_networks[-2].edges[(neighbor, node)]["co_expression"]  # parameters on the edge
                    w = dynamic_networks[-2].edges[(neighbor, node)]["weight"]  # weight od the edge
                    total += max(w * abs(r ** 2) * (dynamic_networks[-2].nodes[neighbor]["expression"] * m + n),
                                 0)  # calculate predicted level of the neighbor on the target node
                    weights += abs(r) * w  # sum the weights
                    c += 1  # sum the number of neighbors
                    neighbors.add(neighbor)  # store the neighbors
    if c > 2:  # check number of neighbors affecting the target node
        deltas.append(node)  # append the altered node
        # avoid essential genes to be inactivated
        if dynamic_networks[-2].nodes[node]["essential"] and static_networks[0].nodes[node]["expression"] > base_exp:
            dynamic_networks[-1].nodes[node]["expression"] = max(
                (dynamic_networks[-1].nodes[node]["expression"] * (step_size - 1) + total / weights) / step_size,
                base_exp + 0.1)
        else:
            dynamic_networks[-1].nodes[node]["expression"] = max(
                (dynamic_networks[-1].nodes[node]["expression"] * (step_size - 1) + total / weights) / step_size, 0)

print("Number of nodes to be updated:", len(deltas))
# Remove not altered and not altering nodes to efficiently use memory
all_nodes = set(dynamic_networks[-1].nodes)
removed_nodes = all_nodes.difference(set(deltas).union(neighbors))
dynamic_networks[-1].remove_nodes_from(removed_nodes)
dynamic_networks[-2].remove_nodes_from(removed_nodes)
# Next iterations
for j in range(number_of_iter):
    dynamic_networks.append(deepcopy(dynamic_networks[-1]))
    for delta in deltas:
        total = 0
        counter = 0
        c = 0
        for neighbor in nx.all_neighbors(dynamic_networks[-2], delta):
            if (neighbor, delta) in dynamic_networks[-2].edges:
                if dynamic_networks[-2].edges[(neighbor, delta)]["co_expression"] is not None:
                    if dynamic_networks[-2].edges[(neighbor, delta)]["interaction_type"] is not None:
                        r, m, n = dynamic_networks[-2].edges[(neighbor, delta)]["co_expression"]
                        w = dynamic_networks[-2].edges[(neighbor, delta)]["weight"]
                        total += max(w * abs(r) * (dynamic_networks[-2].nodes[neighbor]["expression"] * m + n), 0)
                        counter += abs(r) * w
                        c += 1
        if c > 2:
            if dynamic_networks[-2].nodes[delta]["essential"] and static_networks[0].nodes[delta]["expression"] > base_exp:
                dynamic_networks[-1].nodes[delta]["expression"] = max(
                    (dynamic_networks[-1].nodes[delta]["expression"] * (step_size - 1) + total / counter) / step_size,
                    base_exp + 0.1)
            else:
                dynamic_networks[-1].nodes[delta]["expression"] = max(
                    (dynamic_networks[-1].nodes[delta]["expression"] * (step_size - 1) + total / counter) / step_size,
                    0)

    print(f"Iteration: {j + 1}")
print("--- %s seconds for all iterations to finish! ---" % (time.time() - start_time))

# Calculate metrics to evaluate performance
recall = []
precision = []
accuracy = []
changed = []
for i, network in enumerate(dynamic_networks):
    true_positive = 0
    false_negative = 0
    true_negative = 0
    false_positive = 0.0001  # Zero division error!
    for i, delta in enumerate(deltas):
        if static_networks[1].nodes[delta]["expression"] < base_exp < static_networks[0].nodes[delta]["expression"]:
            if network.nodes[delta]["expression"] < base_exp:
                true_positive += 1
                changed.append(delta)
            else:
                false_negative += 1
        elif static_networks[1].nodes[delta]["expression"] > base_exp > static_networks[0].nodes[delta]["expression"]:
            if network.nodes[delta]["expression"] > base_exp:
                true_positive += 1
                changed.append(delta)
            else:
                false_negative += 1
        elif static_networks[1].nodes[delta]["expression"] < base_exp > static_networks[0].nodes[delta]["expression"]:
            if network.nodes[delta]["expression"] < base_exp:
                true_negative += 1
            else:
                false_positive += 1
                changed.append(delta)
        elif static_networks[1].nodes[delta]["expression"] > base_exp < static_networks[0].nodes[delta]["expression"]:
            if network.nodes[delta]["expression"] > base_exp:
                true_negative += 1
            else:
                false_positive += 1
                changed.append(delta)
    accuracy.append(
        100 * (true_positive + true_negative) / (true_positive + true_negative + false_negative + false_positive))
    recall.append(100 * true_positive / (true_positive + false_negative))
    precision.append(100 * true_positive / (true_positive + false_positive))

# Pathway Demo
all_paths = set()
for i, vertice in enumerate(changed):
    if static_networks[0].nodes[vertice]["KEGG_Pathways"] is not None:
        for path in static_networks[0].nodes[vertice]["KEGG_Pathways"]:
            all_paths.add(path)
all_paths = list(all_paths)
counts = [0] * len(all_paths)
for i, vertice in enumerate(changed):
    try:
        for path in static_networks[0].nodes[vertice]["KEGG_Pathways"]:
            counts[all_paths.index(path)] += 1
    except:
        pass
most_altered_path = all_paths[counts.index(max(counts))]
print("Most altered path:", most_altered_path)

scopes = set()
for node in changed:
    if static_networks[0].nodes[node]["KEGG_Pathways"] is not None:
        if most_altered_path in static_networks[0].nodes[node]["KEGG_Pathways"]:
            scopes.add(node)

for i, network in enumerate(dynamic_networks):
    for node in scopes:
        network.nodes[node]["scope"] = 1
        for neighbor in nx.all_neighbors(network, node):
            if network.nodes[neighbor]["KEGG_Pathways"] is not None:
                if most_altered_path in network.nodes[neighbor]["KEGG_Pathways"]:
                    network.nodes[neighbor]["scope"] = 1
    #  Setting some attributes to "None" to save networks.
    nx.set_edge_attributes(network, "None", "co_expression")
    nx.set_edge_attributes(network, "None", "interaction_type")
    nx.set_node_attributes(network, "None", "essential")
    nx.set_node_attributes(network, "None", "KEGG_Pathways")
    nx.write_graphml(network, r"results/yolak_demo/iterasyonlar/" + f"Iteration{i}.graphml")

####################################################################

# Plot the metrics
plt.plot(np.linspace(0, len(accuracy) - 1, len(accuracy)), accuracy, c="blue")
plt.text(len(accuracy) - 1, accuracy[-1] + 2, "%" + str(round(accuracy[-1], 2)), c="blue")
plt.plot(np.linspace(0, len(recall) - 1, len(recall)), recall, c="green")
plt.text(len(recall) - 1, recall[-1] + 2, "%" + str(round(recall[-1], 2)), c="green")
plt.plot(np.linspace(0, len(precision) - 1, len(precision)), precision, c="red")
plt.text(len(precision) - 1, precision[-1] - 2, "%" + str(round(precision[-1], 2)), c="red")
plt.legend([f"Doğruluk", f"Duyarlılık", f"Hassasiyet"])
plt.title(f"GCA Performans Analizi")
plt.ylabel("%")
plt.xlabel("İterasyon")
plt.savefig(r"results/GCA_performans_analizi")
