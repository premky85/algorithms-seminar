import random
import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm
from enum import Enum, auto
from typing import *
import numpy as np

class Compartment(Enum):
    """for compartmental modelling of epidemics"""
    SUSCEPTIBLE = auto()
    INFECTED = auto()
    REMOVED = auto()
    EXPOSED = auto()

    @property
    def plot_color(self) -> str:
        match self.name:
            case "SUSCEPTIBLE": return "#07a0d8"
            case "INFECTED": return "#f77811"
            case "REMOVED": return "#918989"
            case "EXPOSED": return "yellow"

    @property
    def plot_label(self) -> str:
        match self.name:
            case "REMOVED": return "recovered / removed"
            case other: return other.lower()


def plot_disease_spread(ax: plt.Axes, title: str,
                        proportions: Dict[Compartment, List[float]]) -> None:
    ax.set_title(title)
    ax.xaxis.set_label_text("time [steps]")
    ax.yaxis.set_label_text("population percentage [%]")

    num_steps = len(next(iter(proportions.values())))
    ax.stackplot(range(num_steps),
                 *([100 * p for p in props]
                   for props in proportions.values()),
                 labels=[comp.plot_label for comp in proportions.keys()],
                 colors=[comp.plot_color for comp in proportions.keys()])

    ax.set_label("compartments")
    ax.legend()


n = 10000


alpha = 0.2
"""infection probability"""

steps = 10
"""steps for SIR process"""


def avg(l):
    return sum(l) / len(l)


def calc_sir_network_model(G: nx.Graph):
    I = set()
    S = set()
    R = set()

    # select initial infected nodes
    for node in G.nodes:
        if random.random() <= alpha:
            I.add(node)
        else:
            S.add(node)

    It = [len(I)]
    St = [len(S)]
    Rt = [0]
    for t in range(steps):
        _I = set()
        _Ir = set()
        for node in I:

            # infect new neighbors from susceptible set S
            for neighbor in G.neighbors(node):
                if random.random() <= alpha and neighbor in S:
                    S.remove(neighbor)
                    _I.add(neighbor)

            _Ir.add(node)

        I.update(_I)

        # infected nodes are now recovered and moved to set R
        for node in _Ir:
            I.remove(node)
            R.add(node)

        It.append(len(I))
        St.append(len(S))
        Rt.append(len(R))

    lg = len(G.nodes)
    It = [x / lg for x in It]
    St = [x / lg for x in St]
    Rt = [x / lg for x in Rt]

    return It, St, Rt, lg


def calc_sir_network_local(G: nx.Graph, k, q):
    node = random.choice(list(G.nodes))
    _G = nx.ego_graph(G, node, k)

    It = [[] for _ in range(steps + 1)]
    St = [[] for _ in range(steps + 1)]
    Rt = [[] for _ in range(steps + 1)]

    # run SIR process q times
    for _ in range(q):
        I = set()
        S = set()
        R = set()

        # select initial infected nodes
        for node in _G.nodes:
            if random.random() <= alpha:
                I.add(node)
            else:
                S.add(node)

        It[0].append(len(I))
        St[0].append(len(S))
        Rt[0].append(0)

        for t in range(1, steps + 1):
            _I = set()
            _Ir = set()
            for node in I:

                # infect new neighbors from suspectible set S
                for neighbor in _G.neighbors(node):
                    if random.random() <= alpha and neighbor in S:
                        S.remove(neighbor)
                        _I.add(neighbor)

                _Ir.add(node)

            I.update(_I)

            # infected nodes are now recovered and moved to set R
            for node in _Ir:
                I.remove(node)
                R.add(node)

            It[t].append(len(I))
            St[t].append(len(S))
            Rt[t].append(len(R))

    # calculate average values from q runs
    lg = len(_G.nodes)
    It = [avg(l) / lg for l in It]
    St = [avg(l) / lg for l in St]
    Rt = [avg(l) / lg for l in Rt]

    return It, St, Rt, lg


def calc_sir_network_local_multiple(G: nx.Graph, k, q):
    # compartment size counters for each timestep, across ALL graphs
    # (later divided by sum of graph sizes)
    # NOTE: this is biased toward larger subgraphs?

    It, St, Rt = (np.zeros(steps + 1) for _ in range(3))
    total_subgraph_sz = 0

    # run SIR process q times
    for _ in tqdm(range(q), desc=f"SIR simulation on multiple {k}-neighborhoods"):
        node = random.choice(list(G.nodes))
        local_G = nx.ego_graph(G, node, k)
        total_subgraph_sz += len(local_G)

        I = set()
        S = set()
        R = set()

        # select initial infected nodes
        for node in local_G.nodes:
            if random.random() <= alpha:
                I.add(node)
            else:
                S.add(node)

        It[0] += len(I)
        St[0] += len(S)

        for t in range(1, steps + 1):
            _I = set()
            _Ir = set()
            for node in I:

                # infect new neighbors from suspectible set S
                for neighbor in local_G.neighbors(node):
                    if random.random() <= alpha and neighbor in S:
                        S.remove(neighbor)
                        _I.add(neighbor)

                _Ir.add(node)

            I.update(_I)

            # infected nodes are now recovered and moved to set R
            for node in _Ir:
                I.remove(node)
                R.add(node)

            It[t] += len(I)
            St[t] += len(S)
            Rt[t] += len(R)

    # calculate average values from q runs
    It /= total_subgraph_sz
    St /= total_subgraph_sz
    Rt /= total_subgraph_sz
    mean_sz = total_subgraph_sz / q

    return It, St, Rt, mean_sz


def calc_local_outbreak_probability(G: nx.Graph, k, q):
    outbreaks = []  # outbreak happened
    # run SIR process q times
    for _ in range(q):
        node = random.choice(list(G.nodes))
        _G = nx.ego_graph(G, node, k)

        I = {node}
        S = set(_G.nodes) - I
        R = set()

        while len(I) > 0:
            _I = set()
            _Ir = set()
            for node in I:

                # infect new neighbors from suspectible set S
                for neighbor in _G.neighbors(node):
                    if random.random() <= alpha and neighbor in S:
                        S.remove(neighbor)
                        _I.add(neighbor)

                _Ir.add(node)

            I.update(_I)

            # infected nodes are now recovered and moved to set R
            for node in _Ir:
                I.remove(node)
                R.add(node)

        outbreaks.append(len(R) > k)

    return avg(outbreaks)


def calc_outbreak_probability(G: nx.Graph, k, q):
    outbreaks = []  # outbreak happened
    # run SIR process q times
    for _ in range(q):
        node = random.choice(list(G.nodes))

        I = {node}
        S = set(G.nodes) - I
        R = set()

        while len(I) > 0:
            _I = set()
            _Ir = set()
            for node in I:

                # infect new neighbors from suspectible set S
                for neighbor in G.neighbors(node):
                    if random.random() <= alpha and neighbor in S:
                        S.remove(neighbor)
                        _I.add(neighbor)

                _Ir.add(node)

            I.update(_I)

            # infected nodes are now recovered and moved to set R
            for node in _Ir:
                I.remove(node)
                R.add(node)

        outbreaks.append(len(R) > k)

    return avg(outbreaks)


def plot_sir(graph: nx.Graph, k: int, q: int, degree=2, multi_neighborhood=False, save_as: str = None):
    if k > steps:
        print(f"WARNING: infection process can't even reach {k}-neighbors in {steps} steps??")

    if graph == "preferential_attachment":
        G = nx.barabasi_albert_graph(n, degree)
    elif graph == "configuration_model":
        G = None
        while not G:
            try:
                z = [int(random.gammavariate(alpha=9.0, beta=1))
                     for _ in range(n)]
                G = nx.configuration_model(z)
            except Exception:
                continue
    else:
        raise RuntimeError(
            "parameter graph must be one of: preferential_attachment, configuration_model"
        )

    It1, St1, Rt1, n1 = calc_sir_network_model(G)

    if multi_neighborhood:
        It2, St2, Rt2, n2 = calc_sir_network_local_multiple(G, k, q)
    else:
        It2, St2, Rt2, n2 = calc_sir_network_local(G, k, q)

    fig, axs = plt.subplots(2)
    fig.tight_layout(pad=3)

    plot_disease_spread(
        axs[0], title=f"Full graph: $n = {n1}$",
        proportions={
            Compartment.INFECTED: It1,
            Compartment.SUSCEPTIBLE: St1,
            Compartment.REMOVED: Rt1,
        })

    local_title = r"Local graphs: $\bar{{n}} = {}$".format(round(n2)) if multi_neighborhood \
                  else f"Local graph $n = {n2}$"

    plot_disease_spread(
        axs[1], title=local_title,
        proportions={
            Compartment.INFECTED: It2,
            Compartment.SUSCEPTIBLE: St2,
            Compartment.REMOVED: Rt2,
        })

    if save_as is None:
        plt.show()
    else:
        plt.savefig(f"plots/{save_as}.pdf")


def plot_outbreak_probability(graph, q):
    if graph == "gowalla":
        G = nx.read_edgelist(
            "networks/Gowalla_edges.txt",
            delimiter="\t",
            create_using=nx.DiGraph(),
            nodetype=int,
        )
    elif graph == "brightkite":
        G = nx.read_edgelist(
            "networks/Brightkite_edges.txt",
            delimiter="\t",
            create_using=nx.DiGraph(),
            nodetype=int,
        )
    elif graph == "haggle":
        G = nx.read_edgelist(
            "networks/haggle_human_proximity.csv",
            delimiter=",",
            create_using=nx.DiGraph(),
            nodetype=int,
            data=(("weight", float), ("time", float)),
        )
    else:
        raise RuntimeError(
            "parameter graph must be one of: gowalla, brightkite, haggle"
        )

    op1 = []
    op2 = []

    for i in range(10):
        print(i)
        op1.append(calc_outbreak_probability(G, i, q))
        op2.append(calc_local_outbreak_probability(G, i, q))

    plt.title(f"Probability of Outbreak")
    plt.plot(op1, label="Full graph")
    plt.plot(op2, label="Subgraph")

    plt.legend()
    plt.show()


def generate_SIR_plots(k=4):
    global alpha

    # fig. 8, 9
    for alpha_i in [.2, .4]:
        alpha = alpha_i

        for deg in [2,3,4]:
            plot_sir("preferential_attachment", k=k, q=20, degree=deg,
                    multi_neighborhood=False, save_as=f"PA_alfa={alpha_i}_deg={deg}")
            
    # fig. 10
    for alpha_i in [.2, .3, .4]:
        alpha = alpha_i

        plot_sir("configuration_model", k=k, q=20,
                 multi_neighborhood=False, save_as=f"CM_sigle_alfa={alpha_i}")
        
    # fig. 11
    for alpha_i in [.2, .3, .4]:
        alpha = alpha_i

        plot_sir("configuration_model", k=k, q=20,
                 multi_neighborhood=True, save_as=f"CM_multi_alfa={alpha_i}")


if __name__ == "__main__":
    generate_SIR_plots()
    # plot_outbreak_probability("haggle", 20)
