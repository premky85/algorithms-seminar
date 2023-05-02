import random

import matplotlib.pyplot as plt
import networkx as nx

n = 10000

alpha = 0.2  # infection probability
steps = 20  # steps for SIR processs


def avg(l):
    return sum(l) / len(l)


def calc_seir_network_model(G: nx.Graph):
    I = set()
    S = set()
    R = set()
    E = set()
    E_times = {}

    # select initial infected nodes
    for node in G.nodes:
        if random.random() <= alpha:
            I.add(node)
        else:
            S.add(node)

    It = [len(I)]
    St = [len(S)]
    Rt = [0]
    Et = [0]
    for t in range(steps):
        _E = set()
        _e_times = {}
        _Ir = set()
        for node in I:

            # expose new neighbors from suspectible set S
            for neighbor in G.neighbors(node):
                if random.random() <= alpha and neighbor in S:
                    S.remove(neighbor)
                    _E.add(neighbor)
                    _e_times[neighbor] = random.gammavariate(1, 1/alpha)

            _Ir.add(node)
        E.update(_E)
        E_times.update(_e_times)

        # infected nodes are now recovered and moved to set R
        for node in _Ir:
            I.remove(node)
            R.add(node)

        for node in E:
            E_times[node] -= 1
            if E_times[node] <= 0:
                I.add(node)
            
        E = {node for node in E if E_times[node] > 0}
        E_times = {node: E_times[node] for node in E_times if E_times[node] > 0}

        It.append(len(I))
        St.append(len(S))
        Rt.append(len(R))
        Et.append(len(E))

    lg = len(G.nodes)
    It = [x / lg for x in It]
    St = [x / lg for x in St]
    Rt = [x / lg for x in Rt]
    Et = [x / lg for x in Et]

    return It, St, Rt, Et, lg


def calc_seir_network_local_multiple(G: nx.Graph, k, q):
    It = [[] for _ in range(steps + 1)]
    St = [[] for _ in range(steps + 1)]
    Rt = [[] for _ in range(steps + 1)]
    Et = [[] for _ in range(steps + 1)]

    # run SEIR process q times
    for _ in range(q):
        node = random.choice(list(G.nodes))
        _G = nx.ego_graph(G, node, k)

        I = set()
        S = set()
        R = set()
        E = set()
        E_times = {}

        # select initial infected nodes
        for node in _G.nodes:
            if random.random() <= alpha:
                I.add(node)
            else:
                S.add(node)

        It[0].append(len(I))
        St[0].append(len(S))
        Rt[0].append(0)
        Et[0].append(0)

        for t in range(1, steps + 1):
            _E = set()
            _e_times = {}
            _Ir = set()
            for node in I:

                # expose new neighbors from suspectible set S
                for neighbor in G.neighbors(node):
                    if random.random() <= alpha and neighbor in S:
                        S.remove(neighbor)
                        _E.add(neighbor)
                        _e_times[neighbor] = random.gammavariate(1, 1/alpha)

                _Ir.add(node)

            E.update(_E)
            E_times.update(_e_times)

            # infected nodes are now recovered and moved to set R
            for node in _Ir:
                I.remove(node)
                R.add(node)

            for node in E:
                E_times[node] -= 1
                if E_times[node] <= 0:
                    I.add(node)

            E = {node for node in E if E_times[node] > 0}
            E_times = {node: E_times[node] for node in E_times if E_times[node] > 0}

            It[t].append(len(I))
            St[t].append(len(S))
            Rt[t].append(len(R))
            Et[t].append(len(E))

    # calculate average values from q runs
    lg = len(_G.nodes)
    It = [avg(l) / lg for l in It]
    St = [avg(l) / lg for l in St]
    Rt = [avg(l) / lg for l in Rt]
    Et = [avg(l) / lg for l in Et]

    return It, St, Rt, Et, lg


def plot_seir(graph, k, q):
    if graph == "preferential_attachment":
        degree = 2
        G = nx.barabasi_albert_graph(n, degree)
    elif graph == "configuration_model":
        G = None
        while not G:
            try:
                z = [int(random.gammavariate(alpha=9.0, beta=1)) for i in range(n)]
                G = nx.configuration_model(z)
            except Exception:
                continue
    else:
        raise RuntimeError(
            "parameter graph must be one of: preferential_attachment, configuration_model"
        )

    It1, St1, Rr1, Et1, n1 = calc_seir_network_model(G)

    It2, St2, Rt2, Et2, n2 = calc_seir_network_local_multiple(G, k, q)

    fig, axs = plt.subplots(2)
    fig.tight_layout(pad=3)

    axs[0].title.set_text(f"Full graph n: {n1}")
    axs[0].plot(It1, label="I")
    axs[0].plot(St1, label="S")
    axs[0].plot(Rr1, label="R")
    axs[0].plot(Et1, label="E")

    axs[1].title.set_text(f"Local graph n: {n2}")
    axs[1].plot(It2, label="I")
    axs[1].plot(St2, label="S")
    axs[1].plot(Rt2, label="R")
    axs[1].plot(Et2, label="E")

    axs[0].legend()
    axs[1].legend()
    plt.show()


if __name__ == "__main__":
    plot_seir("preferential_attachment", 5, 20)
    plot_seir("configuration_model", 5, 20)