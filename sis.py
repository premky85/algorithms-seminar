import random

import matplotlib.pyplot as plt
import networkx as nx

n = 10000

alpha = 0.2  # infection probability
steps = 20  # steps for SIR processs


def avg(l):
    return sum(l) / len(l)


def calc_sis_network_model(G: nx.Graph):
    I = set()
    S = set()

    # select initial infected nodes
    for node in G.nodes:
        if random.random() <= alpha:
            I.add(node)
        else:
            S.add(node)

    It = [len(I)]
    St = [len(S)]
    for t in range(steps):
        _I = set()
        _Is = set()
        for node in I:

            # infect new neighbors from suspectible set S
            for neighbor in G.neighbors(node):
                if random.random() <= alpha and neighbor in S:
                    S.remove(neighbor)
                    _I.add(neighbor)

            _Is.add(node)

        I.update(_I)

        # infected nodes are now recovered and moved to set S
        for node in _Is:
            I.remove(node)
            S.add(node)

        It.append(len(I))
        St.append(len(S))

    lg = len(G.nodes)
    It = [x / lg for x in It]
    St = [x / lg for x in St]

    return It, St, lg


def calc_sis_network_local_multiple(G: nx.Graph, k, q):
    It = [[] for _ in range(steps + 1)]
    St = [[] for _ in range(steps + 1)]

    # run SIS process q times
    for _ in range(q):
        node = random.choice(list(G.nodes))
        _G = nx.ego_graph(G, node, k)

        I = set()
        S = set()

        # select initial infected nodes
        for node in _G.nodes:
            if random.random() <= alpha:
                I.add(node)
            else:
                S.add(node)

        It[0].append(len(I))
        St[0].append(len(S))

        for t in range(1, steps + 1):
            _I = set()
            _Is = set()
            for node in I:

                # infect new neighbors from suspectible set S
                for neighbor in _G.neighbors(node):
                    if random.random() <= alpha and neighbor in S:
                        S.remove(neighbor)
                        _I.add(neighbor)

                _Is.add(node)

            I.update(_I)

            # infected nodes are now recovered and moved to set R
            for node in _Is:
                I.remove(node)
                S.add(node)

            It[t].append(len(I))
            St[t].append(len(S))

    # calculate average values from q runs
    lg = len(_G.nodes)
    It = [avg(l) / lg for l in It]
    St = [avg(l) / lg for l in St]

    return It, St, lg


def plot_sis(graph, k, q):
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

    It1, St1, n1 = calc_sis_network_model(G)

    It2, St2, n2 = calc_sis_network_local_multiple(G, k, q)

    fig, axs = plt.subplots(2)
    fig.tight_layout(pad=3)

    axs[0].title.set_text(f"Full graph n: {n1}")
    axs[0].plot(It1, label="I")
    axs[0].plot(St1, label="S")

    axs[1].title.set_text(f"Local graph n: {n2}")
    axs[1].plot(It2, label="I")
    axs[1].plot(St2, label="S")

    axs[0].legend()
    axs[1].legend()
    plt.show()


if __name__ == "__main__":
    plot_sis("preferential_attachment", 5, 20)
    plot_sis("configuration_model", 5, 20)