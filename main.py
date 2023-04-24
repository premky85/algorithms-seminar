import random

import matplotlib.pyplot as plt
import networkx as nx

n = 10000
# preferential attachment model
# graph = 'preferential_attachment'
# degree = 2
# G = nx.barabasi_albert_graph(n, degree)

# # configuration model
# graph = 'configuration_model'
# z=[int(random.gammavariate(alpha=9.0,beta=1)) for i in range(n)]
# G=nx.configuration_model(z)



alpha = 0.2# infection probability
steps = 10 # steps for SIR process

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
    It = [[] for _ in range(steps + 1)]
    St = [[] for _ in range(steps + 1)]
    Rt = [[] for _ in range(steps + 1)]

    # run SIR process q times
    for _ in range(q):
        node = random.choice(list(G.nodes))
        _G = nx.ego_graph(G, node, k)

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

def calc_local_outbreak_probability(G: nx.Graph, k, q):
    outbreaks = [] # outbreak happened
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
    outbreaks = [] # outbreak happened
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

def plot_sir(graph, k, q):
    if graph == 'preferential_attachment':
        degree = 2
        G = nx.barabasi_albert_graph(n, degree)
    elif graph == 'configuration_model':
        G = None
        while not G:
            try:
                z=[int(random.gammavariate(alpha=9.0,beta=1)) for i in range(n)]
                G=nx.configuration_model(z)
            except Exception:
                continue
    else:
        raise RuntimeError('parameter graph must be one of: preferential_attachment, configuration_model')


    It1, St1, Rt1, n1 = calc_sir_network_model(G)

    It2, St2, Rt2, n2 = calc_sir_network_local_multiple(G, k, q)

    fig, axs= plt.subplots(2)
    fig.tight_layout(pad=3)

    axs[0].title.set_text(f'Full graph n: {n1}')
    axs[0].plot(It1, label='I')
    axs[0].plot(St1, label='S')
    axs[0].plot(Rt1, label='R')

    axs[1].title.set_text(f'Local graph n: {n2}')
    axs[1].plot(It2, label='I')
    axs[1].plot(St2, label='S')
    axs[1].plot(Rt2, label='R')

    axs[0].legend()
    axs[1].legend()
    plt.show()


def plot_outbreak_probability(graph, q):
    if graph == 'gowalla':
        G = nx.read_edgelist('networks/Gowalla_edges.txt', delimiter="\t", create_using=nx.DiGraph(), nodetype=int)
    elif graph == 'brightkite':
        G = nx.read_edgelist('networks/Brightkite_edges.txt', delimiter="\t", create_using=nx.DiGraph(), nodetype=int)
    elif graph == 'haggle':
        G = nx.read_edgelist('networks/haggle_human_proximity.csv', delimiter=",", create_using=nx.DiGraph(), nodetype=int, data=(('weight', float), ('time', float)))
    else:
        raise RuntimeError('parameter graph must be one of: gowalla, brightkite, haggle')

    op1 = []
    op2 = []

    for i in range(10):
        print(i)
        op1.append(calc_outbreak_probability(G, i, q))
        op2.append(calc_local_outbreak_probability(G, i, q))


    plt.title(f'Probability of Outbreak')
    plt.plot(op1, label='Full graph')
    plt.plot(op2, label='Subgraph')


    plt.legend()
    plt.show()

    
if __name__ == '__main__':
    plot_sir('configuration_model', 5, 20)
    plot_outbreak_probability('haggle', 20)

    

    

    


