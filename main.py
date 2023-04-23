import random

import matplotlib.pyplot as plt
import networkx as nx

# preferential attachment model
G = nx.barabasi_albert_graph(10000, 3)

# # configuration model
# z=[int(random.gammavariate(alpha=9.0,beta=2.0)) for i in range(10000)]
# G=nx.configuration_model(z)


alpha = 0.2 # infection probability
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

    # print(len(_G.nodes))

    
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


                

if __name__ == '__main__':
    It1, St1, Rt1, n1 = calc_sir_network_model(G)

    It2, St2, Rt2, n2 = calc_sir_network_local(G, 4, 20)

    fig, axs= plt.subplots(2)

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
