import numpy as np
import networkx as nx



def compute_and_save_all_shortest_paths(G):
    
    """
    Input: nx.graph, G
    Output: list, paths, where paths[start][end] = [ path1, path2 ]
                where path1 = [node1,nod2] is one of (possibly more than one) shortest paths from 
                start to end
    """
    
    nodes = [node for node in G.nodes()]
    n = len(nodes)
    
    paths = list(np.zeros((n,n)))
    paths = [list(i) for i in paths]

    for start in nodes:
        for end in nodes:
            path = nx.all_shortest_paths(G,start,end)
            path = [p for p in path]
            paths[start][end] = path
    return paths



def random_walk_covertime(G, m = 1,num_trials=10):
    """
    Does a random walk until all nodes have been covered
    at least m times. Does this for num_trials times.
    """
    
    nodes = [node for node in G.nodes()]  #assume labeled 0,1,2,....n-1
    Ts = []
    trial = 0
    while trial < num_trials:
        trial += 1
        counts = np.zeros(len(nodes)) #count[i] count of node is
        time = 0
        current_position = np.random.choice(nodes)
        counts[current_position] += 1
        num_unvisited_nodes = np.sum(counts < m)
        while num_unvisited_nodes > 0:
            time += 1
            neighbours = G[current_position].keys()
            new_position = np.random.choice(neighbours)
            counts[new_position] += 1
            current_position = new_position
            num_unvisited_nodes = np.sum(counts < m)
        Ts.append(time)
    return Ts



def urban_explorer_covertime(G,m,num_trials=10):
    """
    Does an urban explorer until all nodes have been covered
    at least m times. Does this for num_trials times.
    """

    nodes = [node for node in G.nodes()]
    Ts = []
    trial = 0
    
    paths = compute_and_save_all_shortest_paths(G)
    
    while trial < num_trials:
        trial += 1
        time = 0
        counts = np.zeros(len(nodes)) #count[i] count of node is        
        current_position = np.random.choice(nodes)
        counts[current_position] += 1
        num_unvisited_nodes = np.sum(counts < m)
        while num_unvisited_nodes > 0:
            nodes_minus_origin = nodes[:current_position] + nodes[current_position+1:]
            destination = np.random.choice(nodes_minus_origin)
            
            #Pick one of possibly many shortest paths
            all_shortest_paths = paths[current_position][destination]
            temp = np.random.choice(range(len(all_shortest_paths)))
            path = all_shortest_paths[temp]
            path = path[1:] #remove the origin
            
            #Traverse path
            for node in path:
                counts[node] += 1
                num_unvisited_nodes = np.sum(counts < m)
                time += 1
            current_position = destination
        Ts.append(time)
    return Ts



def urban_explorer_stationary_densities(G,counts,T):
    """
    Does an urban explorer process on the nx.graph G for T timesteps.
    Counts[i] = number of times nodes i has been touched. I'll feed this in
    as an arrays of zeros to start. But I can run simulations back-to-back
    if I need to. 
    """

    nodes = [node for node in G.nodes()]
    Ts = []
    trial = 0
    
    paths = compute_and_save_all_shortest_paths(G)
    
    trial += 1
    time = 0
    current_position = np.random.choice(nodes)
    counts[current_position] += 1
    while time < T:
        nodes_minus_origin = nodes[:current_position] + nodes[current_position+1:]
        destination = np.random.choice(nodes_minus_origin)

        #Pick one of possibly many shortest paths
        all_shortest_paths = paths[current_position][destination]
        temp = np.random.choice(range(len(all_shortest_paths)))
        path = all_shortest_paths[temp]
        path = path[1:] #remove the origin

        #Traverse path
        for node in path:
            counts[node] += 1
            time += 1
        current_position = destination
    return counts



def urban_explorer_returntime(G,start_node,num_trials=10):
    """
    Does an urban explorer until all nodes have been covered
    at least m times. Does this for num_trials times.
    """

    nodes = [node for node in G.nodes()]
    Ts = []
    trial = 0
    
    paths = f.compute_and_save_all_shortest_paths(G)
    
    while trial < num_trials:
        trial += 1
        time = 0
        counts = np.zeros(len(nodes)) #count[i] count of node is        
        current_position = start_node
        counts[current_position] += 1
        while True:
            cond = False
            nodes_minus_origin = nodes[:current_position] + nodes[current_position+1:]
            destination = np.random.choice(nodes_minus_origin)
            
            #Pick one of possibly many shortest paths
            all_shortest_paths = paths[current_position][destination]
            temp = np.random.choice(range(len(all_shortest_paths)))
            path = all_shortest_paths[temp]
            path = path[1:] #remove the origin
            
            #Traverse path
            for node in path:
                current_position = node
                time += 1
                if current_position == start_node:
                    cond = True
                    break
            if cond == True:
                break
        Ts.append(time)
    return Ts



def random_walk_returntime(G,start_node,num_trials=10):
    """
    Does a random walk until all nodes have been covered
    at least m times. Does this for num_trials times.
    """
    
    nodes = [node for node in G.nodes()]  #assume labeled 0,1,2,....n-1
    Ts = []
    trial = 0
    while trial < num_trials:
        trial += 1
        time = 0
        current_position = start_node
        while True:
            time += 1
            neighbours = G[current_position].keys()
            new_position = np.random.choice(neighbours)
            current_position = new_position
            if current_position == start_node:
                break
        Ts.append(time)
    return Ts