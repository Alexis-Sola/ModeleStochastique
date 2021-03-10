# Match nodes to their levels
levels = [[1], [2, 3, 4], [5, 6, 7], [8, 9, 10]]

# Which node is linked to the others
nodes = {
    1: [2, 3, 4],
    2: [1, 5, 6],
    3: [1, 5, 6, 7],
    4: [1, 6, 7],
    5: [2, 3, 8, 9],
    6: [2, 3, 4, 8, 9, 10],
    7: [3, 4, 9, 10],
    8: [5, 6],
    9: [5, 6, 7],
    10: [6, 7]
}

# Weight on each link
links = {
    "1-->2": 1, "2-->1": 1,
    "1-->3": 5, "3-->1": 5,
    "1-->4": -3, "4-->1": -3,
    "2-->5": -2, "5-->2": -2,
    "2-->6": 5, "6-->2": 5,
    "3-->5": 2, "5-->3": 2,
    "3-->6": 7, "6-->3": 7,
    "3-->7": 4, "7-->3": 4,
    "4-->6": 7, "6-->4": 7,
    "4-->7": 3, "7-->4": 3,
    "5-->8": 3, "8-->5": 3,
    "5-->9": 5, "9-->5": 5,
    "6-->8": -2, "8-->6": -2,
    "6-->9": 1, "9-->6": 1,
    "6-->10": 2, "10-->6": 2,
    "7-->9": 0, "9-->7": 0,
    "7-->10": 4, "10-->7": 4,
}

# Returns value of weith between 2 nodes if a link exist
def GetWeight(begin, end):
    arrow = "-->"
    if IsLinkExist(begin, end):
        return links[str(begin) + arrow + str(end)]
    else:
        raise Exception("Link does not exist between these two nodes.")

# Check if link exist between 2 nodes
def IsLinkExist(begin, end):
    isExist = False
    for val in nodes[begin]:
        if val == end:
            isExist = True
            break
    return isExist

# Get links values of a specific node
def ValueLink(node):
    weight = []
    for val in nodes[node]:
        if node > val:
            weight.append(GetWeight(node, val))
    return weight

# Returns all the previous linked nodes of a specific node
def GetPreviousNodes(node):
    new_tab = []
    for val in nodes[node]:
        if node > val:
            new_tab.append(val)
    return new_tab

# Calculate the optimum cost of node
def V(node):
    linked_node = GetPreviousNodes(node)
    w = ValueLink(node)

    if len(linked_node) > 1:
        tab = []
        cpt = 0
        for val in w:
            tab.append(val + V(linked_node[cpt]))
            cpt = cpt + 1
        if len(tab) > 0:
            return max(tab)
    else:
        return max(w)

# Calculate the optimum cost of the graph
def ReverseRecur(level):
    if level >= len(levels) or level <= 0:
        raise Exception("Level doest not exist in this context.")

    nodes = levels[level]
    max_per_nodes = []
    for node in nodes:
        max_per_nodes.append(V(node))
    return max(max_per_nodes)


print(ReverseRecur(3))