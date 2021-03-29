import copy
import numpy as np

#Initiliser le plateau
row = 3
col = 4
#Plateau S
S = np.zeros((row, col))
S[0][3] = 1
S[1][3] = -100
forbidden = [1, 1]

V = np.zeros((row, col))
alpha = 0.9
A = ["haut", "bas", "droite", "gauche"]


P = {
    "haut": 0.0,
    "droite": 0.0,
    "gauche": 0.0,
    "bas": 0.0
}

Mirror = {
    "haut": "bas",
    "droite": "gauche",
    "gauche": "droite",
    "bas": "haut"
}

NumberToAction = {
    0: "bas",
    1: "haut",
    2: "gauche",
    3: "droite"
}


def MoveUp(state):
    state_copy = copy.deepcopy(state)
    if state_copy[0] == 0:
        return False
    state_copy[0] = state[0] - 1
    if state_copy != forbidden:
        return True
    return False


def MoveDown(state):
    state_copy = copy.deepcopy(state)
    if state_copy[0] >= len(S) - 1:
        return False
    state_copy[0] = state[0] + 1
    if state_copy != forbidden:
        return True
    return False


def MoveLeft(state):
    state_copy = copy.deepcopy(state)
    if state_copy[1] == 0:
        return False
    state_copy[1] = state[1] - 1
    if state_copy != forbidden:
        return True
    return False

def MoveRight(state):
    state_copy = copy.deepcopy(state)
    if state_copy[1] >= len(S[0]) - 1:
        return False
    state_copy[1] = state[1] + 1
    if state_copy != forbidden:
        return True
    return False

#On vérifie si l'agent est autorisé à bouger
def AuthorizeMove(state, direction):
    if direction == "haut":
        return MoveUp(state)
    elif direction == "bas":
        return MoveDown(state)
    elif direction == "droite":
        return MoveRight(state)
    elif direction == "gauche":
        return MoveLeft(state)
    else:
        raise Exception("Direction {0} does not exist.".format(direction))

#Permet de bouger l'agent
def ConvertActionToState(state, direction):
    new_state = copy.deepcopy(state)
    if direction[1] == "haut" and direction[0] == "move":
        new_state[0] = state[0] - 1
    elif direction[1] == "bas" and direction[0] == "move":
        new_state[0] = state[0] + 1
    elif direction[1] == "droite" and direction[0] == "move":
        new_state[1] = state[1] + 1
    elif direction[1] == "gauche" and direction[0] == "move":
        new_state[1] = state[1] - 1
    elif direction[0] == "rester":
        return state
    return new_state

#Toutes les directions ou l'agent peut aller
def WhichDirections(state):
    tab = []
    for direction in A:
        if AuthorizeMove(state, direction):
            tab.append(["move", direction])
        else:
            tab.append(["rester", direction])
    return tab

#Change les proba en fonction de la direction entrée
def ChangeP(direction):
    P[direction] = 0.8
    for val in A:
        if val != direction:
            P[val] = 0.1
    P[Mirror[direction]] = 0.0

#Calcule la valeur v pour un état
def ComputeValue(state, s, v):
    tab_actions = []
    for action in A:
        tmp = 0
        ChangeP(action)
        for current_action in WhichDirections(state):
            new_state = ConvertActionToState(state, current_action)
            tmp = tmp + P[Mirror[current_action[1]]] * v[new_state[0]][new_state[1]]
        tab_actions.append(s[state[0]][state[1]] + alpha * tmp)

    v_state = max(tab_actions)
    return v_state, tab_actions

#Calcul la valeur v optimale pour tous les états
def ValueIteration(nbIter, s, v):
    v_copy = copy.deepcopy(v)
    mat_actions = []
    for val in range(nbIter):
        for i in range(row):
            for j in range(col):
                if [i, j] == forbidden:
                    continue
                value, tab = ComputeValue([i, j], s, v_copy)
                v_copy[i][j] = value
                if val == nbIter - 1:
                    mat_actions.append(tab)
    return v_copy, mat_actions


def InsertValue(valueTable):
    new_table = []
    cpt = 0
    for val in valueTable:
        if cpt == 5:
            new_table.append([0])
        new_table.append(val)
        cpt = cpt + 1
    return new_table

#Calcule la politique optimale à suivre pour obtenir le reward
def ExtractPolicy(value_table):
    value_table = InsertValue(value_table)
    policy = []
    tmp = []
    cpt = 1
    i = 0
    for val in value_table:
        if cpt % 5 == 0:
            policy.append(tmp)
            tmp = []
            cpt = 1
        tmp.append(np.argmax(np.array(val)))
        cpt = cpt + 1
        i = i + 1
    policy.append(tmp)
    return policy


valueIter, tab = ValueIteration(100, S, V)
bestPolicy = ExtractPolicy(tab)
bestPolicyConverted = []
#Convertit un nombre en une action
for i in range(row):
    tmp = []
    for j in range(col):
         tmp.append(NumberToAction[bestPolicy[i][j]])
    bestPolicyConverted.append(tmp)

bestPolicyConverted[1][1] = "no action"
print(valueIter)
print(np.asarray(bestPolicyConverted))
