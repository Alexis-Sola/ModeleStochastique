import copy
import numpy as np
import random

row = 3
col = 4
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


def ConvertActionToState(state, direction):
    new_state = copy.deepcopy(state)
    if direction[1] == "haut" and direction[0] == "move" and MoveUp(state):
        new_state[0] = state[0] - 1
    elif direction[1] == "bas" and direction[0] == "move" and MoveDown(state):
        new_state[0] = state[0] + 1
    elif direction[1] == "droite" and direction[0] == "move" and MoveRight(state):
        new_state[1] = state[1] + 1
    elif direction[1] == "gauche" and direction[0] == "move" and MoveLeft(state):
        new_state[1] = state[1] - 1
    elif direction[0] == "rester":
        return state
    return new_state


def WhichDirections(state):
    tab = []
    for direction in A:
        if AuthorizeMove(state, direction):
            tab.append(["move", direction])
        else:
            tab.append(["rester", direction])

    return tab


def ChangeP(direction):
    P[direction] = 0.8
    for val in A:
        if val != direction:
            P[val] = 0.1
    P[Mirror[direction]] = 0.0


def Reward(state, s, v):
    tab_actions = []
    for action in A:
        tmp = 0
        ChangeP(action)
        for current_action in WhichDirections(state):
            new_state = ConvertActionToState(state, current_action)
            tmp = tmp + P[Mirror[current_action[1]]] * v[new_state[0]][new_state[1]]
        tab_actions.append(tmp)

    v_state = s[state[0]][state[1]] + alpha * max(tab_actions)
    return v_state


def ValueIteration(nbIter, s, v):
    v_copy = copy.deepcopy(v)
    for val in range(nbIter):
        for i in range(row):
            for j in range(col):
                if [i, j] == forbidden:
                    continue
                v_copy[i][j] = Reward([i, j], s, v_copy)
    return v_copy


def extract_policy(value_table):
    policy = np.zeros((row, col))

    for i in range(row):
        for j in range(col):
            values = []
            for action in A:
                tmp = 0
                ChangeP(action)
                for current_action in WhichDirections([i, j]):
                    new_state = ConvertActionToState([i, j], current_action)
                    tmp = tmp + P[Mirror[current_action[1]]] * (S[new_state[0]][new_state[1]] + alpha * value_table[new_state[0]][new_state[1]])
                values.append(tmp)
                policy[i][j] = np.argmax(np.array(values))
    return policy


bestPolicy = extract_policy(ValueIteration(100, S, V))
bestPolicyConverted = []
for i in range(row):
    tmp = []
    for j in range(col):
         tmp.append(NumberToAction[bestPolicy[i][j]])
    bestPolicyConverted.append(tmp)

print(bestPolicy)
print(np.asarray(bestPolicyConverted))