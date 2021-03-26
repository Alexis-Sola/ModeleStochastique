import copy
import random

S = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, -100], [0, 0, 0, 1]]
V = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

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


def MoveUp(state):
    if state[0] == 0:
        return False
    return True


def MoveDown(state):
    if state[0] >= len(S) - 1:
        return False
    return True


def MoveLeft(state):
    if state[1] == 0:
        return False
    return True


def MoveRight(state):
    if state[1] >= len(S[0]) - 1:
        return False
    return True


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
    if direction[1] == "haut" and direction[0] == "move":
        new_state[0] = state[0] - 1
    elif direction[1] == "bas" and direction[0] == "move":
        new_state[0] = state[0] + 1
    elif direction[1] == "droite" and direction[0] == "move":
        new_state[1] = state[1] + 1
    elif direction[1] == "gauche" and direction[0] == "move":
        new_state[1] = state[1] - 1
    elif direction[0] == "rester":
        new_state = state
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


# def V2(state, s, v):
#     if state == [3, 3] or s[state[0]][state[1]] < 0:
#         return s[state[0]][state[1]], s
#
#     tab_actions = []
#     for action in A:
#         tmp = 0
#         ChangeP(action)
#         for avaible_action in WhichDirections(state):
#             new_state = ConvertActionToState(state, avaible_action)
#
#             if s[new_state[0]][new_state[1]] <= 0:
#                 continue
#
#             tmp = tmp + P[Mirror[avaible_action]] * V(new_state, s, v)[0]
#         tab_actions.append(tmp)
#
#     v_state = s[state[0]][state[1]] + alpha * max(tab_actions)
#     s[state[0]][state[1]] = v_state
#     return v_state, s


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
    for i in range(nbIter):
        for row in range(len(s)):
            for col in range(len(s)):
                v_copy[row][col] = Reward([row, col], s, v_copy)
    return v_copy

print(ValueIteration(100, S, V))