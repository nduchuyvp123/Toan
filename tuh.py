import sys
from pulp import *

sys.stdin = open('input.txt', 'r')
sys.stdout = open('output.txt', 'w')

n = float(input())
m = int(input())
B = []
v = []
for i in range(0, m):
    t = input().split()
    v.append(float(t[0]))
    B.append(int(t[1]))
t = []


def solve(s, sum, v, j, n):
    global t
    if j == len(v) or sum > n:
        return
    if n-sum < min(v):
        t.append(s.copy())
        return
    s.append(v[j])
    solve(s, sum+v[j], v, j, n)
    s.pop()
    solve(s, sum, v, j+1, n)


solve([], 0, v, 0, n)

var = []

for i in range(0, len(t)):
    var.append("x"+str(i+1))

a = []

for i in range(0, m):
    a.append([].copy())
    for j in range(0, len(t)):
        a[i].append(0)

for i in range(0, len(t)):
    for j in range(0, len(t[i])):
        a[v.index(t[i][j])][i] += 1

A = []

for i in range(0, len(a)):
    A.append({})
    for j in range(0, len(var)):
        A[i][var[j]] = a[i][j]

prob = LpProblem("cutting_stock_problem", LpMinimize)

prob_var = LpVariable.dicts("variable", var, 0, cat=LpInteger)

prob += lpSum(prob_var[i] for i in var)

for r in range(0, len(A)):
    prob += lpSum(A[r][i] * prob_var[i] for i in var) >= B[r]

prob.solve()

print("Status: ", LpStatus[prob.status])

# for v in prob.variables():
#     print(v.name, " = ", v.varValue)

for i in range(len(prob.variables())):
    v = prob.variables()[i]
    print(v.name, " = ", v.varValue, '\t', t[i])

print("Optimal value = ", value(prob.objective))

sys.stdout.close()
