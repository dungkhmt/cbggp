#include <bits/stdc++.h>
using namespace std;

typedef std::vector<int> solution;

int k;
solution emptySolution;

struct GraphCDCL
{
    int n;
    vector<vector<int>> adj;
};

struct Literal
{
    int v;
    int c;
    bool negative;
    Literal(int _v = 0, int _c = 0, bool _negative = false) : v(_v), c(_c), negative(_negative) {}
};

struct Clause
{
    vector<Literal> literals;
};

struct Reason
{
    bool isClause;
    int clauseIndex;
    int u;
    int color;
    Reason() : isClause(false), clauseIndex(-1), u(-1), color(-1) {}
    Reason(int clauseIdx) : isClause(true), clauseIndex(clauseIdx), u(-1), color(-1) {}
    Reason(int uu, int ccol) : isClause(false), clauseIndex(-1), u(uu), color(ccol) {}
};

GraphCDCL G;
int numColors;

vector<Clause> implicitConstraints;
vector<int> colorAssigned;
vector<vector<bool>> availableColors;
Clause conflictClause;
vector<Reason> reasonForAssignment;
vector<vector<Reason>> reasonForRemoval;
vector<int> levelOfAssignment;
vector<vector<int>> levelOfRemoval;
vector<bool> isDecisionAssignment;
vector<int> falseCount;
vector<bool> clauseSatisfied;
vector<vector<vector<int>>> clauseIndicesPos;
vector<vector<vector<int>>> clauseIndicesNeg;
int currentLevel;
unsigned long long currentTime;
bool conflictFound;
bool conflictEmptyClause;
int conflictClauseIndex;
int conflictEmptyVertex;

int countFalse(int cIdx)
{
    int countF = 0;
    for (Literal lit : implicitConstraints[cIdx].literals)
    {
        if (!lit.negative)
            countF += (availableColors[lit.v][lit.c] != true);
        else
            countF += (availableColors[lit.v][lit.c] != false);
    }
    return countF;
}

bool isLiteralTrue(const Literal &lit)
{
    if (!lit.negative)
        return (colorAssigned[lit.v] == lit.c);
    else
    {
        if (colorAssigned[lit.v] != -1)
            return (colorAssigned[lit.v] != lit.c);
        else
            return !availableColors[lit.v][lit.c];
    }
}

bool isLiteralFalse(const Literal &lit)
{
    if (!lit.negative)
    {
        if (colorAssigned[lit.v] != -1)
            return (colorAssigned[lit.v] != lit.c);
        else
            return !availableColors[lit.v][lit.c];
    }
    else
        return (colorAssigned[lit.v] == lit.c);
}

int getLiteralLevel(const Literal &lit)
{
    if (!lit.negative)
    {
        if (colorAssigned[lit.v] != -1)
            return levelOfAssignment[lit.v];
        else
            return levelOfRemoval[lit.v][lit.c];
    }
    else
        return levelOfAssignment[lit.v];
}

int addClause(const Clause &cl)
{
    int idx = implicitConstraints.size();
    implicitConstraints.push_back(cl);
    falseCount.push_back(0);
    clauseSatisfied.push_back(false);

    if (implicitConstraints.size() > 100000)
    {
        cerr << "Too many clauses, increase the limit!" << endl;
        exit(1);
    }

    for (auto lit : cl.literals)
    {
        if (!lit.negative)
            clauseIndicesPos[lit.v][lit.c].push_back(idx);
        else
            clauseIndicesNeg[lit.v][lit.c].push_back(idx);

        if (!isLiteralTrue(lit) && isLiteralFalse(lit))
            falseCount[idx]++;

        if (isLiteralTrue(lit))
            clauseSatisfied[idx] = true;
    }
    return idx;
}

void assignColor(int v, int c, int level, bool decision, Reason reason)
{
    // cout << v << " " << c << " " << level << " " << (int)decision << " " << reason.isClause << endl;
    colorAssigned[v] = c;
    levelOfAssignment[v] = level;
    currentTime++;
    isDecisionAssignment[v] = decision;
    if (reason.isClause || reason.u != -1)
        reasonForAssignment[v] = reason;

    for (int c2 = 1; c2 <= numColors; ++c2)
    {
        if (c2 == c)
            continue;
        if (availableColors[v][c2])
        {
            availableColors[v][c2] = false;
            levelOfRemoval[v][c2] = level;
            reasonForRemoval[v][c2] = Reason(v, c);

            for (int clauseIdx : clauseIndicesPos[v][c2])
            {
                if (clauseSatisfied[clauseIdx])
                    continue;
                falseCount[clauseIdx]++;
                if (/*falseCount[clauseIdx]*/ countFalse(clauseIdx) == implicitConstraints[clauseIdx].literals.size())
                {
                    conflictFound = true;
                    conflictEmptyClause = true;
                    conflictClauseIndex = clauseIdx;
                    conflictClause = implicitConstraints[clauseIdx];
                }
            }
        }
    }

    for (int clauseIdx : clauseIndicesPos[v][c])
        clauseSatisfied[clauseIdx] = true;

    for (int clauseIdx : clauseIndicesNeg[v][c])
    {
        if (clauseSatisfied[clauseIdx])
            continue;
        falseCount[clauseIdx]++;
        if (/*falseCount[clauseIdx]*/ countFalse(clauseIdx) == implicitConstraints[clauseIdx].literals.size())
        {
            conflictFound = true;
            conflictEmptyClause = true;
            conflictClauseIndex = clauseIdx;
            conflictClause = implicitConstraints[clauseIdx];
        }
    }

    for (int u : G.adj[v])
    {
        if (availableColors[u][c])
        {
            availableColors[u][c] = false;
            levelOfRemoval[u][c] = currentLevel;
            reasonForRemoval[u][c] = Reason(v, c);
            for (int clauseIdx : clauseIndicesPos[u][c])
            {
                if (clauseSatisfied[clauseIdx])
                    continue;

                if (countFalse(clauseIdx) == implicitConstraints[clauseIdx].literals.size())
                {
                    conflictFound = true;
                    conflictEmptyClause = true;
                    conflictClauseIndex = clauseIdx;
                    conflictClause = implicitConstraints[clauseIdx];
                }
            }
        }
    }
}

void propagateLiteral(const Literal &lit, int level, int clauseIdx = -1)
{
    if (!lit.negative)
    {
        Reason reason;
        if (clauseIdx >= 0)
            reason = Reason(clauseIdx);
        // cout << "pro " << lit.v << " " << lit.c << endl;
        assignColor(lit.v, lit.c, level, false, reason);
    }
    else
    {
        int v = lit.v;
        int c = lit.c;
        if (availableColors[v][c])
        {
            availableColors[v][c] = false;
            levelOfRemoval[v][c] = level;
            Reason reason;
            if (clauseIdx >= 0)
                reason = Reason(clauseIdx);

            reasonForRemoval[v][c] = reason;

            for (int cIdx : clauseIndicesPos[v][c])
            {
                if (clauseSatisfied[cIdx])
                    continue;
                falseCount[cIdx]++;
                if (/*falseCount[cIdx]*/ countFalse(cIdx) == implicitConstraints[cIdx].literals.size())
                {
                    conflictFound = true;
                    conflictEmptyClause = true;
                    conflictClauseIndex = cIdx;
                    conflictClause = implicitConstraints[cIdx];
                    return;
                }
            }

            for (int cIdx : clauseIndicesNeg[v][c])
                clauseSatisfied[cIdx] = true;
        }
    }
}

void unitPropagation(int level)
{
    bool foundUnit = true;

    while (!conflictFound && foundUnit)
    {
        foundUnit = false;

        for (int v = 1; v <= G.n; ++v)
        {
            if (colorAssigned[v] == -1)
            {
                int availableCount = 0;
                int lastColor = -1;
                for (int c = 1; c <= numColors; ++c)
                    if (availableColors[v][c])
                    {
                        availableCount++;
                        lastColor = c;
                        if (availableCount > 1)
                            break;
                    }

                if (availableCount == 0)
                {
                    conflictFound = true;
                    conflictEmptyVertex = v;

                    Clause cc;
                    for (int c = 1; c <= numColors; ++c)
                    {
                        Literal lit(v, c, false);
                        cc.literals.push_back(lit);
                    }
                    conflictClause = cc;
                    return;
                }
                else if (availableCount == 1)
                {
                    int c = lastColor;
                    assignColor(v, c, level, false, Reason());

                    for (int u : G.adj[v])
                    {
                        // cout << "check " << v << " " << c << "   " << u << endl;
                        if (availableColors[u][c])
                        {
                            availableColors[u][c] = false;
                            levelOfRemoval[u][c] = level;
                            reasonForRemoval[u][c] = Reason(v, c);

                            for (int clauseIdx : clauseIndicesPos[u][c])
                            {
                                if (clauseSatisfied[clauseIdx])
                                    continue;
                                falseCount[clauseIdx]++;
                                if (/*falseCount[clauseIdx]*/ countFalse(clauseIdx) == implicitConstraints[clauseIdx].literals.size())
                                {
                                    conflictFound = true;
                                    conflictEmptyClause = true;
                                    conflictClauseIndex = clauseIdx;
                                    conflictClause = implicitConstraints[clauseIdx];
                                    return;
                                }
                                if (/*falseCount[clauseIdx]*/ countFalse(clauseIdx) == implicitConstraints[clauseIdx].literals.size() - 1 && !clauseSatisfied[clauseIdx])
                                {
                                    Literal unitLit;
                                    bool foundLit = false;
                                    for (Literal lit : implicitConstraints[clauseIdx].literals)
                                        if (!isLiteralFalse(lit))
                                        {
                                            unitLit = lit;
                                            foundLit = true;
                                            break;
                                        }

                                    if (foundLit)
                                    {
                                        // cout << "Pre pro " << unitLit.v << " " << unitLit.c << " " << unitLit.negative << endl;
                                        propagateLiteral(unitLit, level);
                                        if (conflictFound)
                                            return;
                                    }
                                }
                            }
                        }
                    }
                    foundUnit = true;
                    break;
                }
            }
        }

        if (conflictFound)
            break;
        if (foundUnit)
            continue;

        for (int ci = 0; ci < implicitConstraints.size(); ++ci)
        {
            if (clauseSatisfied[ci])
                continue;
            int clauseLen = implicitConstraints[ci].literals.size();
            if (/*falseCount[ci]*/ countFalse(ci) == clauseLen)
            {
                conflictFound = true;
                conflictEmptyClause = true;
                conflictClauseIndex = ci;
                conflictClause = implicitConstraints[ci];
                return;
            }
            if (/*falseCount[ci]*/ countFalse(ci) == clauseLen - 1 && !clauseSatisfied[ci])
            {
                Literal unitLit;
                bool foundLit = false;
                for (Literal lit : implicitConstraints[ci].literals)
                    if (!isLiteralFalse(lit))
                    {
                        unitLit = lit;
                        foundLit = true;
                        break;
                    }

                if (foundLit)
                {
                    propagateLiteral(unitLit, level, ci);
                    if (conflictFound)
                        return;
                    foundUnit = true;
                    break;
                }
            }
        }
    }
}

pair<Clause, int> analyzeConflict()
{
    Clause clause = conflictClause;

    while (true)
    {
        int maxLevel = -1;
        int countMaxLevel = 0;
        Literal lastLitAtMax;
        for (Literal lit : clause.literals)
        {
            int lvl = getLiteralLevel(lit);
            if (lvl > maxLevel)
            {
                lastLitAtMax = lit;
                maxLevel = lvl;
                countMaxLevel = 1;
            }
            else if (lvl == maxLevel)
            {
                lastLitAtMax = lit;
                ++countMaxLevel;
            }
        }

        if (countMaxLevel <= 1)
            break;

        Literal li = lastLitAtMax;

        Reason rs;
        if (!li.negative)
            rs = reasonForRemoval[li.v][li.c];
        else
            rs = reasonForAssignment[li.v];

        if (!li.negative && rs.u != -1)
        {
            if (!isDecisionAssignment[rs.u])
            {
                int u = rs.u;
                int assignedColor = rs.color;

                clause.literals.erase(remove_if(clause.literals.begin(), clause.literals.end(),
                                                [&](const Literal &L)
                                                { return L.v == li.v && L.c == li.c && L.negative == li.negative; }),
                                      clause.literals.end());

                for (int c0 = 1; c0 <= numColors; ++c0)
                {
                    if (c0 == assignedColor)
                        continue;
                    Literal newLit(u, c0, false);
                    bool exists = false;
                    for (Literal L2 : clause.literals)
                        if (L2.v == newLit.v && L2.c == newLit.c && L2.negative == newLit.negative)
                        {
                            exists = true;
                            break;
                        }

                    if (!exists)
                        clause.literals.push_back(newLit);
                }
            }
            else
            {
                int u = rs.u;
                int col = rs.color;
                clause.literals.erase(remove_if(clause.literals.begin(), clause.literals.end(),
                                                [&](const Literal &L)
                                                { return L.v == li.v && L.c == li.c && L.negative == li.negative; }),
                                      clause.literals.end());
                Literal newLit(u, col, false);
                bool exists = false;
                for (Literal L2 : clause.literals)
                    if (L2.v == newLit.v && L2.c == newLit.c && L2.negative == newLit.negative)
                    {
                        exists = true;
                        break;
                    }

                if (!exists)
                    clause.literals.push_back(newLit);
            }
        }
        else if (rs.isClause)
        {
            int rClauseIdx = rs.clauseIndex;
            Clause reasonClause = implicitConstraints[rClauseIdx];

            clause.literals.erase(remove_if(clause.literals.begin(), clause.literals.end(),
                                            [&](const Literal &L)
                                            { return L.v == li.v && L.c == li.c && L.negative == li.negative; }),
                                  clause.literals.end());

            for (Literal l2 : reasonClause.literals)
                if (isLiteralFalse(l2))
                {
                    bool exists = false;
                    for (Literal L3 : clause.literals)
                        if (L3.v == l2.v && L3.c == l2.c && L3.negative == l2.negative)
                        {
                            exists = true;
                            break;
                        }

                    if (!exists)
                        clause.literals.push_back(l2);
                }
        }
        else
            break;
    }

    int highestLevel = -1, secondHighest = -1;
    for (Literal lit : clause.literals)
    {
        int lvl = getLiteralLevel(lit);
        if (lvl > secondHighest)
        {
            if (lvl > highestLevel)
            {
                secondHighest = highestLevel;
                highestLevel = lvl;
            }
            else
            {
                secondHighest = lvl;
            }
        }
    }
    if (secondHighest < 0)
        secondHighest = 0;

    return {clause, secondHighest};
}

void backtrack(int targetLevel)
{
    for (int v = 1; v <= G.n; ++v)
    {
        if (levelOfAssignment[v] > targetLevel)
        {
            int col = colorAssigned[v];
            colorAssigned[v] = -1;
            int assignLvl = levelOfAssignment[v];
            levelOfAssignment[v] = -1;
            reasonForAssignment[v] = Reason();

            for (int cIdx : clauseIndicesPos[v][col])
                if (clauseSatisfied[cIdx])
                {
                    bool anyTrue = false;
                    for (Literal L2 : implicitConstraints[cIdx].literals)
                        if (isLiteralTrue(L2))
                        {
                            anyTrue = true;
                            break;
                        }

                    if (!anyTrue)
                        clauseSatisfied[cIdx] = false;
                }

            for (int cIdx : clauseIndicesNeg[v][col])
                if (!clauseSatisfied[cIdx])
                    falseCount[cIdx]--;

            for (int c2 = 1; c2 <= numColors; ++c2)
            {
                if (c2 == col)
                    continue;
                if (levelOfRemoval[v][c2] >= assignLvl)
                {
                    availableColors[v][c2] = true;
                    levelOfRemoval[v][c2] = -1;
                    reasonForRemoval[v][c2] = Reason();

                    for (int cIdx : clauseIndicesPos[v][c2])
                        if (!clauseSatisfied[cIdx])
                            falseCount[cIdx]--;

                    for (int cIdx : clauseIndicesNeg[v][c2])
                        if (clauseSatisfied[cIdx])
                        {
                            bool anyTrue = false;
                            for (Literal L2 : implicitConstraints[cIdx].literals)
                                if (isLiteralTrue(L2))
                                {
                                    anyTrue = true;
                                    break;
                                }

                            if (!anyTrue)
                                clauseSatisfied[cIdx] = false;
                        }
                }
            }
            isDecisionAssignment[v] = false;
        }

        for (int c = 1; c <= numColors; ++c)
            if (levelOfRemoval[v][c] > targetLevel)
            {
                availableColors[v][c] = true;
                levelOfRemoval[v][c] = -1;
                reasonForRemoval[v][c] = Reason();

                for (int cIdx : clauseIndicesPos[v][c])
                    if (!clauseSatisfied[cIdx])
                        falseCount[cIdx]--;

                for (int cIdx : clauseIndicesNeg[v][c])
                    if (clauseSatisfied[cIdx])
                    {
                        bool anyTrue = false;
                        for (Literal L2 : implicitConstraints[cIdx].literals)
                            if (isLiteralTrue(L2))
                            {
                                anyTrue = true;
                                break;
                            }

                        if (!anyTrue)
                            clauseSatisfied[cIdx] = false;
                    }
            }
    }
    currentLevel = targetLevel;
}

bool cdclGCP(int level)
{
    currentLevel = level;

    bool allColored = true;
    for (int v = 1; v <= G.n; ++v)
        if (colorAssigned[v] == -1)
        {
            allColored = false;
            break;
        }

    if (allColored)
        return true;

    unitPropagation(level);
    if (conflictFound)
    {
        if (level == 0)
            return false;

        auto [learnedClause, backtrackLevel] = analyzeConflict();
        int newClauseIdx = addClause(learnedClause);

        backtrack(backtrackLevel);
        conflictFound = false;
        conflictEmptyClause = false;
        conflictEmptyVertex = -1;
        conflictClauseIndex = -1;
        conflictClause = Clause();

        unitPropagation(currentLevel);
        if (conflictFound)
            return false;

        return cdclGCP(currentLevel);
    }

    int chosenV = -1;
    int minDomain = INT_MAX;
    int maxDeg = -1;
    for (int v = 1; v <= G.n; ++v)
        if (colorAssigned[v] == -1)
        {
            int count = 0;
            for (int c = 1; c <= numColors; ++c)
            {
                if (availableColors[v][c])
                    count++;
            }
            if (count < minDomain || (count == minDomain && (int)G.adj[v].size() > maxDeg))
            {
                minDomain = count;
                maxDeg = G.adj[v].size();
                chosenV = v;
            }
        }

    int currentLvl = level;

    for (int c = 1; c <= numColors; ++c)
    {
        if (!availableColors[chosenV][c])
            continue;

        currentLevel = currentLvl + 1;
        assignColor(chosenV, c, currentLevel, true, Reason());

        // for (int u : G.adj[chosenV])
        // {
        //     if (availableColors[u][c])
        //     {
        //         availableColors[u][c] = false;
        //         levelOfRemoval[u][c] = currentLevel;
        //         reasonForRemoval[u][c] = Reason(chosenV, c);
        //         for (int clauseIdx : clauseIndicesPos[u][c])
        //         {
        //             if (clauseSatisfied[clauseIdx])
        //                 continue;
        //             falseCount[clauseIdx]++;
        //             if (/*falseCount[clauseIdx]*/ countFalse(clauseIdx) == implicitConstraints[clauseIdx].literals.size())
        //             {
        //                 conflictFound = true;
        //                 conflictEmptyClause = true;
        //                 conflictClauseIndex = clauseIdx;
        //                 conflictClause = implicitConstraints[clauseIdx];
        //             }
        //         }
        //     }
        // }

        bool result = cdclGCP(currentLevel);
        if (result)
            return true;

        backtrack(currentLvl);
        conflictFound = false;
        conflictEmptyClause = false;
        conflictEmptyVertex = -1;
        conflictClauseIndex = -1;
    }

    return false;
}

int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    int _n, _m, _k;

    cin >> _n >> _m >> _k;

    emptySolution.clear();
    emptySolution.resize(_n + 1, 0);
    vector<int> A = emptySolution;

    k = _k;

    G.n = _n;
    G.adj.resize(G.n + 1);
    numColors = k;

    while (_m--)
    {
        int u, v;
        cin >> u >> v;
        G.adj[u].push_back(v);
        G.adj[v].push_back(u);
    }

    int n = G.n;
    colorAssigned.assign(n + 1, -1);
    availableColors.assign(n + 1, vector<bool>(k + 1, true));
    reasonForAssignment.assign(n + 1, Reason());
    reasonForRemoval.assign(n + 1, vector<Reason>(k + 1));
    levelOfAssignment.assign(n + 1, -1);
    levelOfRemoval.assign(n + 1, vector<int>(k + 1, -1));
    isDecisionAssignment.assign(n + 1, false);
    currentLevel = 0;
    currentTime = 0;
    conflictFound = false;
    conflictEmptyClause = false;
    conflictClauseIndex = -1;
    conflictEmptyVertex = -1;

    clauseIndicesPos.assign(n + 1, vector<vector<int>>(k + 1));
    clauseIndicesNeg.assign(n + 1, vector<vector<int>>(k + 1));

    if (cdclGCP(0))
    {
        A = colorAssigned;
        for (int i = 1; i <= n; ++i)
            cout << A[i] << " ";
    }
    else
        cout << "No solution found." << endl;
}