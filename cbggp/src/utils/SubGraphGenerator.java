package utils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SubGraphGenerator {
    Graph UP;
    List<Edge> E; //  list of edges of the upper bound graph
    Set<Integer> V;// set of nodes
    int nbNodes;

    int nbEdges;
    int[] x; // decision variables x[i] = 1 means that edge i is selected
    int[] appear;// appear[v] is the number of occurences of node v in selected edges
    int nbSelectedEdges;
    Graph result = null;
    int maxNbSolutions = 1;
    List<Graph> solutions = null;
    private boolean check(int v, int k){
        if(nbSelectedEdges + v > nbEdges) return false;
        return true;
    }
    private void solution(){

        ConnectedComponentSolver CCS = new ConnectedComponentSolver();
        List<Edge> SE = new ArrayList<Edge>();
        for(int k = 0; k < E.size(); k++){
            if(x[k] == 1){
                //System.out.print(E.get(k).toString() + " ");
                SE.add(E.get(k));
            }
        }
        //System.out.println();
        int nbCC = CCS.computeNumberConnectedComponents(SE);
        int n = 0;
        for(int v: V){
            if(appear[v] > 0) n+= 1;
        }
        //System.out.println("Solution, nbNodes = " + n + " nbEdges = " + SE.size() + " nbCC = " + nbCC);

            if(n == nbNodes && SE.size() == nbEdges){
                //System.out.println("FOUND");
                result = new UndirectedGraph(SE);
                solutions.add(result);
            }

    }
    private void tryValue(int k){
        if(solutions.size() >= maxNbSolutions) return;
        // try values for x[k]
        for(int v = 1; v >= 0 ; v--){
            if(check(v,k)){
                x[k] = v;
                nbSelectedEdges += v;
                Edge e = E.get(k);
                if(v == 1){
                    appear[e.begin] += 1;
                    appear[e.end] += 1;
                }
                if(k == E.size()-1 || nbSelectedEdges == nbEdges){
                    solution();
                }else{
                    tryValue(k+1);
                }
                nbSelectedEdges -= v;
                if(v == 1){
                    appear[e.begin] -= 1;
                    appear[e.end] -= 1;
                }

            }
        }
    }
    Graph genConnectedSubGraph(int nbNodes, int nbEdges, Graph G, int nbSolutions){
        solutions = new ArrayList<>();
        this.maxNbSolutions = nbSolutions;
        // gen a connected subgraph of G having nbNodes nodes and nbEdges edges
        this.nbEdges = nbEdges; this.nbNodes = nbNodes;
        this.E = G.E;
        V = new HashSet<Integer>();
        int maxNodeId = -1;
        for(Edge e: E){
            V.add(e.begin);
            V.add(e.end);
            maxNodeId = Math.max(maxNodeId,e.begin);
            maxNodeId = Math.max(maxNodeId,e.end);
        }
        appear = new int[maxNodeId + 1];
        x = new int[E.size()];
        result = null;
        tryValue(0);
        return result;
    }
    public static void main(String[] args){
        SubGraphGenerator SGG = new SubGraphGenerator();
        int n = 5;
        UndirectedGraph G = new UndirectedGraph(1,n);
        for(int i = 1; i <= n-1; i++){
            for(int j =  i+1; j <= n; j++)
                G.addEdge(i,j);
        }
        Graph sg = SGG.genConnectedSubGraph(n,n-1,G, 500);
        //sg.print();
        for(Graph g: SGG.solutions){
            g.print();
            System.out.println("-------------");
        }
        System.out.println("number solutions = " + SGG.solutions.size());
        G.print();
    }
}
