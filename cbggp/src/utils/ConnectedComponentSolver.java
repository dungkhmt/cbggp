package utils;

import java.util.*;

public class ConnectedComponentSolver {
    private int[] p;// parent
    private List<Edge> E;
    private Set<Integer> V;
    private Map<Integer, Set<Integer>> A;
    private int nbCC; // number of connected components;
    private void dfs(int u){
        for(int v: A.get(u)){
            if(p[v] == -1){
                p[v] = u;
                dfs(v);
            }
        }
    }
    public int computeNumberConnectedComponents(List<Edge> E){
        if(E == null || E.size() == 0){
            return 0;
        }
        this.E = E;
        V = new HashSet<>();
        int minNode = (int)1e9; int maxNode = 1-minNode;
        for(Edge e: E){
            minNode = Math.min(minNode,e.begin);
            minNode = Math.min(minNode,e.end);
            maxNode = Math.max(maxNode,e.begin);
            maxNode = Math.max(maxNode,e.end);
            V.add(e.begin); V.add(e.end);
        }
        //System.out.println(minNode + " " + maxNode);
        p = new int[maxNode + 1];
        for(int v = minNode; v <= maxNode; v++){
            p[v] = -1;// not visited
        }
        A = new HashMap();
        for(int v: V){
            A.put(v,new HashSet<Integer>());
        }
        for(Edge e: E){
            int u = e.begin; int v = e.end;
            A.get(u).add(v); A.get(v).add(u);
        }
        nbCC = 0;
        for(int v: V){
            if(p[v] == -1){
                nbCC ++; dfs(v);
            }
        }
        return nbCC;
    }
    public static void main(String[] args){

    }
}
