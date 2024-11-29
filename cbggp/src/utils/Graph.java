package utils;

import java.io.PrintWriter;
import java.util.*;

class Edge{
    int begin;
    int end;
    public Edge(int begin, int end){
        this.begin = begin; this.end = end;
    }
    public int getOtherNode(int u){
        if(this.begin == u) return end;
        else if(this.end == u) return this.begin;
        else return -1;
    }
    public String toString(){
        return "(" + begin + "," + end + ")";
    }
}
public abstract class Graph{
    int nbNodes;
    int nbEdges;
    List<Edge> E;
    List<Integer> nodes;
    List<Integer> A[]; // A[v] is the list of indices of adjacent edges (outgoing) of v
    public Graph(List<Edge> edges){
        E = new ArrayList<Edge>();

        nodes = new ArrayList<Integer>();
        Set<Integer> S = new HashSet<Integer>();
        for(Edge e: edges){
            int u = e.begin; int v = e.end;
            E.add(new Edge(u,v)); S.add(u); S.add(v);
        }
        int maxNodeId = -1;
        for(int v: S) maxNodeId = Math.max(maxNodeId, v);
        for(int v: S) nodes.add(v); Collections.sort(nodes);
        nbNodes = nodes.size();
        A = new ArrayList[maxNodeId + 1];
        for(int v = 0; v <= maxNodeId; v++) A[v] = new ArrayList<Integer>();
    }
    public Graph(int fromNodeId, int toNodeId){
        nodes = new ArrayList<Integer>();
        for(int v = fromNodeId; v <= toNodeId; v++) nodes.add(v);
        nbNodes = toNodeId - fromNodeId + 1; nbEdges = 0;
        A = new List[toNodeId + 1];
        for(int v = fromNodeId; v <= toNodeId; v++){
            A[v] = new ArrayList<Integer>();
        }
        E = new ArrayList<Edge>();
    }
    public void print(){
        System.out.println(nbNodes + " " + E.size());
        for(Edge e: E){
            System.out.println(e.begin + " " + e.end);
        }
    }
    abstract public boolean hasEdge(int u, int v);
    abstract public boolean addEdge(int u, int v);
    public int getNodeIdAt(int i){
        return nodes.get(i);
    }
    public void printFile(String filename){
        try{
            PrintWriter out = new PrintWriter(filename);
            out.println(nbNodes + " " + nbEdges);
            for(Edge e: E){
                out.println(e.begin + " " + e.end);
            }
            out.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
}

class UndirectedGraph extends Graph{

    public UndirectedGraph(List<Edge> edges) {
        super(edges);
        for(int i = 0; i < E.size(); i++){
            Edge e = E.get(i);
            int u = e.begin; int v = e.end;
            A[u].add(i);
            A[v].add(i);
        }

    }

    public UndirectedGraph(int fromNodeId, int toNodeId) {
        super(fromNodeId, toNodeId);
    }

    @Override
    public boolean hasEdge(int u, int v) {
        for(int i: A[u]){
            Edge e = E.get(i);
            if(e.getOtherNode(u) == v) return true;
        }
        return false;
    }

    public boolean addEdge(int u, int v){
        if(u == v) return false;
        if(hasEdge(u,v)) return false;
        if(u < v) E.add(new Edge(u,v)); else E.add(new Edge(v,u));
        nbEdges += 1;
        A[u].add(E.size()-1);
        A[v].add(E.size()-1);
        return true;

    }
}

class DirectedGraph extends Graph{

    public DirectedGraph(List<Edge> edges) {
        super(edges);
        for(int i = 0; i < E.size(); i++){
            Edge e = E.get(i);
            int u = e.begin; int v = e.end;
            A[u].add(i);
        }
    }

    public DirectedGraph(int fromNodeId, int toNodeId) {
        super(fromNodeId, toNodeId);
    }

    @Override
    public boolean hasEdge(int u, int v) {
        for(int i: A[u]){
            Edge e = E.get(i);
            if(e.getOtherNode(u)==v) return true;
        }
        return false;
    }

    @Override
    public boolean addEdge(int u, int v) {
        if(u == v) return false;
        if(hasEdge(u,v)) return false;
        E.add(new Edge(u,v));
        A[u].add(E.size()-1);
        return true;
    }
}
