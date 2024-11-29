package utils;


import java.util.*;



public class GraphGenerator {
    Random R = new Random();
    public Graph genStronglyConnectedComponent(int fromNodeId, int toNodeId, int nbEdges){
        Graph G = new DirectedGraph(fromNodeId, toNodeId);
        //G.nbNodes = nbNodes;
        //G.nbEdges = nbEdges;
        Random R = new Random();
        List<Integer> nodes = new ArrayList<Integer>();
        List<Edge> edges = new ArrayList<>();
        Set<Integer>[] A = new Set[toNodeId+1];
        for(int v = fromNodeId; v <= toNodeId; v++){
            A[v] = new HashSet<Integer>();
        }
        int n1 = G.getNodeIdAt(0);
        int n2 = G.getNodeIdAt(1);
        int n3 = G.getNodeIdAt(2);
        Edge e1 = new Edge(n1,n2); A[n1].add(n2);
        Edge e2 = new Edge(n2,n3); A[n2].add(n3);
        Edge e3 = new Edge(n3,n1); A[n3].add(n1);
        G.addEdge(n1,n2); G.addEdge(n2,n3); G.addEdge(n3,n1);
        nodes.add(n1); nodes.add(n2); nodes.add(n3);
        edges.add(e1); edges.add(e2); edges.add(e3);
        for(int k = 3; k < G.nodes.size(); k++){
            int v = G.getNodeIdAt(k);
            while(true){
                int i = R.nextInt(nodes.size()-1);
                int j = R.nextInt(nodes.size()-i-1)+i+1;
                int ni = nodes.get(i);
                int nj = nodes.get(j);
                Edge eiv = new Edge(ni,v); edges.add(eiv); A[ni].add(v);
                Edge evj = new Edge(v,nj); edges.add(evj); A[v].add(nj);
                G.addEdge(ni,v); G.addEdge(v,nj);
                nodes.add(v);
                if(i < j) break;
            }
        }
        while(edges.size() < nbEdges){
            //int v = R.nextInt(nbNodes) + 1;
            int i = R.nextInt(G.nodes.size());
            int v = G.getNodeIdAt(i);
            while(true){
                int idx = R.nextInt(G.nodes.size());
                int u = G.getNodeIdAt(idx);
                if(!A[v].contains(u) && !A[u].contains(v) && u != v){
                    A[v].add(u);
                    edges.add(new Edge(v,u));
                    G.addEdge(v,u);
                    break;
                }
            }
        }

        return G;
    }
    public Graph genNotOrientableGraph(int nbNodes, int nbEdges){
        Graph G = genStronglyConnectedComponent(1,nbNodes-1, nbEdges - 1);
        Graph G1 = new UndirectedGraph(1,nbNodes);
        for(Edge e: G.E){
            G1.addEdge(e.begin,e.end);
        }
        G1.addEdge(1,nbNodes);
        return G1;
    }
    public static void main(String[] args){
        GraphGenerator GG = new GraphGenerator();
        //Graph G = GG.genStronglyConnectedComponent(100000,1000000);
        Graph G = GG.genStronglyConnectedComponent(1,10,20);
        //Graph G = GG.genNotOrientableGraph(1000,20000);

        G.print();
        //G.printFile("C:\\DungPQ\\teaching\\Applied-algorithms\\exercises\\graph-orientation\\test\\5.txt");
    }
}
