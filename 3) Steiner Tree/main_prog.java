import java.util.Scanner;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

class Information {
    int[] dist;
    int[] predecessor;

    public Information(int V) {
        this.dist = new int[V];
        this.predecessor = new int[V];
    }
}

public class main_prog {

    public static final String FILENAME = "input.txt";
    public static Information[] informations;

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        // read the adjacency matrix from the input txt file
        List<List<Integer>> adjMat = readAdjacencyMatrix();

        // store the number of nodes in some variable
        int V = adjMat.size();

        // print the adjacency matrix read from the txt file
        printAdjancencyMatrix(adjMat);

        // read the steiner vertices as an input from the user
        List<Integer> steinerVertices = inputSteinerVertices(V);

        // store the required vertices into another list for future use
        List<Integer> requiredVertices = new ArrayList<Integer>();
        for (int i = 1; i <= V; i++) {
            if (!steinerVertices.contains(i))
                requiredVertices.add(i);
        }

        // convert the Steiner Tree instance I to Metric Steiner Tree instance I'
        List<List<Integer>> mstpAdjMat = convertToMetricSteiner(adjMat, V);
        // print the adjacency matrix read from the txt file
        printAdjancencyMatrix(mstpAdjMat);
        // System.out.println("---------------------------------------------------------");
        // for (int i = 0; i < V; i++) {
        // Information info = informations[i];
        // System.out.println("For node " + (i + 1) + " : ");
        // for (int j = 0; j < V; j++)
        // System.out.println(info.dist[j] + " " + info.predecessor[j]);
        // System.out.println("---------------------------------------------------------");
        // }

        int[] parent = findMinimumSpanningTree(V, mstpAdjMat);

        // System.out.println("Edge \tWeight");
        // for (int i = 1; i < V; i++)
        //     System.out.println((parent[i] + 1) + " - " + (i + 1) + "\t" + mstpAdjMat.get(i).get(parent[i]));

        int[][] tree = new int[V][V];

        System.out.println("=============================================================");
        System.out.println("Edge \tWeight\tPath");        
        for (int i = 1; i < V; i++) {
            // no direct cost path in original steiner tree instance
            if (adjMat.get(parent[i]).get(i) >= mstpAdjMat.get(i).get(parent[i]) || adjMat.get(parent[i]).get(i) == 0) {
                int u = parent[i], v = i;
                List<Integer> path = new ArrayList<Integer>();
                
                Information info = informations[u];
                path.add(v + 1);
                while (u != v) {
                    v = info.predecessor[v] - 1;
                    path.add(v + 1);
                }
                System.out.println((parent[i] + 1) + " - " + (i + 1) + "\t" + mstpAdjMat.get(i).get(parent[i])
                        + "\tPath = " + path);
                for (int j = 0; j < path.size() - 1; j++) {
                    u = path.get(j) - 1;
                    v = path.get(j + 1) - 1;
                    tree[u][v] = adjMat.get(u).get(v);
                    tree[v][u] = adjMat.get(u).get(v);
                }
            }
        }

        // System.out.println(Arrays.deepToString(tree));
        System.out.println(
                "The 2-factor approximate tree we have computed is given below (we describe this tree by listing all the neighbors of all the vertices in the tree): ");
        for (int i = 0; i < V; i++) {
            System.out.printf("Neighbors of Vertex " + (i + 1) + ": ");
            for (int j = 0; j < V; j++) {
                if (tree[i][j] != 0)
                    System.out.printf((j + 1) + " ");
            }
            System.out.println();
        }
        sc.close();
    }

    public static List<List<Integer>> readAdjacencyMatrix() {
        List<List<Integer>> adjMat = new ArrayList<List<Integer>>();
        try {
            File file = new File("./" + FILENAME);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line;
            while ((line = br.readLine()) != null) {
                String[] weights = line.split(" ");
                int n = weights.length;
                List<Integer> list = new ArrayList<Integer>();
                for (int i = 0; i < n; i++) {
                    list.add(Integer.parseInt(weights[i]));
                }
                adjMat.add(list);
            }
            br.close();
        } catch (FileNotFoundException fne) {
            System.out.println("Specified file not found ! " + fne.toString());
        } catch (Exception e) {
            System.out.println("Some exception occured : " + e.toString());
        }
        return adjMat;
    }

    public static void printAdjancencyMatrix(List<List<Integer>> graph) {
        System.out.println("The input matrix A the program read from the file is displayed below: ");
        int n = graph.size();
        for (int i = 0; i < n; i++) {
            List<Integer> list = graph.get(i);
            for (int j = 0; j < n; j++) {
                System.out.printf(list.get(j) + " ");
            }
            System.out.println();
        }
    }

    public static List<Integer> inputSteinerVertices(int V) {
        Scanner sc = new Scanner(System.in);
        List<Integer> steinerVertices = new ArrayList<Integer>();
        System.out.println("List all the Steiner vertices (type * to quit): ");
        String line = "";
        while (!(line = sc.nextLine()).equals("*")) {
            int vertex = Integer.parseInt(line);
            if (vertex > V)
                System.out.println("Vertex number exceeds maximum vertex number of the graph");
            else
                steinerVertices.add(vertex);
        }
        sc.close();
        return steinerVertices;
    }

    public static List<List<Integer>> convertToMetricSteiner(List<List<Integer>> graph, int V) {
        List<List<Integer>> mstpAdjMat = new ArrayList<List<Integer>>();
        informations = new Information[V];
        for (int i = 0; i < V; i++) {
            informations[i] = dijkstra(graph, V, i);
            List<Integer> list = new ArrayList<Integer>();
            for (int k = 0; k < V; k++) {
                list.add(informations[i].dist[k]);
            }
            mstpAdjMat.add(list);
        }
        return mstpAdjMat;
    }

    public static int getMinimumCostVertex(int V, int[] dist, boolean[] visited) {
        int minDist = Integer.MAX_VALUE, minCostVertex = -1;

        for (int v = 0; v < V; v++) {
            if (visited[v] == false && dist[v] <= minDist) {
                minDist = dist[v];
                minCostVertex = v;
            }
        }
        return minCostVertex;
    }

    public static Information dijkstra(List<List<Integer>> graph, int V, int src) {
        Information info = new Information(V);
        boolean[] visited = new boolean[V];
        Arrays.fill(info.dist, Integer.MAX_VALUE);
        Arrays.fill(info.predecessor, -1);
        Arrays.fill(visited, false);

        info.dist[src] = 0;

        // below thing runs in O(V^3) time
        for (int i = 0; i < V - 1; i++) {
            // get minimum cost vertex to be processed
            int u = getMinimumCostVertex(V, info.dist, visited);

            // mark the vertex 'u' as visited, as it will be processed in the current
            // iteration of the Dijsktra's loop
            visited[u] = true;

            for (int v = 0; v < V; v++) {
                int cost = graph.get(u).get(v);
                if (!visited[v] && cost != 0 && info.dist[u] != Integer.MAX_VALUE
                        && info.dist[u] + cost < info.dist[v]) {
                    info.dist[v] = info.dist[u] + cost;
                    info.predecessor[v] = (u + 1);
                }
            }
        }
        return info;
    }

    public static int getMinimumKey(int V, int[] key, boolean[] visited) {
        int minCost = Integer.MAX_VALUE, minCostVertex = -1;

        for (int v = 0; v < V; v++) {
            if (visited[v] == false && key[v] < minCost) {
                minCost = key[v];
                minCostVertex = v;
            }
        }
        return minCostVertex;
    }

    public static int[] findMinimumSpanningTree(int V, List<List<Integer>> graph) {
        int[] parent = new int[V];
        int[] key = new int[V];
        boolean[] visited = new boolean[V];
        Arrays.fill(key, Integer.MIN_VALUE);
        Arrays.fill(visited, false);

        key[0] = 0;
        parent[0] = -1;

        for (int i = 0; i < V - 1; i++) {
            int u = getMinimumKey(V, key, visited);
            visited[u] = true;

            for (int v = 0; v < V; v++) {
                int edgeWt = graph.get(u).get(v);
                if (edgeWt != 0 && visited[v] == false && edgeWt < key[v]) {
                    parent[v] = u;
                    key[v] = edgeWt;
                }
            }
        }

        return parent;
    }
}