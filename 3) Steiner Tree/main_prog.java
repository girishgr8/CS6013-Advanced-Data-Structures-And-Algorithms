import java.util.Scanner;
import java.util.Set;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

// A class which stores the information outputted by Dijkstra 
// which will be used to convert the Metric Steiner Tree Instance back to Steiner Tree Instance
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
        int[][] graph = readAdjacencyMatrix();

        // store the number of nodes in some variable
        int V = graph.length;

        // print the adjacency matrix read from the txt file
        printAdjancencyMatrix(V, graph);

        // read the steiner vertices as an input from the user
        int[] steiner = inputSteinerVertices(V);

        // store the required vertices into another list for future use
        int[] required = new int[V - steiner.length];
        for (int i = 0, r = 0; i < V; i++) {
            boolean isSteiner = false;
            for (int j = 0; j < steiner.length; j++) {
                if (steiner[j] == i)
                    isSteiner = true;
            }
            if (!isSteiner)
                required[r++] = i;
        }

        // convert the Steiner Tree instance I to Metric Steiner Tree instance I'
        int[][] newGraph = convertToMetricSteiner(graph, V);

        // construct a minimum spanning tree on required vertices on the Metric Steiner
        // Tree Instance I'
        int[] parent = findMinimumSpanningTree(required, V, newGraph);

        // System.out.println("parent = " + Arrays.toString(parent));

        // now convert the MST constructed on Metric Steiner Tree back to Steiner Tree
        // this can be done by replacing the direct edge existing in MST but
        // non-existent in the original steiner tree instance
        // this replacement can be done with the help of shortest path between two nodes
        // found by Dijkstra
        int[][] tree = convertBackToSteinerTree(V, required, parent, graph, newGraph);
        int totalCost = 0;

        // check if the MST covers all the required vertices or not
        // this is required as we may have scenario where in the required vertices are
        // in two different connected components
        boolean coversReqdVertices = checkIfMSTCoversReqdVertices(V, tree, required);
        if (!coversReqdVertices) {
            System.out.println(
                    "The Steiner Tree graph instance has disconnected components. The required vertices cannot be connected.");
            System.exit(0);
        }

        // detect and remove cycle on the reverted MST
        detectAndRemoveCycle(V, tree);

        System.out.println(
                "The 2-factor approximate tree we have computed is given below (we describe this tree by listing all the neighbors of all the vertices in the tree): ");

        // print the final node neighbor information from the final MST after doing
        // replacements of non-direct paths
        for (int i = 0; i < V; i++) {
            System.out.printf("Neighbors of Vertex " + (i + 1) + ": ");
            for (int j = 0; j < V; j++) {
                if (tree[i][j] != 0)
                    System.out.printf((j + 1) + " ");
                totalCost += tree[i][j];
            }
            System.out.println();
        }

        // print the cost of the above steiner tree constructed
        System.out.println("\nTotal Cost of given Steiner Tree = " + (totalCost / 2) + "\n");
        sc.close();
    }

    // function which converts MST on Metric Steiner Tree instance to Steiner Tree
    // instance (i.e. original graph)
    public static int[][] convertBackToSteinerTree(int V, int[] required, int[] parent, int[][] steinerGraph,
            int[][] metricSteinerGraph) {

        int[][] tree = new int[V][V];
        for (int p = 1; p < required.length; p++) {
            int i = required[p];
            if (parent[i] == Integer.MIN_VALUE)
                continue;
            // no direct cost path in original steiner tree instance
            // then replace that edge by the stored Dijkstra's shortest path between those
            // two nodes
            if (steinerGraph[parent[i]][i] >= metricSteinerGraph[i][parent[i]] || steinerGraph[parent[i]][i] == 0) {
                int u = parent[i], v = i;
                List<Integer> path = new ArrayList<Integer>();

                Information info = informations[u];
                path.add(v);
                while (u != v) {
                    v = info.predecessor[v];
                    path.add(v);
                }
                int pathSize = path.size();
                for (int j = 0; j < pathSize - 1; j++) {
                    u = path.get(j);
                    v = path.get(j + 1);
                    tree[u][v] = steinerGraph[u][v];
                    tree[v][u] = steinerGraph[u][v];
                }
            }
        }
        return tree;
    }

    private static void detectAndRemoveCycle(int V, int[][] tree) {
        boolean cycle = hasCycle(V, tree);
        if (!cycle)
            return;
    }

    // recursive function to check if cycle exists in the graph or not
    public static boolean hasCycleutil(int V, int[][] tree, boolean[] visited, int u, int parent) {
        visited[u] = true;
        for (int v = 0; v < V; v++) {
            if (tree[u][v] != 0) {
                if (!visited[v]) {
                    if (hasCycleutil(v, tree, visited, v, u))
                        return true;
                }
                if (u != parent)
                    return true;
            }
        }
        return false;
    }

    // utility function to check if cycle exists in the graph or not
    public static boolean hasCycle(int V, int[][] tree) {
        boolean[] visited = new boolean[V];
        for (int u = 0; u < V; u++) {
            if (!visited[u] && hasCycleutil(V, tree, visited, u, -1))
                return true;
        }
        return false;
    }

    // function to check if MST on Steiner Tree covers on required vertices or not
    public static boolean checkIfMSTCoversReqdVertices(int V, int[][] tree, int[] required) {
        Set<Integer> treeVertices = new HashSet<Integer>();

        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (tree[i][j] != 0) {
                    treeVertices.add(i);
                    treeVertices.add(j);
                }
            }
        }

        for (int i = 0; i < required.length; i++) {
            if (!treeVertices.contains(required[i]))
                return false;
        }
        return true;
    }

    // a function to read the adjacency matrix from the input file
    public static int[][] readAdjacencyMatrix() {
        int[][] graph = null;
        try {
            File file = new File("./" + FILENAME);
            BufferedReader br = new BufferedReader(new FileReader(file));

            int V = 0, v = 0;
            while (br.readLine() != null)
                V++;

            br.close();
            br = new BufferedReader(new FileReader(file));

            graph = new int[V][V];
            String line;

            while ((line = br.readLine()) != null) {
                String[] weights = line.split(" ");
                for (int i = 0; i < V; i++)
                    graph[v][i] = Integer.parseInt(weights[i]);
                v++;
            }
            br.close();
        } catch (FileNotFoundException fne) {
            System.out.println("Specified file not found ! " + fne.toString());
        } catch (Exception e) {
            System.out.println("Some exception occured : " + e.toString());
        }
        return graph;
    }

    // a function to print the adjacency matrix
    public static void printAdjancencyMatrix(int V, int[][] graph) {
        System.out.println("The input matrix A the program read from the file is displayed below: ");
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                System.out.printf(graph[i][j] + " ");
            }
            System.out.println();
        }
    }

    // a function to read the steiner vertices from the console until user enters
    // '*'
    public static int[] inputSteinerVertices(int V) {
        Scanner sc = new Scanner(System.in);
        List<Integer> steiner = new ArrayList<Integer>();
        System.out.println("List all the Steiner vertices (type * to quit): ");
        String line = "";
        while (!(line = sc.nextLine()).equals("*")) {
            int vertex = Integer.parseInt(line);
            if (vertex > V)
                System.out.println("Vertex number exceeds maximum vertex number of the graph");
            else
                steiner.add(vertex - 1);
        }
        sc.close();

        int[] arr = new int[steiner.size()];
        for (int i = 0; i < arr.length; i++)
            arr[i] = steiner.get(i);

        return arr;
    }

    public static int[][] convertToMetricSteiner(int[][] graph, int V) {
        int[][] newGraph = new int[V][V];
        informations = new Information[V];
        // since the graph is connected, running dijsktra on every node with that
        // node as source vertex will give make the graph a "Complete Graph"
        // which will obey the triangular inequality which is necessary property in
        // Metric Steiner Tree problem
        for (int i = 0; i < V; i++) {
            informations[i] = dijkstra(graph, V, i);
            for (int k = 0; k < V; k++)
                newGraph[i][k] = informations[i].dist[k];
        }
        return newGraph;
    }

    // a function to get the minimum cost vertex to be selected in the current
    // iteration of Dijkstra
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

    // a function which runs Dijkstra's algorithm which is used to convert Steiner
    // Tree instance to Metric Steiner Tree instance
    public static Information dijkstra(int[][] graph, int V, int src) {
        Information info = new Information(V);
        boolean[] visited = new boolean[V];

        Arrays.fill(info.dist, Integer.MAX_VALUE);
        Arrays.fill(info.predecessor, -1);
        Arrays.fill(visited, false);

        // initialise the distance of src vertex to zero, this node will be picked to
        // start the algorithm
        info.dist[src] = 0;

        // below thing runs in O(V^3) time
        for (int i = 0; i < V - 1; i++) {
            // get minimum cost vertex to be processed
            int u = getMinimumCostVertex(V, info.dist, visited);

            // mark the vertex 'u' as visited, as it will be processed in the current
            // iteration of the Dijsktra's loop
            visited[u] = true;

            for (int v = 0; v < V; v++) {
                int cost = graph[u][v];
                if (!visited[v] && cost != 0 && info.dist[u] != Integer.MAX_VALUE
                        && info.dist[u] + cost < info.dist[v]) {
                    info.dist[v] = info.dist[u] + cost;
                    info.predecessor[v] = u;
                }
            }
        }
        return info;
    }

    // A function to get the minimum key/cost vertex to be selected in the current
    // iteration of Minimum Spanning Tree algorithm
    public static int getMinimumKey(int V, int[] required, int[] key, boolean[] visited) {
        int minCost = Integer.MAX_VALUE, minCostVertex = -1;

        for (int i = 0; i < required.length; i++) {
            int v = required[i];
            if (visited[v] == false && key[v] < minCost) {
                minCost = key[v];
                minCostVertex = v;
            }
        }
        return minCostVertex;
    }

    // a function which finds minimum spanning tree using Kruskal's algorithm
    // but instead of finding a MST on all vertices of I', we will find MST on
    // required vertices of I'
    public static int[] findMinimumSpanningTree(int[] required, int V, int[][] graph) {
        int[] parent = new int[V];
        int[] key = new int[V];
        boolean[] visited = new boolean[V];

        Arrays.fill(key, Integer.MAX_VALUE);
        Arrays.fill(visited, false);
        Arrays.fill(parent, Integer.MIN_VALUE);

        key[required[0]] = 0;
        parent[required[0]] = -1;

        for (int i = 0; i < required.length - 1; i++) {
            int u = getMinimumKey(V, required, key, visited);
            if (u == -1)
                continue;
            visited[u] = true;

            for (int j = 0; j < required.length; j++) {
                int v = required[j];
                int edgeWt = graph[u][v];
                if (edgeWt != 0 && visited[v] == false && edgeWt < key[v]) {
                    parent[v] = u;
                    key[v] = edgeWt;
                }
            }
        }
        return parent;
    }
}