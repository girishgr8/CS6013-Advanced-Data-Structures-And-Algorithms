import java.util.Scanner;
import java.util.ArrayList;
import java.util.Arrays;
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

    public static final String FILENAME = "input4.txt";
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

        int[] parent = findMinimumSpanningTree(required, V, newGraph);
        int[][] tree = new int[V][V];
        int totalCost = 0;

        for (int p = 1; p < required.length; p++) {
            int i = required[p];
            // no direct cost path in original steiner tree instance
            if (graph[parent[i]][i] >= newGraph[i][parent[i]] || graph[parent[i]][i] == 0) {
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
                    tree[u][v] = graph[u][v];
                    tree[v][u] = graph[u][v];
                }
            }
        }

        System.out.println(
                "The 2-factor approximate tree we have computed is given below (we describe this tree by listing all the neighbors of all the vertices in the tree): ");
        for (int i = 0; i < V; i++) {
            System.out.printf("Neighbors of Vertex " + (i + 1) + ": ");
            for (int j = 0; j < V; j++) {
                if (tree[i][j] != 0)
                    System.out.printf((j + 1) + " ");
                totalCost += tree[i][j];
            }
            System.out.println();
        }

        System.out.println("\nTotal Cost of given Steiner Tree = " + (totalCost / 2) + "\n");
        sc.close();
    }

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

    public static void printAdjancencyMatrix(int V, int[][] graph) {
        System.out.println("The input matrix A the program read from the file is displayed below: ");
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                System.out.printf(graph[i][j] + " ");
            }
            System.out.println();
        }
    }

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
        for (int i = 0; i < V; i++) {
            informations[i] = dijkstra(graph, V, i);
            for (int k = 0; k < V; k++)
                newGraph[i][k] = informations[i].dist[k];
        }
        return newGraph;
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

    public static Information dijkstra(int[][] graph, int V, int src) {
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