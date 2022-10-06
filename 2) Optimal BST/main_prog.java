
/* 											CS6013: ADVANCED DATA STRUCTURES AND ALGORITHMS (ADSA)
											Programming Assignment II - OPTIMAL BINARY SEARCH TREE

NAME : 			THATTE GIRISH MAKARAND
ROLL NUMBER : 	CS22MTECH11005
 
*/

// imported required libraries
import java.util.*;

// Java class to store the word and its probaility of occurrence
class Word {
    String word;
    Double prob;

    public Word(String word, Double prob) {
        this.word = word;
        this.prob = prob;
    }
}

class main_prog {
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        // Take the required inputs
        System.out.printf("How many strings do you want to insert in the BST? ");
        int n = sc.nextInt();
        sc.nextLine();
        Word[] words = new Word[n];
        System.out.println("Enter " + n + " strings in sorted dictionary order along with their probabilities:");
        for (int i = 0; i < n; i++) {
            String word = sc.nextLine();
            Double prob = Double.parseDouble(sc.nextLine());
            words[i] = new Word(word, prob);
        }

        // Check if input is valid or not
        // 1. Words should be in sorted(dictionary) order
        // 2. All the probabilities should be distinct
        // 3. The sum of probabilities for all words should add up to 1
        boolean invalidInput = false;
        if (!areWordsInSortedOrder(words)) {
            System.out.println("The strings entered are not in sorted order.");
            invalidInput = true;
        }
        if (!areProbabilitiesDistinct(words)) {
            System.out.println("The probabilities are not distinct.");
            invalidInput = true;
        }
        if (!areProbabilitiesSummingToOne(words)) {
            System.out.println("The probabilities don’t add up to 1.");
            invalidInput = true;
        }

        if (invalidInput)
            System.exit(1);

        // create a dp table to store the intermediatary results and to solve the
        // problem in bottom up fashion
        double[][] dp = new double[n + 1][n + 1];
        int[][] split = new int[n + 1][n + 1];

        // call the optimalBinarySearchTree() function to get the optimal cost
        double minCost = optimalBinarySearchTree(words, dp, split);
        System.out.printf("The minimum expected total access time is %.2f\n", minCost);

        // perform the pre-order traversal of the obtained optimal BST
        System.out.println("Preorder traversal of the BST that provides minimum expected total access time is:");
        doPreorderTraversal(split, words, 0, n);
        System.out.println();

        sc.close();
    }

    // doPreorderTraversal() function does the preorder traversal of the optimal
    // tree.
    // here we can use the split[][] matrix which stored the optimal split for the
    // root node
    private static void doPreorderTraversal(int[][] split, Word[] words, int i, int j) {
        if (i == j)
            return;
        int s = split[i][j];
        System.out.printf(words[s - 1].word + " ");
        doPreorderTraversal(split, words, i, s - 1);
        doPreorderTraversal(split, words, s, j);
    }

    // areWordsInSortedOrder() function checks if all the words inputted are in
    // sorted order or not
    private static boolean areWordsInSortedOrder(Word[] words) {
        for (int i = 0; i < words.length - 1; i++) {
            String word1 = words[i].word;
            String word2 = words[i + 1].word;
            if (word1.compareTo(word2) > 0)
                return false;
        }
        return true;
    }

    // areProbabilitiesDistinct() function checks if all the probabilities are
    // distinct or not
    private static boolean areProbabilitiesDistinct(Word[] words) {
        int n = words.length;
        Set<Double> set = new HashSet<Double>();

        for (Word word : words)
            set.add(word.prob);

        return (set.size() == n);
    }

    // areProbabilitiesSummingToOne() function checks if all the probabilities sum
    // up to 1 or not
    private static boolean areProbabilitiesSummingToOne(Word[] words) {
        double sumOfProb = 0;

        for (Word word : words)
            sumOfProb += word.prob;

        return (Math.abs(sumOfProb - 1f) < 1e-9) ? true : false;
    }

    // optimalBinarySearchTree() function computes the optimal(minimum) total access
    // cost/time
    private static double optimalBinarySearchTree(Word words[], double[][] dp, int[][] split) {
        int n = words.length;

        // dp[i][i + 1] considers all chains of length = 1
        // e.g.: dp[0][1], dp[1][2], dp[2][3], dp[3][4],..... these all means our
        // problem consists of just a single node which will directly be the root node
        for (int i = 0; i < n; i++)
            dp[i][i + 1] = words[i].prob;

        // dp[i][j] considers all chains of length = j - i
        // e.g.: dp[0][2], dp[1][3], dp[2][4],..... these all means our
        // problem consists of two nodes out of which one will be root node and other
        // will be it child node
        // with 2 nodes we have 2 possible configurations, we check which configuration
        // gives minimum cost

        for (int chainLen = 2; chainLen <= n + 1; chainLen++) {
            for (int chainStart = 0; chainStart <= n - chainLen + 1; chainStart++) {

                int chainEnd = chainStart + chainLen - 1;
                dp[chainStart][chainEnd] = Double.MAX_VALUE;
                // The DP formula for computing optimal cost is:
                // c[i][j] = min(i < k <= j){c[i][k-1] + c[k][j] + wt[i,j]}
                // here wt[i,j] is cumulative sum of probabilities from ith word to jth word
                double weight_i_j = sumOverRange(words, chainStart, chainEnd);

                // try out each possible split between i and j,
                // i.e. make kth node as the root node,
                // nodes(i, k-1) as left children and nodes(k, j) as right children
                for (int k = chainStart + 1; k <= chainEnd; k++) {
                    double c = dp[chainStart][k - 1] + dp[k][chainEnd] + weight_i_j;
                    // if we got the lower cost split, then store it in our DP table and also store
                    // the split, i.e. making which node as root gave us lower cost
                    if (c < dp[chainStart][chainEnd]) {
                        dp[chainStart][chainEnd] = c;
                        split[chainStart][chainEnd] = k;
                    }
                }
            }
        }
        // return the final optimal cost of accessing all the words with given
        // probabilities
        return dp[0][n];
    }

    // sumOverRange() function calculates cumulative sum of probabilities of words
    // from index low to high
    private static double sumOverRange(Word[] words, int low, int high) {
        double sum = 0;
        for (int i = low; i < high; i++)
            sum += words[i].prob;
        return sum;
    }
}