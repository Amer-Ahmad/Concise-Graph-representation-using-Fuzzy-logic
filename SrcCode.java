import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Scanner;
import java.util.Stack;

import mdsj.StressMinimization;
import net.sourceforge.jFuzzyLogic.FIS;
import net.sourceforge.jFuzzyLogic.rule.Variable;

public class SrcCode {
    static double[][] X;
    static int[][] PA;
    static double[] r;
    static double[] R;
    static double defCount = 0;
    static int col = -1;
    static int n = 0;
    static int o1, o2;
    static int k = 50;
    static double d[][];
    static int adj[][];
    static FIS fis;
    final static double eps = 0.000001;

    static double countYesAboveHalf = 0;
    static double countFuzzyYes = 0;
    static double sumFuzzyYesInput = 0;
    static double countNoBelowHalf = 0;
    static double countFuzzyNo = 0;
    static double sumFuzzyNoInput = 0;
    static double correctAnswers = 0;

    public static void main(String[] args) throws IOException {
        prepareFuzzyInferenceSystem();
//        String pathname = System.getProperty("input.txt");
//        Scanner scan = new Scanner(new File(pathname));
//        n = scan.nextInt();
//        int e = scan.nextInt();
//        d = new double[n][n];
//        for (int i = 0; i < n; i++)
//            for (int j = 0; j < n; j++) {
//                if (i == j) continue;
//                d[i][j] = n;
//            }
//        for (int i = 0; i < e; i++) {
//            int a = scan.nextInt() - 1;
//            int b = scan.nextInt() - 1;
//            d[a][b] = 1;
//            d[b][a] = 1;
//        }
        
        
       
        
        File input = new File("C:\\Users\\USP\\eclipse-workspace\\BigGraphsResearch\\src\\socfb-Haverford76.txt");
        Scanner scan = new Scanner(input);
        
        n = scan.nextInt();
        int e = scan.nextInt();
        d = new double[n][n];
        adj = new int[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                d[i][j] = adj[i][j] = n;
            }
        while(scan.hasNext())
        {
        	for (int i = 0; i < e; i++) {
        		int a = scan.nextInt() - 1;
                int b = scan.nextInt() - 1;
                //int temp = scan.nextInt() - 1;
                //System.out.println(a + " " + b);
                d[a][b] = 1;
                d[b][a] = 1;
                adj[a][b] = adj[b][a] = 1;
        	}
        }
        
        
        for ( int i=0 ; i<n ; i++ )
        {
        	Queue <Integer> Q = new LinkedList<>();
        	Q.add(i);
        	int vis[] = new int[n];
        	while( !Q.isEmpty() )
        	{
        		int x = Q.poll();
        		if (vis[x]==1) { continue; }
        		vis[x]=1;
        		for ( int j=0 ; j<n ; j++ )
        		{
        			if ( d[x][j]==1 )
        			{
        				d[i][j] = Math.min(d[i][j], 1+d[i][x]);
        				Q.add(j);
        			}
        		}
        	}
        }
        
        r = new double[n];
        R = new double[n];

        long start, end;
        long duration;
        X = new double[n][k];
        PA = new int[k][2];
        start = System.currentTimeMillis();
        FastMap(k, d);
        System.out.println( "after fast map" );
        calculateRs(X);
        end = System.currentTimeMillis();
        duration = end - start;
        System.out.println("FastMap Duration: " + duration);

        double[][] XT = transposeMatrix(X);
        StressMinimization stressMinimization = new StressMinimization(d, XT);
        System.out.println("FastMap stress: " + stressMinimization.getNormalizedStress());
        checkIfError(X);
        printResultsSummary(e);

//        start = System.currentTimeMillis();
//        double[][] res = MDSJ.classicalScaling(d, k);
//        double[][] resT = transposeMatrix(res);
//        calculateRs(resT);
//        end = System.currentTimeMillis();
//        duration = end - start;
//        System.out.println("MDS duration: " + duration);
//
//        countYesAboveHalf = 0;
//        countFuzzyYes = 0;
//        sumFuzzyYesInput = 0;
//        countNoBelowHalf = 0;
//        countFuzzyNo = 0;
//        sumFuzzyNoInput = 0;
//        defCount = 0;
//        StressMinimization stressMinimization2 = new StressMinimization(d, res);
//        System.out.println("MDS stress: " + stressMinimization2.getNormalizedStress());
//        checkIfError(resT);
//        printResultsSummary(e);
    }

    private static void printResultsSummary(int e) {
        System.out.println("*******************************");
        System.out.println( "nodes: " + n + "  Edges: " + e);
        System.out.println("k: " + k);
        System.out.println( "fuzzyYes: " + countFuzzyYes + " fuzzyNo: " + countFuzzyNo );
        System.out.println("Definitive answers: " + ((defCount) / (n * (n - 1) / 2)) * 100 + "%");
        System.out.println("Fuzzy accuracy: " + ((countYesAboveHalf + countNoBelowHalf) / (countFuzzyYes + countFuzzyNo)) * 100 + "%");
        System.out.println("correct Answers (definitive + fuzzy): " + (correctAnswers / ( n * (n-1) / 2 ))*100 + "%" );
        //printNodesProjections();
        System.out.println("*******************************");
    }

    public static double[][] transposeMatrix(double[][] m) {
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }

    public static void prepareFuzzyInferenceSystem() {
        String fileName = "C:\\Users\\USP\\eclipse-workspace\\BigGraphsResearch\\src\\Rules.fcl";
        fis = FIS.load(fileName, true);

        if (fis == null) {
            System.err.println("Can't load file: '" + fileName + "'");
            return;
        }
        
    }

    public static boolean checkIfError(double[][] X) {
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                double status = areNeighbors(X, i, j);
//                System.out.println( i + " ------------ " + j + "\n" );
                if ((adj[i][j] > 1 && status <= 0.5) || (adj[i][j] == 1 && status > 0.5)) {
                	correctAnswers++;
//                    return true;
                }
            }
        return false;
    }

    public static void FastMap(int k, double[][] d) {
        if (k <= 0)
            return;
        else
            col++;
        chooseDistantObjects(d);
        PA[col][0] = o1;
        PA[col][1] = o2;
        if (d[o1][o2] == 0) {
            for (int i = 0; i < n; i++)
                X[i][col] = 0;
            return;
        }
        for (int i = 0; i < n; i++) {
            double x = (Math.pow(d[o1][i], 2) + Math.pow(d[o1][o2], 2) - Math.pow(d[o2][i], 2)) / (2 * d[o1][o2]);
            X[i][col] = x;
        }
        double[][] dprime = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                dprime[i][j] = Math.sqrt(Math.abs(Math.pow(d[i][j], 2) - Math.pow(X[i][col] - X[j][col], 2)));
                dprime[j][i] = dprime[i][j];
            }
        }
        FastMap(k - 1, dprime);
    }
    
    
    public static void printNodesProjections()
    {
    	for ( int i=0 ; i<n ; i++ )
    	{
    		System.out.print( "node " + (i+1) + ": " );
    		for ( int j=0 ; j<k ; j++ )
    		{
    			System.out.print( X[i][j] + " " );
    		}
    		System.out.println();
    	}
    }
    
    
//    static class Point implements Comparable<Point> {
//    	int index;
//    	double x, y;
//    	double[] dist;
//		@Override
//		public int compareTo(Point other) {
//			if (Double.compare(this.x, other.x) == 0)
//				return Double.compare(this.y, other.y);
//			return Double.compare(this.x, other.x);
//		}
//    }
//    
//    void fun (double[][] d, double[] x, double[] y) {
//    	int n = d.length;
//    	Point[] p = new Point[n];
//    	for (int i = 0; i < n; i++) {
//    		p[i] = new Point();
//    		p[i].index = i;
//    		p[i].dist = d[i];
//    		// TODO add point coordinates
//    		p[i].x = 0;
//    		p[i].y = 0;
//    	}
//    	// get ch
//    	// do 2 pointers
//    	List<Point> hull = convexHull(p);
//    	double dist = 0, 
//    }
//    
//    static boolean cw(Point a, Point b, Point c) {
//    	return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) < 0;
//    }
//    
//    static List<Point> convexHull(Point[] p) {
//    	int n = p.length;
//    	if (n <= 3)
//    		return Arrays.asList(p);
//    	Arrays.sort(p);
//    	Stack<Point> q = new Stack<>();
//    	for (int i = 0; i < n; q.push(p[i++])) {
//    		while (q.size() >= 2) {
//    			Point a = q.pop(), b = q.pop();
//    			q.push(b);
//    			if (cw(b, a, p[i])) {
//    				q.push(a);
//    				break;
//    			}
//    		}
//    	}
//        for (int i = n - 2, t = q.size(); i >= 0; q.push(p[i--])) {
//        	while (q.size() > t) {
//        		Point a = q.pop(), b = q.pop();
//    			q.push(b);
//    			if (cw(b, a, p[i])) {
//    				q.push(a);
//    				break;
//    			}
//        	}
//        }
//        return new ArrayList<>(q);
//    }

    public static void chooseDistantObjects(double[][] d) {
        int secondO = (int) (Math.random() * n);
        int firstO = -1;
        for (int c = 0; c < 5; c++) {
            double max = Integer.MIN_VALUE;
            firstO = -1;
            for (int i = 0; i < n; i++)
                if (max < d[secondO][i]) {
                    max = d[secondO][i];
                    firstO = i;
                }
            max = Integer.MIN_VALUE;
            for (int i = 0; i < n; i++)
                if (max < d[firstO][i]) {
                    max = d[firstO][i];
                    secondO = i;
                }
        }
//        double max = Double.NEGATIVE_INFINITY;
//        for ( int i=0 ; i<n ; i++ )
//        {
//        	for ( int j=i+1 ; j<n ; j++ )
//        	{
//        		if ( Double.compare(d[i][j], max)  > 0) {
//        			max = d[i][j];
//        			firstO = i;
//        			secondO = j;
//        		}
//        	}
//        }
        o1 = firstO;
        o2 = secondO;
        System.out.println( "pivotal points are: " + (o1+1) + " -- " + (o2+1) );
    }

    public static void calculateRs(double[][] X) {
        for (int i = 0; i < n; i++) {
            Distance[] distances = new Distance[n];
            for (int j = 0; j < n; j++) {
                distances[j] = new Distance( calculateDistance(X, i, j), j);
            }
            Arrays.sort(distances);
            boolean notAssigned = true;
            int target = -1;
            if (adj[i][distances[0].fromnode] > 1) {
                r[i] = Integer.MIN_VALUE;
            } else {
                for (int j = 1; j < n; j++) {
                    if (adj[i][distances[j].fromnode] > 1) {
                        target = j - 1;
                        double temp = distances[target].dist;
                        if (target == 0 && distances[j].dist == temp) {
                            r[i] = Integer.MIN_VALUE;
                            notAssigned = false;
                        }
                        while (target > 0 && distances[j].dist == temp) {
                            target--;
                            temp = distances[target].dist;
                            if (target == 0 && distances[j].dist == temp) {
                                r[i] = Integer.MIN_VALUE;
                                notAssigned = false;
                            }
                        }
                        break;
                    }
                }
                if (target == -1) {
                    r[i] = Integer.MAX_VALUE;
                } else {
                    if (notAssigned)
                        r[i] = distances[target].dist;
                }
            }
            target = -1;
            notAssigned = true;
            if (adj[i][distances[n - 1].fromnode] == 1) {
                R[i] = Integer.MAX_VALUE;
            } else {
                for (int j = n - 2; j >= 0; j--) {
                    if (adj[i][distances[j].fromnode] == 1) {
                        target = j + 1;
                        double temp = distances[target].dist;
                        if (target == n - 1 && distances[j].dist == temp) {
                            R[i] = Integer.MAX_VALUE;
                            notAssigned = false;
                        }
                        while (target < n - 1 && distances[j].dist == temp) {
                            target++;
                            temp = distances[target].dist;
                            if (target == n - 1 && distances[j].dist == temp) {
                                R[i] = Integer.MAX_VALUE;
                                notAssigned = false;
                            }
                        }
                        break;
                    }
                }
                if (target == -1) {
                    R[i] = Integer.MIN_VALUE;
                } else {
                    if (notAssigned)
                        R[i] = distances[target].dist;
                }
//                if ( r[i] == Integer.MIN_VALUE && R[i]>-1.0 && R[i]<2000 )
//                {
//                	System.out.println( R[i] + "!!" );
//                	for ( int j=0 ; j<n ; j++ )
//                	{
//                		System.out.print( distances[j].dist + " - " );
//                	}
//                	System.out.println();
//                }
            }
        }
    }

    public static double areNeighbors(double[][] X, int a, int b) {
        double dist = calculateDistance(X, a, b);
        if ( dist-r[a]<eps ||  dist-r[b]<eps) {
            defCount++;
            return 1;
        }
        if ( dist-R[a]>eps || dist-R[b]>eps ) {
            defCount++;
            return 0;
        }
        if (r[a] == Integer.MIN_VALUE && r[b] == Integer.MIN_VALUE && R[a] == Integer.MAX_VALUE && R[b] == Integer.MAX_VALUE)
            return 0.5;
        double deffuzzified1 = Integer.MAX_VALUE;
        double input1 = 0;
        if (R[a] != Integer.MAX_VALUE && r[a] != Integer.MIN_VALUE) {
            input1 = Math.abs(((double) R[a]) - ((double) dist)) / Math.abs(((double) R[a]) - ((double) r[a]));
            fis.setVariable("EstimatedDist", input1);
            fis.evaluate();
            Variable result1 = fis.getVariable("NeighborStatus");
            deffuzzified1 = result1.getValue();
        }
        double deffuzzified2 = Integer.MAX_VALUE;
        double input2 = 0;
        if (R[b] != Integer.MAX_VALUE && r[b] != Integer.MIN_VALUE) {
            input2 = Math.abs(((double) R[b]) - ((double) dist)) / Math.abs(((double) R[b]) - ((double) r[b]));
            fis.setVariable("EstimatedDist", input2);
            fis.evaluate();
            Variable result2 = fis.getVariable("NeighborStatus");
            deffuzzified2 = result2.getValue();
        }

        double result = Math.min(deffuzzified1, deffuzzified2);
        result = (result == Integer.MAX_VALUE) ? 0.5 : result;
        if (d[a][b] == 1) {
            countFuzzyYes++;
            countYesAboveHalf += ((result > 0.5) || dist==0.0 ) ? 1 : 0;
            sumFuzzyYesInput += (deffuzzified1 < deffuzzified2) ? input1 : input2;
//            if ( result<=0.5 && dist!=0.0 )
//            {
//            	System.out.println( "-----------------------------------------------------" );
//            	System.out.println( "input1: " + input1 + " ||||||| input2: " + input2 + " ||||||| deffuzzified1: " + deffuzzified1 + " ||||||| deffuzzified2: " + deffuzzified2 );
//            	System.out.println( "nodes :" + a + " --- & --- " + b + " are NEGIHBORS but got wrongly fuzzied!" );
//            	for ( int i=0 ; i<2 ; i++ )
//            	{
//            		System.out.print( "node " + (i+1) + ": " );
//            		int p = (i==0?a:b);
//            		for ( int j=0 ; j<k ; j++ )
//            		{
//            			System.out.print( X[p][j] + " " );
//            		}	
//            		System.out.println();
//            		System.out.println( "r: " + r[p] + " ||| R: " + R[p] );
//            	}
//            	System.out.println( "dist between a & b is :" + dist );
//            	System.out.println( "result: " + result );
//            	System.out.println( "-----------------------------------------------------" );
//            }
//            else
//            {
//            	System.out.println( " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! result is true:: " + result + "!!!!\n"  );
//            }
        } else {
            countFuzzyNo++;
            countNoBelowHalf += (result <= 0.5) ? 1 : 0;
//            if ( result>0.5 )
//            {
//            	System.out.println( "nodes :" + a + " --- & --- " + b + " are *NON*-neighbors but got wrongly fuzzied!" );
//            	for ( int i=0 ; i<2 ; i++ )
//            	{
//            		System.out.print( "node " + (i+1) + ": " );
//            		int p = (i==0?a:b);
//            		for ( int j=0 ; j<k ; j++ )
//            		{
//            			System.out.print( X[p][j] + " " );
//            		}
//            		System.out.println();
//            		System.out.println( "r: " + r[p] + " ||| R: " + R[p] );
//            	}
//            	System.out.println( "dist between a & b is :" + dist );
//            	System.out.println( "result: " + result );
//            }
            sumFuzzyNoInput += (deffuzzified1 < deffuzzified2) ? input1 : input2;
        }
        return result;
    }

    private static double calculateDistance(double[][] X, int a, int b) {
        double dist = 0;
        for (int i = 0; i < k; i++) {
            dist += Math.pow(X[a][i] - X[b][i], 2);
        }
        return Math.sqrt(dist);
    }
}

class Distance implements Comparable<Distance> {
    double dist;
    int fromnode;

    public Distance(double dist, int fromnode) {
        this.dist = dist;
        this.fromnode = fromnode;
    }

    @Override
    public int compareTo(Distance o) {
    	//int diff = (int) ((int)this.dist - o.dist);
        return Double.compare(this.dist , o.dist);
    }
}