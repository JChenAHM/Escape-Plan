import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.Scanner;

public class EscapePlan {
	public int V = 0;
	public int E;
	public int X, S;
	final int capacity = 1;
	public ArrayList<Integer> Xlist = new ArrayList<Integer>(); //X List
	public ArrayList<Integer> Slist = new ArrayList<Integer>(); //S List
	FlowNetwork FN = new FlowNetwork (V);
	
	/**
	 * Bag
	 */
	public static class Bag<Item> implements Iterable<Item> {
	    private int N;               // number of elements in bag
	    private Node<Item> first;    // beginning of bag

	    @SuppressWarnings("hiding")
		private class Node<Item> {
	        private Item item;
	        private Node<Item> next;
	    }

	    public void add(Item item) {
	        Node<Item> oldfirst = first;
	        first = new Node<Item>();
	        first.item = item;
	        first.next = oldfirst;
	        N++;
	    }

	    public Iterator<Item> iterator()  {
	        return new ListIterator<Item>(first);  
	    }

	    @SuppressWarnings("hiding")
		private class ListIterator<Item> implements Iterator<Item> {
	        private Node<Item> current;

	        public ListIterator(Node<Item> first) {
	            current = first;
	        }

	        public boolean hasNext()  { return current != null;                     }
	        public void remove()      { throw new UnsupportedOperationException();  }

	        public Item next() {
	            if (!hasNext()) throw new NoSuchElementException();
	            Item item = current.item;
	            current = current.next; 
	            return item;
	        }
	    }
	}
	
	/**
	 * Queue
	 */
	public class Queue<Item> implements Iterable<Item> {
	    private int N;               // number of elements on queue
	    private Node<Item> first;    // beginning of queue
	    private Node<Item> last;     // end of queue

	    @SuppressWarnings("hiding")
		private class Node<Item> {
	        private Item item;
	        private Node<Item> next;
	    }
	   
	    public boolean isEmpty() {
	        return first == null;
	    }
	    
	    public void enqueue(Item item) {
	        Node<Item> oldlast = last;
	        last = new Node<Item>();
	        last.item = item;
	        last.next = null;
	        if (isEmpty()) first = last;
	        else           oldlast.next = last;
	        N++;
	    }

	    public Item dequeue() {
	        if (isEmpty()) throw new NoSuchElementException("Queue underflow");
	        Item item = first.item;
	        first = first.next;
	        N--;
	        if (isEmpty()) last = null;   // to avoid loitering
	        return item;
	    }

	    public Iterator<Item> iterator()  {
	        return new ListIterator<Item>(first);  
	    }

	    @SuppressWarnings("hiding")
		private class ListIterator<Item> implements Iterator<Item> {
	        private Node<Item> current;

	        public ListIterator(Node<Item> first) {
	            current = first;
	        }

	        public boolean hasNext()  { return current != null;                     }
	        public void remove()      { throw new UnsupportedOperationException();  }

	        public Item next() {
	            if (!hasNext()) throw new NoSuchElementException();
	            Item item = current.item;
	            current = current.next; 
	            return item;
	        }
	    }
	}

	/**
	 * FlowEdge
	 */
	public class FlowEdge {
	    private final int v;             // from
	    private final int w;             // to 
	    private final double capacity;   // capacity
	    private double flow;             // flow

	    public FlowEdge(int v, int w, double capacity) {
	        if (v < 0) throw new IndexOutOfBoundsException("Vertex name must be a nonnegative integer");
	        if (w < 0) throw new IndexOutOfBoundsException("Vertex name must be a nonnegative integer");
	        if (!(capacity >= 0.0)) throw new IllegalArgumentException("Edge capacity must be nonnegaitve");
	        this.v         = v;
	        this.w         = w;  
	        this.capacity  = capacity;
	        this.flow      = 0.0;
	    }

	    public int from() {
	        return v;
	    }  

	    public int to() {
	        return w;
	    }  

	    public double capacity() {
	        return capacity;
	    }

	    public double flow() {
	        return flow;
	    }

	    public int other(int vertex) {
	        if      (vertex == v) return w;
	        else if (vertex == w) return v;
	        else throw new IllegalArgumentException("Illegal endpoint");
	    }

	    public double residualCapacityTo(int vertex) {
	        if      (vertex == v) return flow;              // backward edge
	        else if (vertex == w) return capacity - flow;   // forward edge
	        else throw new IllegalArgumentException("Illegal endpoint");
	    }

	    public void addResidualFlowTo(int vertex, double delta) {
	        if      (vertex == v) flow -= delta;           // backward edge
	        else if (vertex == w) flow += delta;           // forward edge
	        else throw new IllegalArgumentException("Illegal endpoint");
	        if (Double.isNaN(delta)) throw new IllegalArgumentException("Change in flow = NaN");
	        if (!(flow >= 0.0))      throw new IllegalArgumentException("Flow is negative");
	        if (!(flow <= capacity)) throw new IllegalArgumentException("Flow exceeds capacity");
	    }
	}
	
	/**
	 * Flownetwork
	 */
	public class FlowNetwork {
	    private final int V;
	    private int E;
	    private Bag<FlowEdge>[] adj;
	   
	    @SuppressWarnings("unchecked")
		public FlowNetwork(int V) {
	        if (V < 0) throw new IllegalArgumentException("Number of vertices in a Graph must be nonnegative");
	        this.V = V;
	        this.E = 0;
	        adj = (Bag<FlowEdge>[]) new Bag[V];
	        for (int v = 0; v < V; v++)
	            adj[v] = new Bag<FlowEdge>();
	    }

	    public FlowNetwork(Readfile in) {
	        this(in.readInt());
	        int E = in.readInt();
	        if (E < 0) throw new IllegalArgumentException("Number of edges must be nonnegative");
	        for (int i = 0; i < E; i++) {
	            int v = in.readInt();
	            int w = in.readInt();
	            if (v < 0 || v >= V) throw new IndexOutOfBoundsException("vertex " + v + " is not between 0 and " + (V-1));
	            if (w < 0 || w >= V) throw new IndexOutOfBoundsException("vertex " + w + " is not between 0 and " + (V-1));
	            double capacity = in.readDouble();
	            addEdge(new FlowEdge(v, w, capacity));
	        }
	    }

	    public int V() {
	        return V;
	    }
	    
	    public int E() {
	        return E;
	    }

	    private void validateVertex(int v) {
	        if (v < 0 || v >= V)
	            throw new IndexOutOfBoundsException("vertex " + v + " is not between 0 and " + (V-1));
	    }

	    public void addEdge(FlowEdge e) {
	        int v = e.from();
	        int w = e.to();
	        validateVertex(v);
	        validateVertex(w);
	        adj[v].add(e);
	        adj[w].add(e);
	        E++;
	    }

	    public Iterable<FlowEdge> adj(int v) {
	        validateVertex(v);
	        return adj[v];
	    }
	}
	
	/**
	 * FordFulkerson
	 */
	public class FordFulkerson {
	    private boolean[] marked;     // marked[v] = true iff s->v path in residual graph
	    private FlowEdge[] edgeTo;    // edgeTo[v] = last edge on shortest residual s->v path
	    private double value;         // current value of max flow

	    public FordFulkerson(FlowNetwork G, int s, int t) {
	        validate(s, G.V());
	        validate(t, G.V());
	        if (s == t)               throw new IllegalArgumentException("Source equals sink");
	        if (!isFeasible(G, s, t)) throw new IllegalArgumentException("Initial flow is infeasible");

	        // while there exists an augmenting path, use it
	        value = excess(G, t);
	        while (hasAugmentingPath(G, s, t)) {

	            // compute bottleneck capacity
	            double bottle = Double.POSITIVE_INFINITY;
	            for (int v = t; v != s; v = edgeTo[v].other(v)) {
	                bottle = Math.min(bottle, edgeTo[v].residualCapacityTo(v));
	            }

	            // augment flow
	            for (int v = t; v != s; v = edgeTo[v].other(v)) {
	                edgeTo[v].addResidualFlowTo(v, bottle); 
	            }

	            value += bottle;
	        }

	        // check optimality conditions
	        assert check(G, s, t);
	    }

	    public double value()  {
	        return value;
	    }

	    public boolean inCut(int v)  {
	        validate(v, marked.length);
	        return marked[v];
	    }

	    // throw an exception if v is outside prescibed range
	    private void validate(int v, int V)  {
	        if (v < 0 || v >= V)
	            throw new IndexOutOfBoundsException("vertex " + v + " is not between 0 and " + (V-1));
	    }

	    private boolean hasAugmentingPath(FlowNetwork G, int s, int t) {
	        edgeTo = new FlowEdge[G.V()];
	        marked = new boolean[G.V()];

	        // breadth-first search
	        Queue<Integer> queue = new Queue<Integer>();
	        queue.enqueue(s);
	        marked[s] = true;
	        while (!queue.isEmpty() && !marked[t]) {
	            int v = queue.dequeue();

	            for (FlowEdge e : G.adj(v)) {
	                int w = e.other(v);

	                // if residual capacity from v to w
	                if (e.residualCapacityTo(w) > 0) {
	                    if (!marked[w]) {
	                        edgeTo[w] = e;
	                        marked[w] = true;
	                        queue.enqueue(w);
	                    }
	                }
	            }
	        }

	        // is there an augmenting path?
	        return marked[t];
	    }

	    // return excess flow at vertex v
	    private double excess(FlowNetwork G, int v) {
	        double excess = 0.0;
	        for (FlowEdge e : G.adj(v)) {
	            if (v == e.from()) excess -= e.flow();
	            else               excess += e.flow();
	        }
	        return excess;
	    }

	    // return excess flow at vertex v
	    private boolean isFeasible(FlowNetwork G, int s, int t) {
	        double EPSILON = 1E-11;

	        // check that capacity constraints are satisfied
	        for (int v = 0; v < G.V(); v++) {
	            for (FlowEdge e : G.adj(v)) {
	                if (e.flow() < -EPSILON || e.flow() > e.capacity() + EPSILON) {
	                    System.err.println("Edge does not satisfy capacity constraints: " + e);
	                    return false;
	                }
	            }
	        }

	        // check that net flow into a vertex equals zero, except at source and sink
	        if (Math.abs(value + excess(G, s)) > EPSILON) {
	            System.err.println("Excess at source = " + excess(G, s));
	            System.err.println("Max flow         = " + value);
	            return false;
	        }
	        if (Math.abs(value - excess(G, t)) > EPSILON) {
	            System.err.println("Excess at sink   = " + excess(G, t));
	            System.err.println("Max flow         = " + value);
	            return false;
	        }
	        for (int v = 0; v < G.V(); v++) {
	            if (v == s || v == t) continue;
	            else if (Math.abs(excess(G, v)) > EPSILON) {
	                System.err.println("Net flow out of " + v + " doesn't equal zero");
	                return false;
	            }
	        }
	        return true;
	    }

	    // check optimality conditions
	    private boolean check(FlowNetwork G, int s, int t) {

	        // check that flow is feasible
	        if (!isFeasible(G, s, t)) {
	            System.err.println("Flow is infeasible");
	            return false;
	        }

	        // check that s is on the source side of min cut and that t is not on source side
	        if (!inCut(s)) {
	            System.err.println("source " + s + " is not on source side of min cut");
	            return false;
	        }
	        if (inCut(t)) {
	            System.err.println("sink " + t + " is on source side of min cut");
	            return false;
	        }

	        // check that value of min cut = value of max flow
	        double mincutValue = 0.0;
	        for (int v = 0; v < G.V(); v++) {
	            for (FlowEdge e : G.adj(v)) {
	                if ((v == e.from()) && inCut(e.from()) && !inCut(e.to()))
	                    mincutValue += e.capacity();
	            }
	        }

	        double EPSILON = 1E-11;
	        if (Math.abs(mincutValue - value) > EPSILON) {
	            System.err.println("Max flow value = " + value + ", min cut value = " + mincutValue);
	            return false;
	        }

	        return true;
	    }
	}
	
	/**
	 * Escapeplan
	 */
	public EscapePlan (Readfile file){
		this.V = file.readInt();
		this.E = file.readInt();
		this.X = file.readInt();
		this.S = file.readInt();
		FN = new FlowNetwork (V+2);
		
		for (int i=0; i < X; i++){
			int Xtemp=file.readInt();
			Xlist.add(Xtemp);
			FlowEdge e = new FlowEdge(V,Xtemp,capacity);
			FN.addEdge(e);
		}
		
		for (int i=0; i < S; i++){
			int Stemp=file.readInt();
			Slist.add(Stemp);
			FlowEdge e = new FlowEdge(Stemp,V+1,X);
			FN.addEdge(e);
		}
	
		for (int i = 0; i < E; i++){
			int No1 = file.readInt();
        	int No2 = file.readInt();
        	FlowEdge e = new FlowEdge (No1, No2, capacity);	
        	FN.addEdge(e);
		}		
	}
	
	/**
	 * maxflow
	 */
	public int maxflow(){
		int s = V;
		int t = V+1;
		int fv;
		FordFulkerson mf = new FordFulkerson(FN, s, t);
		fv = (int)mf.value();
		// StdOut.println("Maxflow: "+fv);
		if (fv == X) System.out.println("YES");
		else System.out.println("NO");
		return fv;
	}
	
	/**
	 * readfile
	 */
    public static final class Readfile {
	    
    	private static Random random;
	    private Scanner scanner;
	    private static final String CHARSET_NAME = "UTF-8";
	    private final Locale LOCALE = Locale.US;
	    
	    public Readfile(String s) {
	        try {
	            // first try to read file from local file system
	            File file = new File(s);
	            if (file.exists()) {
	                scanner = new Scanner(file, CHARSET_NAME);
	                scanner.useLocale(LOCALE);
	                return;
	            }
	            
	            URL url = getClass().getResource(s);

	            if (url == null) { url = new URL(s); }

	            URLConnection site = url.openConnection();

	            InputStream is     = site.getInputStream();
	            scanner            = new Scanner(new BufferedInputStream(is), CHARSET_NAME);
	            scanner.useLocale(LOCALE);
	        }
	        catch (IOException ioe) {
	            System.err.println("Could not open " + s);
	        }
	    } 
	    
	    public int readInt() {
	        return scanner.nextInt();
	    }
	    
	    public double readDouble() {
	        return scanner.nextDouble();
	    }
	    
	    public static int uniform(int N) {
	        if (N <= 0) throw new IllegalArgumentException("Parameter N must be positive");
	        return random.nextInt(N);
	    }

	}

	public static void main(String[] args) {
		Readfile file1 = new Readfile("input.txt");
		EscapePlan escape1 = new EscapePlan (file1);
		escape1.maxflow();
		/*Readfile file2 = new Readfile("input2.txt");
		EscapePlan escape2 = new EscapePlan (file2);
		escape2.maxflow();
		Readfile file3 = new Readfile("input3.txt");
		EscapePlan escape3 = new EscapePlan (file3);
		escape3.maxflow();*/
	}

}

