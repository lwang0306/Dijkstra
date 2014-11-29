import java.lang.Math;
import java.util.Hashtable;

public class WeightedGraph {
	private String[] names; // 1-d array to store the labels of each vertex
	private Hashtable<String, Integer> nameIndex;
	private double[] x; // 1-d array to store x-coordinate of each vertex
	private double[] y; // 1-d array to store y-coordinate of each vertex
	private EdgeLinkList[] Edges; // 1-d array to store adjacencies between
									// vertices
	private int numVertices;
	private int numEdges;
	private Hashtable<String, String[]> computedPath;
	private Hashtable<String, Double> computedCost;

	// Default constructor. Sets aside enough capacity for one vertex
	public WeightedGraph() {
		this(1);
	}

	// Constructor that sets aside as much capacity as specified by the user
	public WeightedGraph(int capacity) {
		names = new String[capacity];
		nameIndex = new Hashtable<String, Integer>();
		x = new double[capacity];
		y = new double[capacity];
		Edges = new EdgeLinkList[capacity];
		for (int i = 0; i < capacity; i++) {
			Edges[i] = new EdgeLinkList();
		}
		computedPath = new Hashtable<String, String[]>();
		computedCost = new Hashtable<String, Double>();

	}

	public int numberOfVertices() {
		return numVertices;
	}

	public int numberOfEdges() {
		return numEdges;
	}

	// Finds the location at which a vertex is stored in Vertices.
	// Returns -1 if vertex not found
	public int getIndex(String vertex) {
		for (int i = 0; i < numVertices; i++)
			if (vertex.equals(names[i]))
				return i;

		return -1;
	}

	// Resizes the array of vertices. Can make it larger or smaller, depending
	// on what newSize is.
	private String[] resize(String[] array, int newSize) {
		String[] temp = new String[newSize];

		int smallerSize = newSize;
		if (array.length < smallerSize)
			smallerSize = array.length;

		for (int i = 0; i < smallerSize; i++)
			temp[i] = array[i];

		return temp;
	}

	// Resizes the array of coordinates. Can make it larger or smaller,
	// depending
	// on what newSize is.
	private double[] resize(double[] array, int newSize) {
		double[] temp = new double[newSize];

		int smallerSize = newSize;
		if (array.length < smallerSize)
			smallerSize = array.length;

		for (int i = 0; i < smallerSize; i++)
			temp[i] = array[i];

		return temp;
	}

	// Resizes the array of Edges. Can make it larger or smaller, depending
	// on what newSize is.
	private EdgeLinkList[] resize(EdgeLinkList[] array, int newSize) {
		EdgeLinkList[] temp = new EdgeLinkList[newSize];

		int smallerSize = newSize;
		if (array.length < smallerSize)
			smallerSize = array.length;

		for (int i = 0; i < smallerSize; i++)
			temp[i] = array[i];

		for (int i = smallerSize; i < temp.length; i++)
			temp[i] = new EdgeLinkList();

		return temp;
	}

	// Adds a new vertex
	public void addVertex(String newVertex, double xcoord, double ycoord) {
		if (nameIndex.containsKey(newVertex)) {
			System.out.print("addVertex: ");
			System.out.print(newVertex);
			System.out.println(" failed -- vertex already exists.");
			return;
		}

		// if array of vertices is full, we need to expand it and
		// also expand Edges
		if (names.length == numVertices) {
			names = resize(names, 2 * numVertices + 1);
			x = resize(x, 2 * numVertices + 1);
			y = resize(y, 2 * numVertices + 1);
			Edges = resize(Edges, 2 * numVertices + 1);
		}

		names[numVertices] = newVertex;
		nameIndex.put(newVertex, numVertices);
		x[numVertices] = xcoord;
		y[numVertices++] = ycoord;
	}

	// Adds a new weighted edge. The edge is specified by
	// indices of endpoints and the weight equals the Euclidean distance.
	public void addWeightedEdge(String vertex1, String vertex2) {
		int i = nameIndex.get(vertex1);
		int j = nameIndex.get(vertex2);

		if ((i < 0) || (i > numVertices - 1)) {
			System.out.print("addEdge failed: ");
			System.out.print("index " + i);
			System.out.println(" is out of bounds.");
			return;
		}

		if ((j < 0) || (j > numVertices - 1)) {
			System.out.print("addEdge failed: ");
			System.out.print("index " + j);
			System.out.println(" is out of bounds.");
			return;
		}

		double x1 = x[i];
		double y1 = y[i];
		double x2 = x[j];
		double y2 = y[j];
		double diffX = x1 - x2;
		double diffY = y1 - y2;
		double dist = Math.sqrt(diffX * diffX + diffY * diffY);
		addEdge(vertex1, vertex2, dist);
	}

	// Adds a new edge with weight w. The edge is specified by
	// indices of endpoints
	public void addEdge(int i, int j, double w) {
		if ((i < 0) || (i > numVertices - 1)) {
			System.out.print("addEdge failed: ");
			System.out.print("index " + i);
			System.out.println(" is out of bounds.");
			return;
		}

		if ((j < 0) || (j > numVertices - 1)) {
			System.out.print("addEdge failed: ");
			System.out.print("index " + j);
			System.out.println(" is out of bounds.");
			return;
		}

		Edges[i].insertFirst(names[j], w);
		Edges[j].insertFirst(names[i], w);

		numEdges++;
	}

	// Adds a new edge with weight w
	public void addEdge(String vertex1, String vertex2, double w) {
		int i = nameIndex.get(vertex1);
		int j = nameIndex.get(vertex2);

		if (i == -1) {
			System.out.print("addEdge failed: ");
			System.out.print(vertex1);
			System.out.println(" does not exist.");
			return;
		}

		if (j == -1) {
			System.out.print("addEdge failed: ");
			System.out.print(vertex2);
			System.out.println(" does not exist.");
			return;
		}

		addEdge(i, j, w);
	}

	// returns the names of all the neighbors of a given vertex in a
	// String array
	private String[] getNeighbors(String vertex) {
		int source = nameIndex.get(vertex);
		if (source == -1) {
			System.out.print("getNeighbors failed: Vertex ");
			System.out.print(vertex);
			System.out.println(" does not exist.");
			return null;
		}

		return Edges[source].copyIntoArray();
	}

	// returns the indices of all the neighbors of a given vertex. The
	// vertex is specified as an index and the neighbors are returned
	// in an int array
	private int[] getNeighbors(int index) {
		if ((index < 0) || (index >= numVertices)) {
			System.out.print("getNeighbors failed: Index");
			System.out.print(index);
			System.out.println(" is out of bounds.");
			return null;
		}

		// Call the earlier getNeighbors function to get the names of
		// neighbors
		String[] nbrNames = getNeighbors(names[index]);

		// Turn the array of neighbor names into an array
		// of neighbor indices
		int[] nbrIndices = new int[nbrNames.length];
		for (int i = 0; i < nbrIndices.length; i++)
			nbrIndices[i] = nameIndex.get(nbrNames[i]);

		return nbrIndices;
	}

	// Gets the weight of the edge connecting a pair of vertices
	// whose indices are given
	private Double getWeight(int i, int j) {
		if ((i < 0) || (i > numVertices - 1)) {
			System.out.print("getWeight failed: ");
			System.out.print("index " + i);
			System.out.println(" out of bounds.");
			return null;
		}

		if ((j < 0) || (j > numVertices - 1)) {
			System.out.print("getWeight failed: ");
			System.out.print("index " + j);
			System.out.println(" out of bounds.");
			return null;
		}

		// Look for vertex j in Edges[i]
		EdgeLink e = Edges[i].find(names[j]);

		// If vertex j is found in Edges[i] then return the weight of
		// the edge, otherwise return null
		if (e != null)
			return new Double(e.weight);
		else
			return null;
	}

	/*
	 * Calculates the shortest path from source to dest using Dijkstra's
	 * Shortest Path. The return value contains the labels of the vertices that
	 * are visited along the path.
	 */
	public String[] shortestPath(String source, String dest) {
		// Get index of source
		int sourceIndex = nameIndex.get(source);
		// Get index of destination
		int destIndex = nameIndex.get(dest);
		// form a key for (sourceIndex, destIndex)
		String key = "(" + sourceIndex + ", " + destIndex + ")";

		if (computedPath.containsKey(key))
			return computedPath.get(key);
		else {
			if (sourceIndex == -1) {
				System.out.print("shortestPath failed: ");
				System.out.print(source);
				System.out.println(" does not exist.");
				return null;
			}
			if (destIndex == -1) {
				System.out.print("shortestPath failed: ");
				System.out.print(dest);
				System.out.println(" does not exist.");
				return null;
			}

			// Perform DSP from destination
			int[] spTree = DSP(dest, source);

			// If source is unreachable from destination
			if (spTree[sourceIndex] == -1)
				return null;

			// Define a String[] for shortest path and place the source vertex
			// in it
			String[] path = new String[numVertices];
			path[0] = names[sourceIndex];

			// Start following parent pointers and store each new vertex
			// encountered, in the path array. The while-loop executes
			// until the root of the tree is encountered
			int currentIndex = sourceIndex;
			int pathLength = 0;
			// int spCurrent = spTree[currentIndex];
			// nameIndex.put(names[spCurrent], currentIndex);
			while (currentIndex != spTree[currentIndex]) {
				currentIndex = spTree[currentIndex];
				if (names[currentIndex] == null)
					break;
				// System.out.println("-------- currentIndex = " + currentIndex
				// + " -------------");
				pathLength++;
				path[pathLength] = names[currentIndex];
			}

			// Resize the path array to be exactly of the correct size
			String[] newPath = new String[pathLength + 1];
			for (int i = 0; i < newPath.length; i++)
				newPath[i] = path[i];

			// save the result of existing query into Hashtable
			computedPath.put(key, newPath);

			return newPath;
		}
	}

	/*
	 * Calculates the cost of the shortest path from source to dest.
	 */
	public double shortestPathCost(String source, String dest) {
		// Get index of source
		int sourceIndex = nameIndex.get(source);
		// Get index of destination
		int destIndex = nameIndex.get(dest);
		// form a key for (sourceIndex, destIndex)
		String key = "(" + sourceIndex + ", " + destIndex + ")";
		String[] path;

		if (computedCost.containsKey(key))
			return computedCost.get(key);
		else if (computedPath.containsKey(key))
			path = computedPath.get(key);
		else {
			if (sourceIndex == -1) {
				System.out.print("shortestPathCost failed: ");
				System.out.print(source);
				System.out.println(" does not exist.");
				return Double.MAX_VALUE;
			}

			if (destIndex == -1) {
				System.out.print("shortestPathCost failed: ");
				System.out.print(dest);
				System.out.println(" does not exist.");
				return Double.MAX_VALUE;
			}
			path = shortestPath(source, dest);
			computedPath.put(key, path);
		}
		double pathCost = 0;
		for (int i = 0; i < path.length - 1; i++) {
			int uIndex = nameIndex.get(path[i]);
			int vIndex = nameIndex.get(path[i + 1]);
			pathCost += getWeight(uIndex, vIndex);
		}
		computedCost.put(key, new Double(pathCost));
		return pathCost;
	}

	/*
	 * Implementation of Dijkstra's shortest path algorithm. The return value is
	 * a representation of the shortest path tree, where the k-th value of the
	 * array gives the index of the vertex that is the parent of vertex #k in
	 * the tree with the source as the root. That is, the k-th value tells which
	 * vertex number to visit next when traversing from vertex #k back to the
	 * source.
	 */

	private int[] DSP(String source, String dest) {

		int sourceIndex = nameIndex.get(source);

		// Declarations
		double[] dist = new double[numVertices];
		int[] previous = new int[numVertices];
		VertexHeap Q = new VertexHeap(numVertices);

		// Initializations
		for (int i = 0; i < numVertices; i++) // Initializations
		{
			dist[i] = Double.MAX_VALUE; // Unknown distance function from s to v
			previous[i] = -1; // the array previous stores parent info.

			// Create a Vertex object corresponding to vertex i and insert into
			// VertexHeap
			Node v = new Node(dist[i], i, previous[i]);
			Q.insert(v);
		}

		dist[sourceIndex] = 0; // Distance from source to source is set to 0
		previous[sourceIndex] = sourceIndex; // parent of sourceIndex is itself
		Q.change(Q.getIndex(sourceIndex), 0, sourceIndex); // and this is
															// updated in the
															// priority queue as
															// well

		// The main loop
		while (!Q.isEmpty()) {
			Node u = Q.delete(); // Remove best vertex from priority queue;
									// returns source on first iteration
			int uIndex = u.getIdentity();

			// get the neighbors of u
			int[] nbrs = getNeighbors(uIndex);
			// if (nbrs.length == 2) {
			// Edges[nbrs[0]].delete(names[uIndex]);
			// Edges[nbrs[1]].delete(names[uIndex]);
			// }

			for (int j = 0; j < nbrs.length; j++) {
				int vIndex = nbrs[j];
				int heapVIndex = Q.getIndex(vIndex);
				if (heapVIndex != -1) { // this line is added to eliminate
										// unnecessary work
					double alt = dist[uIndex] + getWeight(uIndex, vIndex);
					if (alt < dist[vIndex]) // Relax (u,v)
					{
						dist[vIndex] = alt;
						previous[vIndex] = uIndex;
						Q.change(heapVIndex, dist[vIndex], uIndex);
					} // end of if alt < dist[vIndex]
				} // end of if heapVIndex != -1
			} // end of for-loop that scans the neighbors
			if (nameIndex.get(dest) == uIndex)
				break;
		} // end of while-Q-is-not-empty

		return previous;

	} // end of function

	public void replaceTwoEdgesWithOneEdge() {

	}

} // end of class