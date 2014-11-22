import java.io.File;
import java.util.Scanner;

/*********************************************
 * 
 * Do not change this file!!
 * 
 *********************************************/


public class Main {

	public static void main(String[] args) {

		// make sure all the runtime arguments are properly specified
		if (args.length < 4) {
			System.out.println("Usage: java Main [name_of_data_file] [label_of_start_vertex] [label_of_dest_vertex_1] [label_of_dest_vertex_2]");
			System.exit(0);
		}
		String fileName = args[0];
		String startLabel = args[1];
		String destLabel1 = args[2];
		String destLabel2 = args[3];
		
		// build the graph
		long start = System.currentTimeMillis();
		WeightedGraph g = buildGraph(fileName);
		long end = System.currentTimeMillis();
		System.out.println("Time to build graph: " + (end-start) + "ms");
		if (g == null) {
			System.out.println("Uh oh, graph couldn't be built. I give up");
			System.exit(0);
		}
		
		// calculate the shortest path from the start to the first destination
		start = System.currentTimeMillis();
		String[] path = g.shortestPath(startLabel, destLabel1);
		end = System.currentTimeMillis();
		System.out.println("Time to find shortest path (first time) from " + startLabel + " to " + destLabel1 + ": " + (end-start) + "ms");
		
		// use this to make sure the output is still correct
		for (String edge : path) {
			System.out.print(edge + "-");
		}
		System.out.println();

		// calculate the shortest path from the start to the first destination -- again!
		start = System.currentTimeMillis();
		path = g.shortestPath(startLabel, destLabel1);
		end = System.currentTimeMillis();
		System.out.println("Time to find shortest path (second time) from " + startLabel + " to " + destLabel1 + ": " + (end-start) + "ms");

		// calculate the cost of the shortest path from the start to the first destination
		start = System.currentTimeMillis();
		double cost = g.shortestPathCost(startLabel, destLabel1);
		end = System.currentTimeMillis();
		System.out.println("Time to find cost of shortest path from " + startLabel + " to " + destLabel1 + ": " + (end-start) + "ms");
		System.out.println("Cost is " + cost);

		// calculate the shortest path from the start to the SECOND destination
		start = System.currentTimeMillis();
		path = g.shortestPath(startLabel, destLabel2);
		end = System.currentTimeMillis();
		System.out.println("Time to find shortest path from " + startLabel + " to " + destLabel2 + ": " + (end-start) + "ms");

		// use this to make sure the output is still correct
		for (String edge : path) {
			System.out.print(edge + "-");
		}
		System.out.println();
		
	}
	
	public static WeightedGraph buildGraph(String fileName) {
		
		WeightedGraph g = null;
		
		try {
			Scanner in = new Scanner(new File(fileName));
			
			// first token is number of vertices
			int numVertices = in.nextInt();
			System.out.println("Number of vertices = " + numVertices);
			// then the number of edges
			int numEdges = in.nextInt();
			System.out.println("Number of edges = " + numEdges);
			// consume the next line
			in.nextLine();

			// create the graph with the capacity equal to the number of vertices
			g = new WeightedGraph(numVertices);

			// read the vertices
			for (int i = 0; i < numVertices; i++) {
				if (in.hasNext() == false) {
					// uh-oh, someone lied to us
					throw new Exception("Error reading input file: seems like number of vertices is wrong");
				}
				String label = in.next();
				int x = in.nextInt();
				int y = in.nextInt();
				g.addVertex(label, x, y);
				//System.out.println("Added vertex " + label + " at (" + x + ", " + y + ")");
			}

			// read the edges
			for (int i = 0; i < numEdges; i++) {
				if (in.hasNext() == false) {
					// uh-oh, someone lied to us
					throw new Exception("Error reading input file: seems like number of edges is wrong");
				}
				String a = in.next();
				String b = in.next();
				g.addWeightedEdge(a, b);
				//System.out.println("Added edge from " + a + " to " + b);
			}
			
			// consumes the rest of the file
			while (in.hasNext()) in.nextLine();

		}
			
		catch (Exception e) {
			e.printStackTrace();
		}
		
		return g;
		
	}
}
