public class VertexHeap
{
	public Node[] heapArray;
	private int maxSize;           // size of array
	private int currentSize;       // number of nodes in array
	public int[] map; 		  // holds a map of the vertices to locations in the heap
	// -------------------------------------------------------------
	public VertexHeap(int mx)            // constructor
	{
		maxSize = mx;
		currentSize = 0;

		// create the heap array
		heapArray = new Node[maxSize];  

		// Create and initialize the map array
		map = new int[maxSize];
		// Sriram: reduced this also to maxSize
		for(int i=0; i< maxSize;i++)
			map[i] = -1;
	}
	// -------------------------------------------------------------
	public boolean isEmpty()
	{ return currentSize==0; }
	// -------------------------------------------------------------
	// A method to print the map array; is useful for testing. It only prints
	// elements not -1.
	public void printMap()
	{
		for(int i = 0; i < maxSize; i++)
			if(map[i] != -1)
				System.out.print("("  + i + "," +  map[i] + "," + heapArray[map[i]].getPriority() + ") ");
		System.out.println(" ");
	} 
	// -------------------------------------------------------------
	public boolean insert(Node newNode)
	{
		if(currentSize==maxSize)
			return false;

		heapArray[currentSize] = newNode;
		map[newNode.getIdentity()] = currentSize;
		trickleUp(currentSize++);
		return true;
	}  // end insert()
	// -------------------------------------------------------------
	public void trickleUp(int index)
	{
		int parent = (index-1) / 2;
		Node bottom = heapArray[index];

		// Main loop for trickling up
		while( index > 0 && heapArray[parent].getPriority() > bottom.getPriority() )
		{
			heapArray[index] = heapArray[parent];  // move it down
			map[heapArray[parent].getIdentity()] = index;
			index = parent;        
			parent = (parent-1) / 2;
		}

		heapArray[index] = bottom;
		map[bottom.getIdentity()] = index;
	}  // end trickleUp()
	// -------------------------------------------------------------
	public Node delete()           // delete item with max key
	{  						// (assumes non-empty list)
		if(!isEmpty())
		{
			Node root = heapArray[0];
			map[root.getIdentity()] = -1;
			heapArray[0] = heapArray[--currentSize];
			trickleDown(0);
			return root;
		}
		return null;
	}  // end remove()
	// -------------------------------------------------------------
	public void trickleDown(int index)
	{
		int smallerChild;
		Node top = heapArray[index];       // save root

		while(index < currentSize/2)       // while node has at least one child
		{                               
			int leftChild = 2*index+1;
			int rightChild = leftChild+1;

			// find larger child
			if(rightChild < currentSize &&  heapArray[leftChild].getPriority() > heapArray[rightChild].getPriority())
				smallerChild = rightChild;
			else
				smallerChild = leftChild;

			// top >= largerChild?
			if( top.getPriority() <= heapArray[smallerChild].getPriority() )
				break;

			// shift child up
			heapArray[index] = heapArray[smallerChild];
			map[heapArray[smallerChild].getIdentity()] = index;
			index = smallerChild;            // go down
		}  // end while

		heapArray[index] = top;            // root to index
		map[top.getIdentity()] = index;

	}  // end trickleDown()
	// -------------------------------------------------------------
	public boolean change(int index, double newPriority, int newParent)
	{
		if(index<0 || index>=currentSize)
			return false;

		double oldValue = heapArray[index].getPriority(); // remember old
		heapArray[index].setPriority(newPriority);  // change to new
//		heapArray[index].setParent(newParent);  // change to new

		if(oldValue > newPriority)             // if raised,
			trickleUp(index);                // trickle it up
		else                                // if lowered,
			trickleDown(index);              // trickle it down
		return true;
	}  // end change()
	// -------------------------------------------------------------
	public void displayHeap()
	{
		System.out.print("heapArray: ");    // array format
		for(int m=0; m<currentSize; m++)
			if(heapArray[m] != null)
				System.out.print( heapArray[m].getPriority() + " ");
			else
				System.out.print( "-- ");
		System.out.println();
		// heap format
		int nBlanks = 32;
		int itemsPerRow = 1;
		int column = 0;
		int j = 0;                          // current item
		String dots = "...............................";
		System.out.println(dots+dots);      // dotted top line

		while(currentSize > 0)              // for each heap item
		{
			if(column == 0)                  // first item in row?
				for(int k=0; k<nBlanks; k++)  // preceding blanks
					System.out.print(' ');
			// display item
			System.out.print(heapArray[j].getPriority());

			if(++j == currentSize)           // done?
				break;

			if(++column==itemsPerRow)        // end of row?
			{
				nBlanks /= 2;                 // half the blanks
				itemsPerRow *= 2;             // twice the items
				column = 0;                   // start over on
				System.out.println();         //    new row
			}
			else                             // next item on row
				for(int k=0; k<nBlanks*2-2; k++)
					System.out.print(' ');     // interim blanks
		}  // end for
		System.out.println("\n"+dots+dots); // dotted bottom line
	}  // end displayHeap()
	// -------------------------------------------------------------
	public int getIndex(int vertex)
	{
		int placeInTree = map[vertex];

		return placeInTree;
	}
	// -------------------------------------------------------------
	public double getPriority(int index)
	{
		double x = heapArray[index].getPriority(); 
		return x;
	}
	// -------------------------------------------------------------
	public void clearHeap(VertexHeap heap)
	{
		while(!heap.isEmpty())
		{
			heap.delete();
		}
	}
	// -------------------------------------------------------------
}  // end class VertexHeap