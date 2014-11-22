public class Node
{
	private double priority;     // data item (key)
	private int identity;	// data item (place in array)
	private int parent; 		// data item (place of parent in array)
	// -------------------------------------------------------------
	public Node(double p, int ID, int P)           // constructor
	{ priority = p;
	identity = ID;
	parent = P;}
	// -------------------------------------------------------------
	public double getPriority()
	{
		return priority;
	}
	// -------------------------------------------------------------
	public void setPriority(double pr)
	{
		priority = pr;
	}
	// -------------------------------------------------------------
	public int getParent()
	{
		return parent;
	}
	// -------------------------------------------------------------
	public void setParent(int pa)
	{
		parent = pa;
	}
	// -------------------------------------------------------------
	public int getIdentity()
	{
		return identity;
	}
	// -------------------------------------------------------------
	public void setIdentity(int id)
	{
		identity = id;
	}
}  // end class Node
////////////////////////////////////////////////////////////////