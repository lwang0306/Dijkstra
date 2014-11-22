public class EdgeLink
{
	public String sData;   	// name of neighbor
	public double weight; 	// weight of edge connecting to neighbor
	public EdgeLink next;   // next link in list
	// -------------------------------------------------------------
	public EdgeLink(String id, double w,  EdgeLink newNext) // constructor
	{
		sData = id;
		weight = w;
		next = newNext;
	}
	// -------------------------------------------------------------
	public void displayLink()      // display ourself
	{
		System.out.print("{" + sData + " " + weight + "} ");
	}
}  // end class EdgeLink