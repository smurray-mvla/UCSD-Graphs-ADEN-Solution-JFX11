package roadgraph;

import geography.GeographicPoint;

/**
 * @author Scott
 * 
 * A class to represent a segment (edge) of road that connects two vertices in a graph. 
 */
public class MapEdge {

	//private GeographicPoint to;
	private String roadName;
	private String roadType;
	private double roadLength;
	
	/**
	 * @param to       - the endpoint of this segment of the road
	 * @param name     - the name of the road
	 * @param type	   - the type of the road
	 * @param length   - the length of the road
	 */
	public MapEdge(/*GeographicPoint to,*/ String name, String type, double length) {
		//this.to = to;
		this.roadName = name;
		this.roadType = type;
		this.roadLength = length;
	}
	
	/* getters */ 
	
	/**
	 * @return The GeographicPoint at the end of the road.
	 */
//	public GeographicPoint getTo() {
//		return to;
//	}
	

	/**
	 * @return The name of the road
	 */
	public String getRoadName() {
		return roadName;
	}

	/**
	 * @return The type of road
	 */
	public String getRoadType() {
		return roadType;
	}

	/**
	 * @return The length of the road 
	 */
	public double getRoadLength() {
		return roadLength;
	}
	
}
