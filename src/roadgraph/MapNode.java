package roadgraph;

import java.util.Set;
import java.util.HashMap;
import geography.GeographicPoint;


public class MapNode {
	HashMap<GeographicPoint,MapEdge> edges;
	
	public MapNode() 
	{
		edges = new HashMap<GeographicPoint, MapEdge>();
	}
	
	public HashMap<GeographicPoint, MapEdge> getEdges()
	{
		return edges;
	}
	
	public Set<GeographicPoint> getNeighborPoints() 
	{
		return edges.keySet();
	}
	
	public MapEdge getEdge(GeographicPoint to)
	{
		if (to == null) {
			return null;
		} else if (!edges.containsKey(to)) {
			return null;
		} else {
			return edges.get(to);
		}
	}
	
	public void addEdge(GeographicPoint to, String roadName, String roadType, double length)
	{
		edges.put(to,new MapEdge(roadName,roadType,length));
	}

}
