package roadgraph;

import geography.GeographicPoint;
import java.util.List;
import java.util.Set;
import java.util.HashMap;
import java.util.HashSet;

public class TSP_Node {
	
	private GeographicPoint point;
	
	private HashMap<GeographicPoint,Double> distanceTo;
	private HashMap<GeographicPoint, List<GeographicPoint>> pathMap;
	
	public TSP_Node(GeographicPoint point)
	{
		this.point = point;
		distanceTo = new HashMap<GeographicPoint,Double>();
		pathMap = new  HashMap<GeographicPoint, List<GeographicPoint>>();
	}

	public GeographicPoint getPoint() 
	{
		return point;
	}
	
	public List<GeographicPoint> getPath(GeographicPoint to)
	{
		if (pathMap.containsKey(to)) {
			return pathMap.get(to);
		}
		return null;
	}
	
	public HashMap<GeographicPoint,Double> getPointDistanceMap() 
	{
		return distanceTo;
	}
	
	public Double getDistToPoint(GeographicPoint to)
	{
		if (distanceTo.containsKey(to)) {
			return distanceTo.get(to);
		}
		return(Double.POSITIVE_INFINITY);
}
	
	public void addPath(GeographicPoint to, List<GeographicPoint> path)
	{
		pathMap.put(to, path);
	}
	
	public void setDistanceToPoint(GeographicPoint to, Double distance) 
	{
		distanceTo.put(to, distance);
	}
	
	public GeographicPoint findMinPath(Set<GeographicPoint> visited) 
	{
		double distance = Double.POSITIVE_INFINITY;
		GeographicPoint minPoint = null;
		for (GeographicPoint next : pathMap.keySet()) {
			if (!visited.contains(next)) {
				double distPt = getDistToPoint(next);
				if (distPt < distance) {
					distance = distPt;
					minPoint = next;
				}
			}
		}		
		return minPoint;
	}
}
