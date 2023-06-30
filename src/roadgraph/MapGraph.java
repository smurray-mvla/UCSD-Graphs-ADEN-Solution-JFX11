/**
 * @author UCSD MOOC development team and YOU
 * 
 * A class which reprsents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
package roadgraph;


import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.HashMap;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Queue;
import java.util.PriorityQueue;
import java.util.Comparator;


import geography.GeographicPoint;
import util.GraphLoader;
import week3example.MazeNode;

/**
 * @author UCSD MOOC development team and YOU
 * 
 * A class which represents a graph of geographic locations
 * Nodes in the graph are intersections between segments
 *
 */
public class MapGraph {
	private int numVertices;
	private int numEdges;
//	public int numVisited;
	private HashMap<GeographicPoint,MapNode> vertices;
	
	/** 
	 * Create a new empty MapGraph 
	 */
	public MapGraph()
	{
		numVertices = 0;
		numEdges = 0;
		// vertices will contain all vertices that exist in the graph
		vertices = new HashMap<GeographicPoint,MapNode>();
		// nodeMap maps all edges (ie, roads) and their related information that connect OUT from each vertice
	}
	
	/**
	 * Get the number of vertices (road intersections) in the graph
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices()
	{
		return numVertices;
	}
	
	/**
	 * Return the intersections, which are the vertices in this graph.
	 * @return The vertices in this graph as GeographicPoints
	 */
	public Set<GeographicPoint> getVertices()
	{
		return vertices.keySet();
	}
	
	/**
	 * Get the number of road segments in the graph
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges()
	{
		return numEdges;
	}

	
	
	/** Add a node corresponding to an intersection at a Geographic Point
	 * If the location is already in the graph or null, this method does 
	 * not change the graph.
	 * @param location  The location of the intersection
	 * @return true if a node was added, false if it was not (the node
	 * was already in the graph, or the parameter is null).
	 */
	public boolean addVertex(GeographicPoint location)
	{
		if (location == null) {
			return false;
		} else if (vertices.containsKey(location)) {
			return false;
		}
		
		/*
		System.out.printf("Adding vertice %3d (%f,%f)\n",numVertices,location.getX(),location.getY());
		*/
		
		MapNode mapNode = new MapNode();
		vertices.put(location,mapNode);
		numVertices++;
		
		return true;
	}
	
	/**
	 * Adds a directed edge to the graph from pt1 to pt2.  
	 * Precondition: Both GeographicPoints have already been added to the graph
	 * @param from The starting point of the edge
	 * @param to The ending point of the edge
	 * @param roadName The name of the road
	 * @param roadType The type of the road
	 * @param length The length of the road, in km
	 * @throws IllegalArgumentException If the points have not already been
	 *   added as nodes to the graph, if any of the arguments is null,
	 *   or if the length is less than 0.
	 */
	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName,
			String roadType, double length) throws IllegalArgumentException {

		// Error checking the input..
        if ((from == null) || (!vertices.containsKey(from)))
        	throw new IllegalArgumentException("\'from\' node cannot be null and must exist in graph");
        if ((to == null) || (!vertices.containsKey(to)))
        	throw new IllegalArgumentException("\'to\' node cannot be null and must exist in graph");
        if (roadName == null)
        	throw new IllegalArgumentException("roadName cannot be null");
        if (roadType == null)
        	throw new IllegalArgumentException("roadType cannot be null");
        if (length < 0) 
        	throw new IllegalArgumentException("length cannot be negative");

        /* for debug only
        System.out.printf("Adding Edge %3d from (%f,%f) to ((%f,%f) Name: %15s Type: %15s Length: %f\n",
        		numEdges,from.getX(),from.getY(),to.getX(),to.getY(),roadName,roadType,length);
        */
        
        // Since all parameters are valid, create a new segment and add it to the HashMap for the "from" 
		vertices.get(from).addEdge(to,roadName,roadType,length);
		numEdges ++;
	}
	

	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		//pointSearched = new ArrayList<GeographicPoint>();
        //Consumer<GeographicPoint> temp = (x) -> {pointSearched.add(x);};
        Consumer<GeographicPoint> temp = (x) -> {};
        return bfs(start, goal, temp);
	}
	
	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal). Return null if either start or goal locations
	 *   do not exist, or no path was found that connected start to goal
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, 
			 					     GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		if ((start==null) || (goal==null)) {
			System.out.println("Either start or goal points are null -- no route exists.");
			return null;   
		} 
		HashMap<GeographicPoint,GeographicPoint> parentMap = new HashMap<GeographicPoint,GeographicPoint>();
		
		boolean found = bfsSearch(start, goal, parentMap, nodeSearched);
		
		if (!found) {
			System.out.println("No path exists");
			return null;
		}

		return reconstructPath(start, goal, parentMap);

	}
	
	/** Reconstructs the path from start to goal and returns as a list of intersections in order.
	 * 
	 * @param start
	 * @param goal
	 * @param parentMap
	 * @return  path - a list of points of the shortest path in the order traversed from start to finish
	 */
	private List<GeographicPoint> reconstructPath(GeographicPoint start, GeographicPoint goal, HashMap<GeographicPoint, GeographicPoint> parentMap)
	{
		LinkedList<GeographicPoint> path = new LinkedList<GeographicPoint>();
		GeographicPoint curr = goal;
		while (curr != start) {
			path.addFirst(curr);
			curr = parentMap.get(curr);
		}
		path.addFirst(start);
		return path;

	}
	
	/** Execute a Breadth-First search attempting to find a path from start to goal. 
	 *  Returns true if a path was found or false if not. If a path is found, parentMap
	 *  will reflect the connectivity of the path
	 *  
	 * @param start - Starting location of the search
	 * @param goal  - Desired end of the path
	 * @param parentMap - Map that reflects the connectivity of the path from start to goal (if one exists)
	 * @param nodeSearched - parameter used for visualization of the search (MapApp)
	 * @return found - true if path found, false otherwise.
	 */
	private boolean bfsSearch(GeographicPoint start, GeographicPoint goal, 
			HashMap<GeographicPoint,GeographicPoint> parentMap, Consumer<GeographicPoint> nodeSearched) 
	{
		HashSet<GeographicPoint> visited = new HashSet<GeographicPoint>();
		Queue<GeographicPoint> toExplore = new LinkedList<GeographicPoint>();
		System.out.println("\nStarting BFS Search from " + start + " to "+goal+":");
		toExplore.add(start);
		visited.add(start);
///		numVisited++;
		boolean found = false;
		while (!toExplore.isEmpty()) {
			GeographicPoint curr = toExplore.remove();
			// Hook for visualization.  See writeup.
			nodeSearched.accept(curr);  // Note that this is only listing those nodes that actually *DO* get explored!

			/* debug
			System.out.println("Exploring Vertice @ ("+curr.getX()+","+curr.getY()+")");
			*/
			
			if (curr.equals(goal)) {
				found = true;
				break;
			}
			Set<GeographicPoint> neighbors = vertices.get(curr).getNeighborPoints();
			if (neighbors != null ) {
				for (GeographicPoint next : neighbors) {
					if (!visited.contains(next)) {
						visited.add(next);
	//					numVisited++;
					    parentMap.put(next,curr);
					    toExplore.add(next);
					}

				}
			}
		}
		return found;
	}
	

	/** Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		// You do not need to change this method.
        Consumer<GeographicPoint> temp = (x) -> {};
        return dijkstra(start, goal, temp);
	}
	
	/** Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, 
										  GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		if ((start==null) || (goal==null)) {
			System.out.println("Either start or goal points are null -- no route exists.");
			return null;   
		} 
		HashMap<GeographicPoint,GeographicPoint> parentMap = new HashMap<GeographicPoint,GeographicPoint>();

		boolean found = dijkstraSearch(start, goal, parentMap, nodeSearched);

		if (!found) {
			System.out.println("No path exists");
			return null;
		}

		return reconstructPath(start, goal, parentMap);
		
	}

	/** Print detailed information about an element encountered during a Dijkstra or A* search
	 * 
	 * @param type - string indicating the type of search - for information only
	 * @param elem - an actual QueueElement....
	 */
	private void PrintQueueElement (String type, QueueElement elem )
	{
		
	    GeographicPoint curr = elem.getPoint();
	    Set<GeographicPoint> neighbors = vertices.get(curr).getNeighborPoints();
		System.out.println("Search: " + type + "   Node@ Lat: " + elem.getPoint().getX()+", "+elem.getPoint().getY());
		System.out.print("        Intersections: ");
		
		for (GeographicPoint to : neighbors) {
			MapEdge neighbor = vertices.get(curr).getEdge(to);
			System.out.print(neighbor.getRoadName()+", ");
		}
		System.out.println();
		System.out.println("        DistToNode: " + elem.getDistFromStart() + "  Priority: " + elem.priorityDistance());
	}
	
	/** The actual search engine implementing the Dijkstra search. Note that handling neighbors is implemented as a helper method
	 *  for readability purposes
	 *  
	 * @param start - the starting point
	 * @param goal - the desired destination
	 * @param parentMap - Map that reflects the connectivity of the path from start to goal (if one exists)
	 * @param nodeSearched - hook for visualization
	 * @return boolean indicating if the search found a path
	 */
	private boolean dijkstraSearch(GeographicPoint start, GeographicPoint goal,
								  HashMap<GeographicPoint,GeographicPoint> parentMap, Consumer<GeographicPoint> nodeSearched)
	{
		boolean found = false;
		QueueElement qElem;
		GeographicPoint curr;
		
		HashMap<GeographicPoint,Double> distToNode = new HashMap<GeographicPoint,Double>();
		for (GeographicPoint pt: getVertices()) {
			distToNode.put(pt, Double.POSITIVE_INFINITY);
		}
		PriorityQueue<QueueElement> pQueue = new PriorityQueue<>();	
		HashSet<GeographicPoint> visited = new HashSet<GeographicPoint>();

		distToNode.put(start,0.0);
		pQueue.offer(new QueueElement(start,distToNode.get(start)));

		
	    while (!pQueue.isEmpty()) {
	    	qElem=pQueue.poll();
	    	curr = qElem.getPoint();
	        if (!visited.contains(curr)) {
//	        	PrintQueueElement("DIJKSTRA", qElem);
	        	visited.add(curr);
//	    		numVisited++;
	        	nodeSearched.accept(curr);
	        	if (curr.equals(goal)) {
	        		found = true;
	        		break;
	        	}
	        	dijkstraEnqueueNeighbors(qElem,goal,distToNode,visited,parentMap,pQueue);        	
	        }
	    }
//	    System.out.println("Nodes Visited: " + visited.size()+"\n");
		
		return found;
		
	}
	
	/** Helper method to evaluate each of the neighbors of the current queue element, and enqueue them if necessary
	 * 
	 * @param qElem - The queue element (point and distance info) most recently removed from the Priority Queue
	 * @param goal - destination point
	 * @param distToNode - Map of the current shortest known distance from the start of the path to the given point
	 * @param visited - list of points that have already been visited
	 * @param parentMap - Map that reflects the connectivity of the path from start to goal (if one exists)
	 * @param pQueue - The PriorityQueue.
	 */
	private void dijkstraEnqueueNeighbors(QueueElement qElem, GeographicPoint goal, HashMap<GeographicPoint, Double> distToNode,
									      HashSet<GeographicPoint> visited, 
									      HashMap<GeographicPoint,GeographicPoint> parentMap, PriorityQueue<QueueElement> pQueue)
	{
        GeographicPoint curr = qElem.getPoint();
		double distFromStart = qElem.getDistFromStart();
		double distToNext;
		QueueElement nextElem;

		Set<GeographicPoint> neighbors = vertices.get(curr).getNeighborPoints();
		if (neighbors != null ) {
			for (GeographicPoint next : neighbors) {
				if (!visited.contains(next)) {
					MapEdge neighbor = vertices.get(curr).getEdge(next);
					distToNext = distFromStart + neighbor.getRoadLength();
					if (distToNext < distToNode.get(next)) {
						parentMap.put(next, curr);
						nextElem = new QueueElement(next,distToNext);
						pQueue.offer(nextElem);
						distToNode.put(next, nextElem.priorityDistance());
					}
				}
			}
		}
	}
	
	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return aStarSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, 
											 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		if ((start==null) || (goal==null)) {
			System.out.println("Either start or goal points are null -- no route exists.");
			return null;   
		} 
		HashMap<GeographicPoint,GeographicPoint> parentMap = new HashMap<GeographicPoint,GeographicPoint>();

		boolean found = false;
		QueueElement qElem;
		GeographicPoint curr;
		
		HashMap<GeographicPoint,Double> distPriority = new HashMap<GeographicPoint,Double>();
		for (GeographicPoint pt: getVertices()) {
			distPriority.put(pt, Double.POSITIVE_INFINITY);
		}
		PriorityQueue<QueueElement> pQueue = new PriorityQueue<>();	
		HashSet<GeographicPoint> visited = new HashSet<GeographicPoint>();

		
		pQueue.offer(new QueueElement(start,0.0,goal));


	    while (!pQueue.isEmpty()) {
	    	qElem=pQueue.poll();
        	PrintQueueElement("A*      ", qElem);
	    	curr = qElem.getPoint();
	        if (!visited.contains(curr)) {
//	        	System.out.println("Visited: " + curr);
	        	visited.add(curr);
//	    		numVisited++;
	        	nodeSearched.accept(curr);
	        	if (curr.equals(goal)) {
	        		found = true;
	        		break;
	        	}
	        	aStarEnqueueNeighbors(qElem,goal,distPriority,visited,parentMap,pQueue);        	        	
	        }
	    }

//	    System.out.println("Nodes Visited: " + visited.size()+"\n");

		if (!found) {
			System.out.println("No path exists");
			return null;
		}

		return reconstructPath(start, goal, parentMap);	
	}

	
	/** Helper code to process the neighbors of the current element - added for readability
	 * 
	 * @param qElem - The queue element (point and distance info) most recently removed from the Priority Queue
	 * @param goal - desitnation point
	 * @param distPriority - Map of the current shortest predicted distance from the start of the path to the goal through the given point
	 * @param visited - list of points that have already been visited
	 * @param parentMap - Map that reflects the connectivity of the path from start to goal (if one exists)
	 * @param pQueue - The PriorityQueue.
	 */
	private void aStarEnqueueNeighbors(QueueElement qElem, GeographicPoint goal, HashMap<GeographicPoint, Double> distPriority,
									      HashSet<GeographicPoint> visited, 
									      HashMap<GeographicPoint,GeographicPoint> parentMap, PriorityQueue<QueueElement> pQueue)
	{
        GeographicPoint curr = qElem.getPoint();
		double distFromStart = qElem.getDistFromStart();
		double distToNext;
		QueueElement nextElem;

		Set<GeographicPoint> neighbors = vertices.get(curr).getNeighborPoints();
		if (neighbors != null ) {
			for (GeographicPoint next: neighbors) {
				if (!visited.contains(next)) {
					MapEdge neighbor = vertices.get(curr).getEdge(next);
					distToNext = distFromStart + neighbor.getRoadLength();
					nextElem = new QueueElement(next,distToNext,goal);
					if (nextElem.priorityDistance() < distPriority.get(next)) {
						parentMap.put(next, curr);
						pQueue.offer(nextElem);
						distPriority.put(next,nextElem.priorityDistance());
					}
				}
			}
		}
		
	}
	
	/** Given an in-order list of points representing each intersection in a path, find the length of the path
	 *  by traversing the each edge of the graph.
	 *  
	 * @param path - The in-order list of points
	 * @return - the total distance as a double
	 */
	private Double calcPathDistance(List<GeographicPoint> path) 
	{
		double distance = 0.0;
		GeographicPoint curr;
		GeographicPoint next;
		
		//System.out.println("CalcDistance: Path = "+path);
		
		ListIterator<GeographicPoint> it = path.listIterator();
        curr = it.next(); 
        while (it.hasNext()) {
        	next = it.next();
        	distance += vertices.get(curr).getEdge(next).getRoadLength();
        	curr = next;
        }
	   
		return distance;
	}
	
	/** Given the a list of points traversed in order, return the total distance traveled by looking up
	 *  distance between points (calculated by explicit traversal previously) from a distanceMap
	 *  
	 * @param path  The in-order list of points 
	 * @param distanceMap - data structure that contains the distance for shortest routes between two points.
	 * @return - the total distance
	 */
	private Double calcTotalPathDistance(List<GeographicPoint> path, HashMap<GeographicPoint,TSP_Node> distanceMap)
	{
		double distance = 0.0;
	
		GeographicPoint curr;
		GeographicPoint next;
		
		//System.out.println("CalcDistance: Path = "+path);
		
		ListIterator<GeographicPoint> it = path.listIterator();
        curr = it.next(); 
        while (it.hasNext()) {
        	next = it.next();
        	distance += distanceMap.get(curr).getDistToPoint(next);
        	curr = next;
        }
	   
		return distance;		
	}
	
	
	/** Given a list of points to visit, use aStar search to find the shortest path between the two points for all possible 
	 *  pairs of points from a list of points to visit. Store the distance of each path in an HashMap 
	 *  Throws an IllegalArgumentException if there is a pair of points that are unreachable (no path exists between them).
	 *  
	 * @param verticesTSP 
	 * @return - HashMap of TSP_Nodes - that contain the lengths of all paths from the key node to all other nodes
	 */
	private HashMap<GeographicPoint,TSP_Node> TSPDataInit(List<GeographicPoint> verticesTSP) 
											  throws IllegalArgumentException
	{
		HashMap<GeographicPoint,TSP_Node> mapInit = new HashMap<GeographicPoint,TSP_Node>();
		TSP_Node node_TSP;
		List<GeographicPoint> path;
		
		for (GeographicPoint from: verticesTSP) {
			node_TSP = new TSP_Node(from);
			for (GeographicPoint to : verticesTSP) {
			   if (!from.equals(to)) {
				   path = aStarSearch(from,to);
				   if (path == null) {
					   throw new IllegalArgumentException("No path from "+from+" to "+to+" found. TSP Search failed.");
				   }
				   double pathDistance = calcPathDistance(path);
				   node_TSP.addPath(to,path);
				   node_TSP.setDistanceToPoint(to,pathDistance);
				   mapInit.put(from,node_TSP);
				 //  System.out.println("Path from "+from+" to "+to+" ==>" + path);
				 //  System.out.println("Distance = " + pathDistance);
				   
			   }
			}
		}
		return mapInit;
	}
	
	/** Helper method to determine whether we have visited all the nodes. Returns true if all nodes have been visited 
	 * 
	 * @param mapNodeDistance - used compare all nodes vs the list of nodes visited
	 * @param visited  - the list of nodes that have been visited
	 * @return  - true if all nodes have been visited; false otherwise
	 */
	private boolean visitedAllTSPNodes(HashMap<GeographicPoint,TSP_Node> mapNodeDistance,Set<GeographicPoint> visited)
	{
		for (GeographicPoint pt : mapNodeDistance.keySet()) {
			if (!visited.contains(pt)) return false;
		}
		return (mapNodeDistance.keySet().size()==visited.size());
	}
	
	/** Executes the greedy search algorithm for the TSP problem: From any given node, find an adjacent node that
	 *  a) has not been visited and
	 *  b) with the shortest distance.
	 *  A node that has been visited already cannot be used again. 
	 *  
	 * @param start - the starting (and ending) node of the search
	 * @param mapNodeDistance - data structure containing the distance from this point to all other points to possibly
	 *                          be visited
	 * @param visited - the list of nodes that have been visited - so that no repeats happen
	 * @return - an ArrayList of the identified path in order.
	 */
	private ArrayList<GeographicPoint> greedySearch(GeographicPoint start, 
			                                         HashMap<GeographicPoint,TSP_Node> mapNodeDistance, 
											         Set<GeographicPoint> visited) 
	{
        TSP_Node curr_node;
		ArrayList<GeographicPoint> greedyPath = new ArrayList<GeographicPoint>();
//		debugPathList = new ArrayList<ArrayList <GeographicPoint>>();
		GeographicPoint curr_pt = start;
		visited.add(curr_pt);
//		numVisited++;
		greedyPath.add(curr_pt);
		while (!visitedAllTSPNodes(mapNodeDistance,visited)) {
			curr_node = mapNodeDistance.get(curr_pt);
			GeographicPoint next_pt = curr_node.findMinPath(visited);
		    if (next_pt == null) {
		    	System.out.println("Could not find path from "+curr_pt+" to another node");
		    	return null;
		    }
		    visited.add(next_pt);
//			numVisited++;
		    greedyPath.add(next_pt);
		    curr_pt = next_pt;
		}
		// got here - need to add path back to start;
		curr_node = mapNodeDistance.get(curr_pt);
        double distToStart = curr_node.getDistToPoint(start);
		if (distToStart == Double.POSITIVE_INFINITY) {
			System.out.println("Could not from path from "+curr_pt+" back to "+start);
			return null;
		}
		greedyPath.add(start);
//		debugPathList.add(new ArrayList<GeographicPoint>(greedyPath));
		
		return greedyPath;
	}
	
	/** Helper method to implement the 2-opt optimization. It generates a proposed path optimization
	 *  by switching the nodes in locations n1 and n2, and reversing the order of all nodes between them
	 * @param n1 - index of node to be swapped - guaranteed to be < n2
	 * @param n2 - index of node to be swapped - guaranteed to be > n1
	 * @param path - the existing path
	 * @return - proposed new path optimization
	 */
	private ArrayList<GeographicPoint> swap2opt (int n1, int n2, ArrayList<GeographicPoint> path)
	{
		ArrayList<GeographicPoint> propPath = new ArrayList<GeographicPoint>();
		for (int i = 0; i < n1; i++) 
			propPath.add(path.get(i));
        
		for (int i=n2; i>= n1; i--)
			propPath.add(path.get(i));
		
		for (int i = n2+1; i < path.size(); i++) 
			propPath.add(path.get(i));
		
		return propPath;
		
	}
	
	/** implements the 2opt optimization of a path found by the greedy algorithm
	 * @param path - the initial path found by the greedy search
	 * @param pathLength - the overall length of the initial path
	 * @param mapDistance - the mapping of distances between nodes.
	 * @return - the best path found after completing the optimization - no worse than the original path.
	 */
	private ArrayList<GeographicPoint> optimizeGreedyPath(ArrayList<GeographicPoint> path, 
            										 double pathLength, 
            										 HashMap<GeographicPoint,TSP_Node> mapDistance) 
	{	

		// 2-opt can only run on vertices between the first and last nodes...
		int minIndex = 1;
		int maxIndex = path.size()-2; // can't go further than next to last node.
        int i,j;
        double greedyPathDist = pathLength;
        double propPathDist;
        boolean betterPathFound = false;
		ArrayList<GeographicPoint> propPath;
		
		
        i= minIndex;
		while (i < maxIndex) {
			j= i+1;
			while (j <= maxIndex) {
				propPath = swap2opt(i,j,path);
				propPathDist = calcTotalPathDistance(propPath,mapDistance);
				if (propPathDist < pathLength) {
					System.out.println("Found New path: "+propPath + "  Distance: "+propPathDist);
					betterPathFound = true;
					path = propPath;
					pathLength = propPathDist;
//					debugPathList.add(new ArrayList<GeographicPoint>(propPath));
					//reset counters and start again...
					i = minIndex;
					j = i+1;
				} else {
					j++;
				}
			}
			i++;
		}
		
		if (betterPathFound) {
			System.out.println("Greedy 2-opt found a better path: "+path);
			System.out.println("Greedy Path Length = " + greedyPathDist);
			System.out.println("2-opt  Path Length = " + pathLength);
			System.out.println("Distance Delta     = " + (greedyPathDist-pathLength));
		}
		
		return path;
	
	}

	/** Helper function for debugging to show the path connectivity
	 * @param path - the in-order points in the path 
	 */
	private void printPath(List<GeographicPoint> path)
	{
		boolean first = true;
		for (GeographicPoint pt: path) {
			if (!first) System.out.print(" ==> ");
			System.out.print(pt);
			first = false;
		}
		System.out.println();
	}
	
	/** implements the TSP search - takes a list of points and then tries to find a path that connects them
	 *  such that each point is visited only once before returning to the start point. The start point is assumed to be
	 *  the first point of the list;
	 *  
	 * @param verticesTSP - list of points to find a connecting path through
	 * @return - an optimized path connecting the points such that each point is visited once.
	 */
	public ArrayList<GeographicPoint> TSPSearch(List<GeographicPoint> verticesTSP) 
									  throws IllegalArgumentException
	{
		double pathLength;
		
		for (int i = 0; i<verticesTSP.size();i++) {
			if (!vertices.containsKey(verticesTSP.get(i))) {
				throw new IllegalArgumentException("Point at " + verticesTSP.get(i) + "not found in graph.");
			}
		}
		
		
		// 0. Initialization data structure for distances between the points
		Set<GeographicPoint> explicitlyVisited = new HashSet<GeographicPoint>();
		HashMap<GeographicPoint,TSP_Node> mapTwoVerticesDist = TSPDataInit(verticesTSP);

		ArrayList<GeographicPoint> solution = greedySearch(verticesTSP.get(0),mapTwoVerticesDist,explicitlyVisited);
		printPath(solution);
		pathLength = calcTotalPathDistance(solution,mapTwoVerticesDist);
		System.out.println("Total Greedy Path Length = " + pathLength);
		if (solution != null) 
			solution = optimizeGreedyPath(solution,pathLength,mapTwoVerticesDist);
		// 3. execute greedy algorithm to find a path. Ignore inflight vertices
		//4 return path.
		return solution;
	}

	/** Helper Method to all the user to specify the start point of the TSP search. Reorders the list of nodes to 
	 *  search such that the desired node is at the top of the list, and then calls the default TSPSearch method
	 *  
	 * @param verticesTSP - list of points to search for a path
	 * @param start - the point to start at
	 * @return - the optimized path (if found).
	 */
	public ArrayList<GeographicPoint> TSPSearch(List<GeographicPoint> verticesTSP, GeographicPoint start) 
	{

		List<GeographicPoint> orderedList = new ArrayList<GeographicPoint>();

		orderedList.add(start);
		for (GeographicPoint pt: verticesTSP ) {
			if (!pt.equals(start)) {
				orderedList.add(pt);
			}
		}
		
		ArrayList<GeographicPoint> solution = TSPSearch(orderedList);
		return solution;
	}

	/** takes a map with String Keys and GeographicPoint values and creates the inverse map
	 *  This allows mapping a point to named location (such as a City in the capital searches)
	 * @param aMap
	 * @return
	 */
	HashMap<GeographicPoint,String> reverseMap(HashMap<String,GeographicPoint> aMap)
	{
		HashMap<GeographicPoint,String> revMap = new HashMap<GeographicPoint,String>();
		
		for (String key: aMap.keySet()) {
			revMap.put(aMap.get(key), key);
		}
		return revMap;
	}
	
	
	/** Method to visualize the graph and connections in a text format. Used primarily for debug
	 * 
	 */
	public void printGraph()
	{
		System.out.println("Printing Graph information:\n");
		
		for (GeographicPoint pt : getVertices()) {
			System.out.println("Vertice at ("+pt+"):");
			for (GeographicPoint next : vertices.get(pt).getNeighborPoints()) {
				MapEdge segment = vertices.get(pt).getEdge(next);
				System.out.println("   Connected by "+ segment.getRoadName()+" to ("+next+") Distance =" + segment.getRoadLength() );
			}
		}
	}
	
	public static void main(String[] args)
	{
		System.out.print("Making a new map...");
		MapGraph firstMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...\n");
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", firstMap);
		System.out.println("Printing all vertices");
		System.out.println(firstMap.getVertices());
		firstMap.printGraph();
		System.out.println("DONE.");
		GeographicPoint testStart = new GeographicPoint(1.0, 1.0);
		GeographicPoint testEnd = new GeographicPoint(8.0, -1.0);
		List<GeographicPoint> path = firstMap.bfs(testStart,testEnd);
		System.out.println("Path: "+path);
		
		System.out.print("Making a new map...");
		MapGraph secondMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...\n");
		GraphLoader.loadRoadMap("data/graders/mod2/map2.txt", secondMap);
		secondMap.printGraph();
		System.out.println("DONE.");
		testStart = new GeographicPoint(6, 6);
		testEnd = new GeographicPoint(0, 0);
		path = secondMap.bfs(testStart,testEnd);
		System.out.println("Path: " + path);
		
	
		
		// You can use this method for testing.  
		/* Here are some test cases you should try before you attempt 
		 * the Week 3 End of Week Quiz, EVEN IF you score 100% on the 
		 * programming assignment.
		 */
		
		MapGraph simpleTestMap = new MapGraph();
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);
		
        testStart = new GeographicPoint(1.0, 1.0);
		testEnd = new GeographicPoint(8.0, -1.0);
		
		System.out.println("Test 1 using simpletest: Dijkstra should be 9 and AStar should be 5");
		List<GeographicPoint> testroute = simpleTestMap.dijkstra(testStart,testEnd);
		System.out.println(testroute);
		
		List<GeographicPoint> testroute2 = simpleTestMap.aStarSearch(testStart,testEnd);
		System.out.println(testroute2);		
		MapGraph testMap = new MapGraph();
		GraphLoader.loadRoadMap("data/maps/utc.map", testMap);
		
		// A very simple test using real data
		testStart = new GeographicPoint(32.869423, -117.220917);
		testEnd = new GeographicPoint(32.869255, -117.216927);
		System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5");
		testroute = testMap.dijkstra(testStart,testEnd);
		System.out.println(testroute);
		
		testroute2 = testMap.aStarSearch(testStart,testEnd);
		System.out.println(testroute2);
		
		
		// A slightly more complex test using real data
		testStart = new GeographicPoint(32.8674388, -117.2190213);
		testEnd = new GeographicPoint(32.8697828, -117.2244506);
		System.out.println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10");
		testroute = testMap.dijkstra(testStart,testEnd);
		System.out.println(testroute);
		testroute2 = testMap.aStarSearch(testStart,testEnd);
		System.out.println(testroute2);		
		
		
		/* Use this code in Week 3 End of Week Quiz */
		MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		System.out.println("DONE.");

		GeographicPoint start = new GeographicPoint(32.8648772, -117.2254046);
		GeographicPoint end = new GeographicPoint(32.8660691, -117.217393);
		
		
		List<GeographicPoint> route = theMap.dijkstra(start,end);
		System.out.println("Dijkstra: "+route);
		List<GeographicPoint> route2 = theMap.aStarSearch(start,end);
		System.out.println("A*      : "+route2);

//		System.out.println("\n\nTesting Extension now - simpletest.map:");
//		//firstMap.printGraph();
//		route.clear();
//		route.add(new GeographicPoint(4,1));
//		route.add(new GeographicPoint(5,1));
//		route.add(new GeographicPoint(7,3));
//		route.add(new GeographicPoint(6.5,0));
//		route.add(new GeographicPoint(8,-1));
//        firstMap.TSPSearch(route);		
//		
//        //second test
//        //theMap.printGraph();
//		System.out.println("\n\nTesting utc.map - 8 random points (must be routeable:");
//        route.clear();
//        route.add(new GeographicPoint(32.8649272,-117.229097));
//        route.add(new GeographicPoint(32.8676506,-117.2271116));        
//        route.add(new GeographicPoint(32.8701728,-117.2258286));        
//        route.add(new GeographicPoint(32.868629,-117.215393));        
//        route.add(new GeographicPoint(32.8665949,-117.2164918));        
//        route.add(new GeographicPoint(32.8660691,-117.217393));        
//        route.add(new GeographicPoint(32.864117,-117.216476));        
//        route.add(new GeographicPoint(32.86348,-117.221451));        
//        theMap.TSPSearch(route);
//		
//        System.out.println("\n\nTesting State Capitals");
//        MapGraph capitalsGraph = new MapGraph();
//        HashMap<String,GeographicPoint> capitalsMap = new HashMap<String,GeographicPoint>();
//        GraphLoader.loadCapitals("data/testdata/statecapitals.map",capitalsGraph, capitalsMap);
//        HashMap<GeographicPoint,String> reverseCapMap = capitalsGraph.reverseMap(capitalsMap);
//        //capitalsGraph.printGraph();
//        
//        route.clear();
//        for (GeographicPoint cap : capitalsGraph.getVertices()) {
//        	route.add(cap);
//        }
//        ArrayList<GeographicPoint> capPath = capitalsGraph.TSPSearch(route,capitalsMap.get("Sacramento"));
//       System.out.println("Greedy 2-opt path through all state capitals, starting in Sacramento");
//       for (int i = 0; i < capPath.size(); i++) {
//    	   System.out.println("   "+reverseCapMap.get(capPath.get(i)));
//       }
//       
//       route.clear();
//       System.out.println("\n\nTesting Greedy 2-opt with subset of capitals");
//       String[] capList = {"Phoenix", "Albany","Columbus","Cheyenne","Baton Rouge","St. Paul", "Jackson","Helena","Springfield","Honolulu","Topeka"};
//       for (int i = 0; i < capList.length; i++) {
//    	   route.add(capitalsMap.get(capList[i]));
//       }
//       capPath = capitalsGraph.TSPSearch(route);
//       for (int i = 0; i < capPath.size(); i++) {
//    	   System.out.println("   "+reverseCapMap.get(capPath.get(i)));
//       }
//      
//       System.out.println("\n\nTesting Greedy 2-opt with subset of capitals - start in St. Paul");
//       capPath = capitalsGraph.TSPSearch(route,capitalsMap.get("St. Paul"));
//       for (int i = 0; i < capPath.size(); i++) {
//    	   System.out.println("   "+reverseCapMap.get(capPath.get(i)));
//       }    
//       
//       System.out.println("\n\nTesting Intel Cities North America");
//       MapGraph intelGraph = new MapGraph();
//       HashMap<String,GeographicPoint> intelMap = new HashMap<String,GeographicPoint>();
//       GraphLoader.loadCapitals("data/testdata/intelcities.map",intelGraph, intelMap);
//       HashMap<GeographicPoint,String> reverseIntelMap = intelGraph.reverseMap(intelMap);
//       route.clear();
//       for (GeographicPoint city : intelGraph.getVertices()) {
//       	route.add(city);
//    	   System.out.println("   "+reverseIntelMap.get(city)+"  "+city);
//       }
//      ArrayList<GeographicPoint> intelPath = intelGraph.TSPSearch(route, intelMap.get("Santa Clara"));
//      for (int i = 0; i < intelPath.size(); i++) {
//   	   System.out.println("   "+reverseIntelMap.get(intelPath.get(i)));
//      }    
//      
//      String city;
//      System.out.println("printing greedy paths + successful optimizations");
//      for (int i = 0; i < intelGraph.debugPathList.size(); i++) {
//    	  for (int j = 0; j < intelGraph.debugPathList.get(i).size(); j++) {
//    		  city = reverseIntelMap.get(intelGraph.debugPathList.get(i).get(j));
//    		  if (j > 0) city = " -> "+city;
//    		  System.out.print(city);
//    	  }
//    	  System.out.println();
//      }
      
	}
}
