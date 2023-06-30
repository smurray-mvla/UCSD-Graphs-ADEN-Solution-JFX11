package roadgraph;

import geography.GeographicPoint;


/** Class to hold an element of the PriorityQueue. Relates a point to distance from start and straightline
 * distance to goal.
 * 
 * @author Scott
 *
 */
/**
 * @author Scott
 *
 */
/**
 * @author Scott
 *
 */
/**
 * @author Scott
 *
 */
public class QueueElement implements Comparable<QueueElement> {
	
	private GeographicPoint point;
    private double distFromStart;
    private double distToGoal;
    
    
    /** Constructor for Dijsktra search
     * @param pt - current point
     * @param d1 - distance from start to this pt on current path
     */
    public QueueElement(GeographicPoint pt, double d1) 
    {
    	this.point = pt;
    	this.distFromStart = d1;
    	this.distToGoal = 0.0;
    }
    
    /** Constructor for A* search
     * @param pt - current point
     * @param d1 - distance from start to this point on current path
     * @param goal - straightline distance from this point to goal point
     */
    public QueueElement(GeographicPoint pt, double d1, GeographicPoint goal) 
    {
    	this.point = pt;
    	this.distFromStart = d1;
    	this.distToGoal = pt.distance(goal);
    }

    /** gets the current point
     * @return
     */
    public GeographicPoint getPoint()
    {
    	return point;
    }
    
    /** Gets the distance from start
     * @return
     */
    public double getDistFromStart() 
    {
    	return distFromStart;
    }
    
    /** gets the distance to the goal
     * @return
     */
    public double getDistToGoal()
    {
    	return distToGoal;
    }
    

    /** gets the distance used in the priority queue (sums distance from start with distance to goal).
     * If Dijkstra - then distance to goal is 0;
     * @return
     */
    public double priorityDistance()
    {
    	return distFromStart + distToGoal;
    }

    /** Comparator function for the priority queue. Note same function for both algorithms
     *
     */
    public int compareTo(QueueElement o )
    {
    	return Double.compare(this.distFromStart+this.distToGoal,o.getDistFromStart()+o.getDistToGoal());
    }

}
