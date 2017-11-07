package structures;

import java.util.ArrayList;

/**
 * Encapsulates an interval tree.
 * 
 * @author runb-cs112
 */
public class IntervalTree {
	
	/**
	 * The root of the interval tree
	 */
	IntervalTreeNode root;
	
	/**
	 * Constructs entire interval tree from set of input intervals. Constructing the tree
	 * means building the interval tree structure and mapping the intervals to the nodes.
	 * 
	 * @param intervals Array list of intervals for which the tree is constructed
	 */
	public IntervalTree(ArrayList<Interval> intervals) {
		
		// make a copy of intervals to use for right sorting
		ArrayList<Interval> intervalsRight = new ArrayList<Interval>(intervals.size());
		for (Interval iv : intervals) {
			intervalsRight.add(iv);
		}
		
		// rename input intervals for left sorting
		ArrayList<Interval> intervalsLeft = intervals;
		
		// sort intervals on left and right end points
		sortIntervals(intervalsLeft, 'l');
		sortIntervals(intervalsRight,'r');
		
		// get sorted list of end points without duplicates
		ArrayList<Integer> sortedEndPoints = 
							getSortedEndPoints(intervalsLeft, intervalsRight);
		
		// build the tree nodes
		root = buildTreeNodes(sortedEndPoints);
		
		// map intervals to the tree nodes
		mapIntervalsToTree(intervalsLeft, intervalsRight);
	}
	
	/**
	 * Returns the root of this interval tree.
	 * 
	 * @return Root of interval tree.
	 */
	public IntervalTreeNode getRoot() {
		return root;
	}
	
	/**
	 * Sorts a set of intervals in place, according to left or right endpoints.  
	 * At the end of the method, the parameter array list is a sorted list. 
	 * 
	 * @param intervals Array list of intervals to be sorted.
	 * @param lr If 'l', then sort is on left endpoints; if 'r', sort is on right endpoints
	 */
	public static void sortIntervals(ArrayList<Interval> intervals, char lr) {
		System.out.println("Sorting...");
		Queue<Interval> sort = new Queue<Interval>();
		ArrayList<Interval> temp = new ArrayList<Interval>();
		int num = 0;
		if(lr == 'l')
		{	
			while(intervals.size()!=temp.size())
			{
				for(int i=0; i<intervals.size(); i++)
				{
					if(intervals.get(i).leftEndPoint == num)
					{
						temp.add(intervals.get(i));
					}
				}
				num++;
			}
		while(!intervals.isEmpty()) { intervals.remove(0); }
		while(!temp.isEmpty()) { sort.enqueue(temp.get(0)); temp.remove(0); }
		while(!sort.isEmpty()) { intervals.add(sort.dequeue()); }
		}
		else if(lr == 'r')
		{
			while(intervals.size()!=temp.size())
			{
				for(int i=0; i<intervals.size(); i++)
				{
					if(intervals.get(i).rightEndPoint == num)
					{
						temp.add(intervals.get(i));
					}
				}
				num++;
			}
			while(!intervals.isEmpty()) { intervals.remove(0); }
			while(!temp.isEmpty()) { sort.enqueue(temp.get(0)); temp.remove(0); }
			while(!sort.isEmpty()) { intervals.add(sort.dequeue()); }
		}
		for(int i=0; i<intervals.size(); i++)
		{
			System.out.println(intervals.get(i));
		}
	}
	
	/**
	 * Given a set of intervals (left sorted and right sorted), extracts the left and right end points,
	 * and returns a sorted list of the combined end points without duplicates.
	 * 
	 * @param leftSortedIntervals Array list of intervals sorted according to left endpoints
	 * @param rightSortedIntervals Array list of intervals sorted according to right endpoints
	 * @return Sorted array list of all endpoints without duplicates
	 */
	public static ArrayList<Integer> getSortedEndPoints(ArrayList<Interval> leftSortedIntervals, ArrayList<Interval> rightSortedIntervals) {
		ArrayList<Integer> points = new ArrayList<Integer>();
		ArrayList<Integer> tempPoints = new ArrayList<Integer>();
		int positionL = 0, positionR =0;
			while(positionL < leftSortedIntervals.size())
			{
				int leftNum = leftSortedIntervals.get(positionL).leftEndPoint;
				points.add(leftNum);
				for(int i=0; i<leftSortedIntervals.size(); i++)
				{
					if(leftNum == leftSortedIntervals.get(i).leftEndPoint)
					{
						positionL++;
					}
				}
			}
			while(positionR < rightSortedIntervals.size())
			{
				int leftNum = rightSortedIntervals.get(positionR).rightEndPoint;
				tempPoints.add(leftNum);
				for(int i=0; i<rightSortedIntervals.size(); i++)
				{
					if(leftNum == leftSortedIntervals.get(i).rightEndPoint)
					{
						positionR++;
					}
				}
			}
			
			for(int i=0; i<points.size(); i++)
			{
				for(int j=0; j<tempPoints.size(); j++)
				{
					if(points.get(i)==tempPoints.get(j))
					{
						tempPoints.remove(j);
					}
				}
			}
			
			while(tempPoints.size()!=0)
			{
				int count = 0;
				int spot = 0;
				int num = tempPoints.get(0);
				for(spot=0; spot<points.size(); spot++)
				{
					if(points.get(spot)>num) { break; }
					count++;
				}
				points.add(count, num);
				tempPoints.remove(0);
			}
		/*System.out.print("Points: [");
			for(int k=0; k<points.size()-1; k++)
			{
				System.out.print(points.get(k) + ", ");
			}	
			System.out.println(points.get(points.size()-1) + "]");*/
		
		return points;
	}
	
	/**
	 * Builds the interval tree structure given a sorted array list of end points
	 * without duplicates.
	 * 
	 * @param endPoints Sorted array list of end points
	 * @return Root of the tree structure
	 */
	public static IntervalTreeNode buildTreeNodes(ArrayList<Integer> endPoints) {
		Queue<IntervalTreeNode> tree = new Queue<IntervalTreeNode>();
		IntervalTreeNode T;
		for(int i=0; i<endPoints.size(); i++)
		{
			T = new IntervalTreeNode(endPoints.get(i), endPoints.get(i), endPoints.get(i)); 
			T.leftIntervals = new ArrayList<Interval>();
			T.rightIntervals = new ArrayList<Interval>();
			tree.enqueue(T);
		}
		
		IntervalTreeNode newTree = null;
		while(tree.size()>0){
		if(tree.size()==1)
		{
			newTree = tree.dequeue();
			return newTree;
		}
		else 
		{
			int tempSize = tree.size();	
			while(tempSize>1)
			{
				IntervalTreeNode T1 = tree.dequeue();
				IntervalTreeNode T2 = tree.dequeue();
				float v1 = T1.maxSplitValue;
				float v2 = T2.minSplitValue;
				float x = (v1+v2)/2;
				IntervalTreeNode N = new IntervalTreeNode(x, T1.minSplitValue, T2.maxSplitValue);
				N.leftIntervals = new ArrayList<Interval>();
				N.rightIntervals = new ArrayList<Interval>();
				N.leftChild = T1;
				N.rightChild = T2;
				tree.enqueue(N);
				tempSize = tempSize-2;
				}
			if(tempSize==1)
			{
				IntervalTreeNode temp = tree.dequeue();
				tree.enqueue(temp);
			}
		}
		}
		return newTree;
	}
	
	/**
	 * Maps a set of intervals to the nodes of this interval tree. 
	 * 
	 * @param leftSortedIntervals Array list of intervals sorted according to left endpoints
	 * @param rightSortedIntervals Array list of intervals sorted according to right endpoints
	 */
	public void mapIntervalsToTree(ArrayList<Interval> leftSortedIntervals, ArrayList<Interval> rightSortedIntervals) {
		IntervalTreeNode temp = root;

		for(int i=0; i<leftSortedIntervals.size(); i++)
		{
			while(!leftSortedIntervals.get(i).contains(root.splitValue))
			{	
				if(leftSortedIntervals.get(i).rightEndPoint < root.splitValue)
				{
					root = root.leftChild;
				}
				else if(leftSortedIntervals.get(i).leftEndPoint > root.splitValue)
				{
					root = root.rightChild;
				}
			}
			root.leftIntervals.add(leftSortedIntervals.get(i));
			root = temp;	
		}
		
		for(int j=0; j<rightSortedIntervals.size(); j++)
		{
			while(!rightSortedIntervals.get(j).contains(root.splitValue))
			{	
				if(rightSortedIntervals.get(j).rightEndPoint < root.splitValue)
				{
					root = root.leftChild;
				}
				else if(rightSortedIntervals.get(j).leftEndPoint > root.splitValue)
				{
					root = root.rightChild;
				}
			}
			root.rightIntervals.add(rightSortedIntervals.get(j));
			root = temp;	
		}
	}
	
	/**
	 * Gets all intervals in this interval tree that intersect with a given interval.
	 * 
	 * @param q The query interval for which intersections are to be found
	 * @return Array list of all intersecting intervals; size is 0 if there are no intersections
	 */
	public ArrayList<Interval> findIntersectingIntervals(Interval q) {
		ArrayList<Interval> ResultList = new ArrayList<Interval>();
		
		if(q.contains(root.splitValue))
		{
			for(int i=0; i<root.leftIntervals.size(); i++)
			{
				ResultList.add(root.leftIntervals.get(i));
			}
			queryTree(root.leftChild, q, ResultList);
			queryTree(root.rightChild, q, ResultList);
		}
		else if(root.splitValue < q.leftEndPoint)
		{
			root.matchRight(q, ResultList);
			queryTree(root.rightChild, q, ResultList);
		}
		else if(root.splitValue > q.rightEndPoint)
		{
			root.matchLeft(q, ResultList);
			queryTree(root.leftChild, q, ResultList);
		}
		return ResultList;
	}
	
	private void queryTree(IntervalTreeNode root, Interval q, ArrayList<Interval> ResultList)
	{
		if(root!=null){
		for(int i=0; i<root.leftIntervals.size(); i++)
		{
			if(q.intersects(root.leftIntervals.get(i)))
			{
				ResultList.add(root.leftIntervals.get(i));
			}
		}
		
		queryTree(root.leftChild, q, ResultList);
		queryTree(root.rightChild, q, ResultList);
		}
		return;
	}
}

