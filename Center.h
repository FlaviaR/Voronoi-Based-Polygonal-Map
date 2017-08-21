#include <list>
#include <CGAL/Cartesian.h>

#include "Corner.h"

typedef CGAL::Cartesian<double> K;  ///< Kernel
typedef CGAL::Point_2<K>  Point;  ///< Point in 2D
typedef std::list<Point>  PointList;
typedef std::list<Corner> CornerList;

// The Center object is responsible to keep track of the relationships between Voronoi faces
class Center {
	public:
		// this center's index within centerList
		int index;
	    // site of voronoi face
		Point point;
		// elevation of this center object
		float elevation;
	
		// is this tile water?
		bool isWater; 
		// water can be either an ocean or a lake
		// if ocean is false but isWater is set, then its assumed to be a lake
		bool ocean;    
		// depictor for the border of the map
		bool isBorder;

	
		// List of the corners of this voronoi face
		CornerList corners;
		// List of the neighbors of this voronoi face
		std::list<Center> neighbors;
		// List of the edges of this voronoi face
		//EdgeList edges;
	
	
	
	// (A friend function of a class is defined outside that class' scope but it has the right to access all private and protected members of the class. Even though the prototypes for friend functions appear in the class definition, friends are not member functions.)
	// Print the data from center
	friend std::ostream& operator <<(std::ostream& out, const Center& center ) {
		out << "\n Site " << center.index;
		out << "\n location: ("	<< center.point[0]  << " , " << center.point[1] << ")";
		out << "\n isWater: "	<< center.isWater;
		out << "\n ocean: "		<< center.ocean;
		out << "\n isBorder: "	<< center.isBorder;
		
		CornerList corners = center.corners;
		CornerList::iterator corner_iter = corners.begin();
		
		while (corner_iter != corners.end()) {
			Corner corner = *corner_iter++;
			std::cout << corner << std::endl ;

		}
		
		std::list<Center> neighbors = center.neighbors;
		std::list<Center>::iterator neigh_iter = neighbors.begin();
		
		while (neigh_iter != neighbors.end()) {
			Center n = *neigh_iter++;
			std::cout << "\n Neighbor: (" << n.point[0] << ", " << n.point[1] << ")" << std::endl ;
			
		}
		
		return out ;
	}
};
