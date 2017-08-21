#include <list>
#include <CGAL/Cartesian.h>

typedef CGAL::Cartesian<double> K;  ///< Kernel
typedef CGAL::Point_2<K>  Point;  ///< Point in 2D

// The Corner object contains information regarding the corners 
// of a specific voronoi face
class Corner {
	public:
		// index in respect to the corner list			  
		int index;
	    // Corner location
		Point point;  

		// is this tile water?
		bool isWater; 
		// water can be either an ocean or a lake
		// if ocean is false but the is water, then its assumed to be a lake
		bool ocean;    

		// is this corner on the border (edges) of the map?
		bool isBorder;
	
		float elevation; // between 0.0 and 1.0
	
	// Print the data from center
	friend std::ostream& operator <<(std::ostream& out, const Corner& corner ) {
		out << "\n \t Corner " << corner.index;
		out << "\n \t point: ("		<< corner.point[0]  << " , " << corner.point[1] << ")";
		out << "\n \t isWater: "	<< corner.isWater;
		out << "\n \t ocean: "		<< corner.ocean;
		out << "\n \t isBorder: "	<< corner.isBorder;
		out << "\n \t Elevation: "	<< corner.elevation;

		return out ;
	}
};
