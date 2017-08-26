/*
 Generates randomized maps based off of voronoi diagrams
 In mathematics, a Voronoi diagram is a partitioning of a plane into regions based on distance
 to points in a specific subset of the plane. The set of points is specified beforehand, and
 for each seed there is a corresponding region consisting of all points closer to that seed than
 to any other -> Voronoi cells.
 @see https://en.wikipedia.org/wiki/Voronoi_diagram
 In mathematics, a Delaunay triangulation for a set of points P in a plane is a triangulation DT (P)
 such that no point in P is inside the circumcircle of any triangle in DT(P).
 @see https://en.wikipedia.org/wiki/Delaunay_triangulation
 Delaunay triangulation is the dual of Voronoi diagrams, meaning that for every triangulation there
 exists one and only one corresponding Voronoi tesselation. And vice versa.
 The exception being the degenerative cases, eg. when the Voronoi diagrams result in a square grid.
 @author Flavia Cavalcanti
 @since 17/07/2017
 @see http://www-cs-students.stanford.edu/~amitp/game-programming/polygon-map-generation/
 */

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GL/glut.h>
#include <CGAL/Cartesian.h>
#include <list>
#include <algorithm>
#include <time.h>       /* time */
#include <string>
#include <queue>        // std::queue
#include <tuple>
#include <map>			// std::map
#include <cmath>
#include <random>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */


// includes for the definition of the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

#include "Center.h"

#define NULL __null
#define PI 3.14159265

// typedefs for defining the voronoi adaptor
typedef CGAL::Cartesian<double> 											 K;  ///< Kernel
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;


typedef CGAL::Point_2<K> 													 Point;  ///< Point in 2D
typedef CGAL::Segment_2<K> 													 Segment;  ///< Line Segment in 2D
typedef CGAL::Ray_2<K> 													 	 Ray;  ///< Ray (halfline) in 2D
typedef CGAL::Vector_2<K>  													 Vector;  ///< Vector in 2D
typedef CGAL::Line_2<K>  													 Line;  ///< Line in 2D
typedef CGAL::Delaunay_triangulation_2<K>  									 Delaunay; ///< A Delaunay Triangulation
typedef CGAL::Voronoi_diagram_2<DT, AT, AP>  								 Voronoi;
typedef std::list<Point>  													 PointList;
typedef std::list<Center*>  												 CenterList;
typedef std::list<Corner*>  												 CornerList;
typedef VD::Halfedge_handle  												 Halfedge_handle;
typedef Delaunay::Vertex_handle												 Vertex_handle;
typedef Delaunay::Face_circulator											 Face_circulator;
typedef Delaunay::Vertex_circulator											 Vertex_circulator;
typedef boost::tuple<std::map <Point, CornerList>, std::map <Point, Corner*> > TupleMaps;

// See: http://www-cs-students.stanford.edu/~amitp/game-programming/polygon-map-generation/

// Storing two different versions of each generation since any direct modification would
// require recalculations from scratch. These include randomized lists and lists that use the relaxed version of the points.
PointList points;
PointList pointsRand;

Delaunay triang;
Delaunay triangRand;

Voronoi voronoi;
Voronoi voronoiRand;

CenterList centerList;
CornerList cornerList;


// Be sure to use a square map - to make my life easier ( .^. )
int SIZE = 500;
// Width of canvas
int width = SIZE;
// Height of canvas
int height = SIZE;
// Number of generated randomized points
int numOfPoints = 300;
// Number of times to run lloyd's relaxation
int numOfLloydIterations = 2;
// Whether to add new points when clicking the left mouse button - probably will not use this
bool addPointsWithClick = false;

// Whether to display the Voronoi polygons generated via randomely generated points
bool displayRandom = true;
// Whether to display the generated map biomes
bool displayBiome = false;
// Whether to display a radially generated map
bool displayRadial = false;
// Whether to display the Delaunay triangulation of the currently displayed Voronoi polygons
bool displayDelaunay = false;
// Whether to display the Voronoi polygons
bool displayVoronoi = true;

// Whether to run testing methods
bool debug = true;

// --------------------------------------------------------------------------------------


// It is expected that 'points' is filled before calling this function
/// Draw Delaunay Triangulation
/// @param: a delauney triangle object
void drawDelaunayTriang(Delaunay triangle) {
	
	glColor3f (0.0, 0.0, 0.0);
	glBegin (GL_LINES);
	
	Delaunay::Edge_iterator ie = triangle.edges_begin ();
	while (ie != triangle.edges_end ()) {
		Delaunay::Edge e = *ie++;
		if (triangle.is_infinite (e)) continue;
		Delaunay::Segment s = triangle.segment (e);
		for (int i = 0; i < 2; ++i) {
			Point p = s[i];
			glVertex2f (p [0], p [1]);
		}
	}
	glEnd ();
}

// It is expected that var points is filled before calling this function
/// Draw the corresponding voronoi diagram by using its dual delauney triangulation
/// @param: a delauney triangle object
void drawVoronoiGivenDelaunay(Delaunay triangle) {
	if (displayBiome) glColor3f (0.0, 0.0, 0.0);
	else glColor3f (1.0, 1.0, 1.0);
	glBegin (GL_LINES);
	
	Delaunay::Edge_iterator ie = triangle.edges_begin ();
	ie = triangle.edges_begin ();
	while (ie != triangle.edges_end ()) {
		Delaunay::Edge e = *ie++;
		CGAL::Object obj = triangle.dual (e);
		
		Segment s;
		Ray r;
		if (CGAL::assign (s, obj)) {
			glVertex2f (s[0][0], s[0][1]);
			glVertex2f (s[1][0], s[1][1]);
		} else if (CGAL::assign (r, obj)) {
			Point p = r.source();
			Vector v = r.direction().vector();
			glVertex2f (p[0], p[1]);
			p = p + 1000 * v;
			glVertex2f (p[0], p[1]);
		}
	}
	
	glEnd();
	
}

// It is expected that points is filled before calling this function
/// Draws the voronoi diagram by using Cgal's voronoi diagram adaptor
void drawVoronoi(Voronoi voro) {
	voro.insert(points.begin(), points.end());
	Voronoi::Edge_iterator ve = voro.edges_begin ();
	
	/// Draw Voronoi Diagram
	glColor3f (1.0, 0.7, 0.3); // weird orange
	glBegin (GL_LINES);
	
	while (ve != voronoi.edges_end ()) {
		Voronoi::Halfedge v = *ve++;
		Point ps;
		Point pt;
		if (v.has_source()) {
			// v.source() returns a vertex handler
			ps = v.source() -> point();
		}
		if (v.has_target()) {
			pt = v.target() -> point();
		}
		glVertex2f (ps[0], ps[1]);
		glVertex2f (pt[0], pt[1]);
	}
	
	glEnd ();
	
}

// --------------------------------------------------------------------------------------

// Return the distance from (0, 0) to this point
double pointLength (Point p) {
	return sqrt(pow(p[0], 2.0) + pow(p[1], 2.0));
}

// Gernerates a circular map using overlapping sine waves (ensue black magic)
double ISLAND_FACTOR = 1.15; // Factor should ideally be between 1.0 and 2.0
// lower bounds lead to more constricted islands
bool generateRadialMap(Point p) {
	srand (p[0] + p[1] + 0.5);
	int bumps = rand()%5 + 1; // random between 1 and 6
	
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dist(0, 2*PI);
	double dipAngle = dist(gen); // random between 0 and 2PI
	double startAngle = dist(gen); // random between 0 and 2PI
	
	double dipWidth = ((double)(rand()%5)/10) + 0.2; // random between 0.2 and 0.7
	double angle = atan(p[1]/p[0]);
	double length = 0.5 * (std::max(std::abs(p[0]), std::abs(p[1]))) + pointLength(p);
	double r1 = 0.5 + 0.40*sin(startAngle + bumps*angle + cos((bumps + 3) * angle));
	double r2 = 0.7 - 0.20*sin(startAngle + bumps*angle - sin((bumps + 2) * angle));
	
	if (std::abs(angle - dipAngle) < dipWidth or std::abs(angle - dipAngle + 2*PI) < dipWidth or std::abs(angle - dipAngle - 2*PI) < dipWidth) {
		r1 = r2 = 0.2;
	}
	
	bool isInside = (length < r1 or (length > r1 * ISLAND_FACTOR && length < r2));
	return isInside;
}

// Generates a randomized map -- very sparse
bool generateRandomMap(Point p) {
	srand (p[0] + p[1] + 0.5);
	int num = rand();
	if (num % 2 == 0) return false;
	else return true;
}

// CHeck to see whether a given point is inside land or in water
bool isInside(Point p) {
	if (displayRadial) {
		int SIZE = width;
		// The póint has to be normalized to be between -1.0 and +1.0
		Point pt = Point(1.5 * (p[0]/SIZE - 0.5), 1.5 * (p[1]/SIZE - 0.5));
		return generateRadialMap(pt);
	}
	else return generateRandomMap(p);
}

void assignFaceElevations(Center* center) {
	CornerListObj corners = center->corners;
	float sum = 0.0;
	
	for ( CornerListObj::iterator corner_iter = corners.begin();
		 corner_iter != corners.end(); ++corner_iter) {
		sum += (*corner_iter).elevation;
	}
	
	center->elevation = sum/corners.size();
}

// Determine elevation for each Voronoi polygon corner
void assignCornerElevations (void) {
	CornerList::iterator corner_iter = cornerList.begin();
	
	while (corner_iter != cornerList.end()) {
		Corner* corner = *corner_iter++;
		
		corner->isWater = !isInside(corner->point);
		
		if (corner->isBorder) {
			corner->elevation = 0.0;
			// push to queue
		} else {
			// elevations are supposed to be between 0.0 and 1.0 -
			// just add a large number here as a flag
			corner->elevation = 10.0;
			
		}
		
	}
	
}

// If a given center contains a corner with the 'border' attribute,
// then return a modified center object set as a border and ocean
// All border faces are expected to be oceans
void setOceanBorders (Center* center) {
	CornerListObj corners = center->corners;
	CornerListObj::iterator ic = corners.begin ();
	
	while (ic != corners.end ()) {
		Corner corner = *ic++;
		
		if (corner.isBorder)  {
			center->isWater = true;
			center->ocean = true;
			center->isBorder =  true;
		}
		
		if(center->ocean) center->isWater = true;
	}
	
}

// Given a Voronoi site, return a list of Corner objects composed
// from the corners of the voronoi face that corresponds to the given site
// This is performed by using the Delaunay triangulation to fetch the faces incident to the given site
CornerList initCorners (Point site) {

	CornerList corners;
	Vertex_handle handle = triang.nearest_vertex(site);
	Face_circulator f_circ = triang.incident_faces(handle),
	done(f_circ);
	
	do
	{
		// Given an incident face to a given vertex, fetch the dual of this face,
		// which corresponds to the center of the circle circumscribed to face f,
		// AKA: the corresponding corner of a voronoi face
		Point location = triang.dual(f_circ);
		Corner* corner = new Corner();
		corner->point = location;
		corner->index = cornerList.size();
		// In openGL, point (0.0, 0.0) is located at the top left corner of the window
		// Cgal continues generating Voronoi points outside of window space
		// Check to see if any of the points lie outside of the window coordinates
		corner->isBorder = (corner->point[0] <= 0 or corner->point[0] >= SIZE or corner->point[1] <= 0 or corner->point[1] >= SIZE);
		
		corners.push_back(corner);
	} while(++f_circ != done);
	
	return corners;
}

// Initialize Center objects based off of the point list
	TupleMaps initCornersAndMap (void) {
//
	std::map <Point, CornerList> cornerMapV; // key: voronoi site, value: corresponding corners
	std::map <Point, Corner*> cornerMap;    // key: corner site, value: corresponding corner pointer
	boost::tuple<std::map <Point, CornerList>, std::map <Point, Corner*> > retValues;
	
	PointList::iterator pl = points.begin();
	while (pl != points.end()) {
		Point p = *pl++;
		CornerList corners = initCorners (p);

		cornerMapV.insert(std::pair<Point, CornerList> (p, corners));
		
		for (CornerList::iterator corner_iter = corners.begin();
			 corner_iter != corners.end(); ++corner_iter) {
			
			Corner* corner = *corner_iter;
			cornerMap.insert(std::pair<Point, Corner*> (corner->point, corner));
			cornerList.push_back(corner);
		}
		
	}
	
	return boost::make_tuple(cornerMapV, cornerMap);
}

std::list<Corner> fetchAdjacentCorners(Point corner, std::map<Point, Corner*> map) {
	;
}

// Given a Voronoi site, return a list of Corner objects composed
CornerListObj fetchCorners (Center* center) {
	CornerListObj cornerList;
	Vertex_handle handle = triang.nearest_vertex(center->point);
	Face_circulator f_circ = triang.incident_faces(handle),
	done(f_circ);
	
	do
	{
		// Given an incident face to a given vertex, fetch the dual of this face,
		// which corresponds to the center of the circle circumscribed to face f,
		// AKA: the corresponding corner of a voronoi face
		Point location = triang.dual(f_circ);
		Corner corner;
		corner.index = cornerList.size();
		corner.point = location;
		// In openGL, point (0.0, 0.0) is located at the top left corner of the window
		// Cgal continues generating Voronoi points outside of window space
		// Check to see if any of the points lie outside of the window coordinates
		corner.isBorder = (corner.point[0] <= 0 or corner.point[0] >= SIZE or corner.point[1] <= 0 or corner.point[1] >= SIZE);
		
		cornerList.push_back(corner);
		
	} while(++f_circ != done);
	
	return cornerList;
}


// Given a Voronoi site, return a Center object with the following initialized:
// location, isWater, ocean, isBorder, and corners
Center* initCenters (Point site, std::map <Point, CornerList> cornerMapV) {
//Center* initCenters (Point site) {

	Center* center = new Center();
	// CornerList cornerList = fetchCorners(site);
	
	//center->index = centerList.size();
	center->point = site;
	center->isWater = !isInside(center->point);
	CornerList corners = cornerMapV[center->point];
//	CornerListObj corners = fetchCorners(center);
	CornerList::iterator corner_iter = corners.begin();

	// center.corners.insert(corner_iter, cornerList.begin(), cornerList.end());
	
	while (corner_iter != corners.end ()) {
		Corner* c = *corner_iter++;
		Corner corner = *c;
		center->corners.push_back(corner);
	}
	
	setOceanBorders(center);
	
	return center;
}

// Given a Voronoi site, return a list of Centers of the corresponding neighboring faces
// This is performed by using the Delaunay triangulation to fetch the vertices incident to the given site
// Note: these centers will not of yet have their neighboring lists filled
//std::list<Center> fetchNeighbors(Point site, std::map<Point, Center*> map) {
std::list<Center> fetchNeighbors(Point site, std::map<Point, Center*> map) {
	std::list<Center> neighboringCenters;
	Vertex_handle handle = triang.nearest_vertex(site);
	Vertex_circulator vc = triang.incident_vertices(handle),
	done(vc);
	
	if (vc != 0) {
		do {
			
			Point v_site = vc -> point();
			
			if (map.find(v_site) != map.end()) {
				Center c = *map[v_site];
				//std::cout << v_site << " -> " << *map.find(v_site)->second << std::endl;
				
				neighboringCenters.push_back (c);
			}
			
		} while(++vc != done);
	}
	
	return neighboringCenters;
}

// Initialize Center objects based off of the point list
std::map <Point, Center*> initCentersAndMap (std::map <Point, CornerList> cornerMapV) {
//std::map <Point, Center*> initCentersAndMap () {
	int index = 0;
	std::map <Point, Center*> centerMap;
	
	PointList::iterator pl = points.begin();
	while (pl != points.end()) {
		Point p = *pl++;
		Center* c = initCenters (p, cornerMapV);
//		Center* c = initCenters (p);
		
		c->index = index;
		
		centerMap.insert(std::pair<Point, Center*> (c->point, c));
		centerList.push_back(c);
		
		index++;
	}
	
	return centerMap;
}

// Generates relationship graph
void buildRelationshipLists () {
	int index = 0;
	
	TupleMaps cornerMaps = initCornersAndMap();
	std::map <Point, CornerList> cornerMapV = boost::get<0>(cornerMaps);
	std::map <Point, Corner*> cornerMap  = boost::get<1>(cornerMaps);
	std::map <Point, Center*> centerMap  = initCentersAndMap(cornerMapV);
//	std::map <Point, Center*> centerMap  = initCentersAndMap();

	// Establish the neighbors of each polygonal face
	for ( CenterList::iterator center_iter = centerList.begin();
		 center_iter != centerList.end(); ++center_iter) {
		
		Center* center = *center_iter;
		center -> neighbors = fetchNeighbors(center->point, centerMap);
		
		// if (debug) std::cout << *center << std::endl;
	}
	
}

// As randomely generated points lead to more deformed polygons,
// semi-random numbers or quasi-randomely generated points will be used.
// This will be approximated by using a variation of Lloyd Relaxation to make
// the points more evenly distributed.
// Lloyd Relaxation replaces each point by the centroid of the polygon encapsulating it.
void performPointDistribution() {
	Voronoi voronoiRedestributed;
	Delaunay triangRedestributed;
	
	for (int i = 0; i < numOfLloydIterations; i++) {
		PointList toRet;
		Voronoi voronoiTemp;
		voronoiTemp.insert(points.begin(), points.end());
		
		Voronoi::Face_iterator vf = voronoiTemp.faces_begin();
		
		while (vf != voronoiTemp.faces_end()) {
			Voronoi::Face f = *vf++;
			
			Voronoi::Ccb_halfedge_circulator ccb_start = f.outer_ccb();
			Voronoi::Ccb_halfedge_circulator halfEdgeCirc = ccb_start;
			
			
			int length = 0;
			double px = 0.0; // the x coordinate of the new point
			double py = 0.0; // the y coordinate of the new point
			
			do {
				if (halfEdgeCirc->has_source()) {
					Point ps = halfEdgeCirc->source() -> point();
					px += ps[0];
					py += ps[1];
					length += 1;
				}
			} while (++halfEdgeCirc != ccb_start);
			
			px /= length;
			py /= length;
			
			Point p (px, py);
			toRet.push_back(p);
		}
		
		points = toRet;
	}
	
	voronoiRedestributed.insert(points.begin(), points.end());
	triangRedestributed.insert(points.begin(), points.end());
	
	voronoi = voronoiRedestributed;
	triang = triangRedestributed;
}



// Set points array to contain a series of randomely generated points
// Notice how this leads to lumpy and irregular polygon sizes if drawn without the relaxation.
void generateRandPoints() {
	PointList pointList;
	PointList randomList;
	Delaunay  triangle;
	
	for (int i = 0; i < numOfPoints ; i ++) {
		
		int x = rand() % width  + 1;
		int y = rand() % height + 1;
		
		Point p ((double) x, (double) y);
		pointList.push_back(p);
		randomList.push_back(p);
		
		triangle.push_back(p);
	}
	points = pointList;
	pointsRand = randomList;
	triangRand = triangle;
}

// Deletes all pointers from centerList
void resetPointerLists(void) {
	for ( CenterList::iterator center_iter = centerList.begin();
		 center_iter != centerList.end(); ++center_iter)
		if (*center_iter != NULL) {
			
			delete *center_iter;
		}
	// The list must also be cleared to get rid of NULL pointing pointers
	// Failing to do so will lead to Segmentation Faults
	centerList.clear();
}

// --------------------------------------------------------------------------------------

// Draws polygons according to the corners of each voronoi face
void drawBiomes(void) {
	
	CenterList::iterator ic = centerList.begin ();
	while (ic != centerList.end ()) {
		Center* c = *ic++;
		
		if (c->isWater and !c->isBorder) glColor3f (0.0, 0.23, 0.45);
		if (c->isBorder) glColor3f (1.0, 0.23, 0.45);
		if (!c->isBorder and !c->isWater) glColor3f (0.41, 0.58, 0.35);
		
		CornerListObj corners = c->corners;
		
		CornerListObj::iterator corner_iter = corners.begin();
		
		glBegin(GL_POLYGON);
		while (corner_iter != corners.end()) {
			Corner corner = *corner_iter++;
			Point pCorner = corner.point;
			glVertex2f (pCorner[0], pCorner[1]);
		}
		glEnd();
	}
}

// Draw the generated map
// If displayDelaunay is set to true then each voronoi polygon
// will be displayed along with its corresponding delaunay triangulation
void drawMap(void) {
	Delaunay toDrawTriang = ((displayRandom) ? triangRand : triang);
	// Voronoi toDrawVoro = ((displayRandom) ? voronoiRand : voronoi);
	
	if (displayBiome) drawBiomes();
	if (displayVoronoi) drawVoronoiGivenDelaunay(toDrawTriang);
	if (displayDelaunay) drawDelaunayTriang(toDrawTriang);
	glutPostRedisplay();
	
}

// Draws the points in the generated point list
void drawPoints (void) {
	PointList toDrawPts = ((displayRandom) ? pointsRand : points);
	
	/// Draw points
	if (displayBiome) glColor3f (0.0, 0.0, 0.0);
	else glColor3f (0.7, 0.0, 0.0);
	glPointSize (3);
	glBegin(GL_POINTS);
	PointList::iterator ip = toDrawPts.begin ();
	while (ip != toDrawPts.end ()) {
		Point p = *ip++;
		glVertex2f (p [0], p [1]);
	}
	glEnd();
}

void display(void) {
	glClear (GL_COLOR_BUFFER_BIT);
	glClearColor(0.6, 0.6, 0.7, 1.0);
	drawMap();
	if (!displayBiome) drawPoints();
	
	/// Finish
	glutSwapBuffers ();
}

void reshape (int wid, int hgt)
{
	glViewport(0,0,width=wid,height=hgt);
	glMatrixMode (GL_PROJECTION_MATRIX);
	glLoadIdentity ();
	gluOrtho2D (0, width, height, 0);
}

// Performed once at startup and whenever 'r' is pressed
// Generates new voronoi diagrams and relationship graphs
void generateNewGraphs () {
	generateRandPoints();
	performPointDistribution();
	buildRelationshipLists();
	
}

// --------------------------------------------------------------------------------------

// Mouse interactivity
void mouse (int button, int state, int x, int y)
{
	if (state != GLUT_DOWN) return;
	
	// Add new points on left mouse button click
	// Will only be performed if addPointsWithClick is set
	if (button == GLUT_LEFT_BUTTON and addPointsWithClick) {
		Point p ((double) x, (double) y);
		points.push_back (p);
		
		triang.push_back (p);
		glutPostRedisplay ();
	}
	
	if (button == GLUT_MIDDLE_BUTTON) {
		displayDelaunay = not displayDelaunay;
	}
	
	if (button == GLUT_RIGHT_BUTTON) {
		displayVoronoi = not displayVoronoi;
	}
	
	else return;
}

// Keyboard interactivity
void keyboard (unsigned char key, int x, int y) {
	std::string help;
	switch (key) {
			// display Voronoi diagram using only randomized points
		case '1':
			displayRandom = true;
			displayBiome = false;
			displayRadial = false;
			break;
			
			// display Voronoi diagram with applied lloyd relaxation
		case '2':
			displayRandom = false;
			displayBiome = false;
			displayRadial = false;
			break;
			
			// display radial
		case '3':
			displayRandom = false;
			displayBiome = true;
			displayRadial = true;
			buildRelationshipLists();
			break;
			
			// generate new diagrams
		case 'r':
			// Be sure to delete previous pointers in list as new ones will be generated
			resetPointerLists();
			generateNewGraphs();
			break;
			
		case 'h':
			help = "\n 1 - Voronoi diagram with random points \n "
			"2 - Voronoi diagram with applied LLoyd relaxation \n "
			"3 - Randomely generated map \n "
			"4 - Radial map generarion \n "
			"r - generate a new diagram \n "
			"Middle Mouse Button - display Delaunay triangulation \n "
			"Left Mouse Button - display Voronoi diagram \n";
			std::cout << help << std::endl; // flush the stream
			break;
			
		default:
			break;
	}
	
	glutPostRedisplay();
}


// Test related methods
// --------------------------------------------------------------------------------------

void printCornerListInfo () {
	CenterList::iterator ic = centerList.begin ();
	while (ic != centerList.end ()) {
		Center* c = *ic++;
		std::cout << c->point;
	}
	
}

// Prints the endpoint of a given halfedge
// @see: http://i.cs.hku.hk/~wchu/cgal_manual/doc_html/cgal_manual/Voronoi_diagram_2/Chapter_main.html
void print_endpoint(Halfedge_handle e, bool is_src) {
	std::cout << "\t";
	
	if ( is_src ) {
  // print the source vertex of the halfedge
		if ( e->has_source() )  std::cout << e->source()->point() << std::endl;
		else  std::cout << "point at infinity" << std::endl;
	} else {
  // print the end point of the halfedge
		if ( e->has_target() )  std::cout << e->target()->point() << std::endl;
		else  std::cout << "point at infinity" << std::endl;
	}
}

// --------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glClearColor(0.7, 0.7, 0.7, 1.0);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize (width, height);
	glutInitWindowPosition (100, 400);
	glutCreateWindow ("Map Generator");
	generateNewGraphs();
	glutDisplayFunc(display);
	glutMouseFunc (mouse);
	glutKeyboardFunc (keyboard);
	glutReshapeFunc (reshape);
	glutMainLoop();
	return 0;
}
