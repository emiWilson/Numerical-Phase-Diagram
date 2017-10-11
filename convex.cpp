#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <functional>
#include <vector>
#include <set>
#include <stdlib.h>
#include <time.h>
#include <tuple>
#include <deque>

using namespace std;


////Classes for geometry/////

class point {
public:
    double x, y, z;
    point(double x, double y, double z) : x(x), y(y), z(z) {};
};

class gridPoint {
public:
    int x, y;
    double z;
    gridPoint(int x, int y, double z) : x(x), y(y), z(z) {};
    operator point() {return point(x,y,z);};
    bool operator ==(gridPoint p) {return x == p.x && y == p.y && z == p.z;};
};

class edge {
public:
    gridPoint p1, p2;
    tuple<int,int,int,int> sortedXY;
    edge(gridPoint p1, gridPoint p2) : p1(p1), p2(p2) {
        //Implement so that edges can be sorted/tested for equality. We want edges to be considered equal regardless of the order
        //the endpoints are specified, hence getSortedXY, which orders the endpoints into a tuple. z-values are ignored. 
        
        //Sort endpoints into a tuple
        tuple<int,int> t1(p1.x, p1.y), t2(p2.x, p2.y);
    
        //sort the tuples by which is "lower"
        tuple<int,int> tLower, tUpper;
        if(t1 > t2) {
            tLower = t2;
            tUpper = t1;
        } else {
            tLower = t1;
            tUpper = t2;            
        }
    
        sortedXY = tuple<int,int,int,int>(get<0>(tLower), get<1>(tLower), get<0>(tUpper), get<1>(tUpper));
    };
    
};

class triangle {
public:
    gridPoint p1, p2, p3;
    triangle(gridPoint p1, gridPoint p2, gridPoint p3) : p1(p1), p2(p2), p3(p3) {};
    triangle(gridPoint p, edge e) : p1(p), p2(e.p1), p3(e.p2) {};
};


////Required for sorting edges/////

//For use with std::sort. Equivalent edges will be sorted next to each other
bool sortEdges (edge lhs, edge rhs) {
    return lhs.sortedXY < rhs.sortedXY;//leftTuple < rightTuple;
}

//Equality check
bool operator==(const edge &lhs, const edge &rhs) {
    return lhs.sortedXY == rhs.sortedXY;
}






////Printing/////

string printPoint(point p) {
    return "(" + to_string(p.x) + "," + to_string(p.y) + "," + to_string(p.z) + ")";
}

string printPoint(gridPoint p) {
    //return "(" + to_string(p.x) + "," + to_string(p.y) + ",-)";
    return "{" + to_string(p.x) + "," + to_string(p.y) + "," + to_string(p.z) + "}";
}

string printPointNoZ(gridPoint p) {
    //return "(" + to_string(p.x) + "," + to_string(p.y) + ",-)";
    return "{" + to_string(p.x) + "," + to_string(p.y) + "}";
}

string printEdge(edge e) {
    return "[" + printPoint(e.p1) + "," + printPoint(e.p2) + "]";
}

string printEdgeNoZ(edge e) {
    return "{{" + to_string(e.p1.x) + "," + to_string(e.p1.y) + "},{" + to_string(e.p2.x) + "," + to_string(e.p2.y) + "}}";
}

string printTri(triangle tri) {
    return "{" + printPoint(tri.p1) + "," + printPoint(tri.p2) + "," + printPoint(tri.p3) + "}";
}

string printPath(vector<gridPoint> path) {
    string result =  "{";
    for(auto p = path.begin(); p != path.end(); p++) {
        result += printPointNoZ(*p);
    }
    result += "}";
    return result;
}





/////Basic geometry/////

double dotProduct(point v1, point v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double dotProduct2d(point v1, point v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

point crossProduct(point v1, point v2) {
    return point(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

gridPoint subtract(gridPoint p1, gridPoint p2) {
    return gridPoint(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

point subtract(point p1, point p2) {
    return point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}


/////Geometry predicates/////

//is the 2d projection of the triangle onto the xy plane degenerate (a line or point)?
bool isDegenerate2d(triangle tri) {
    gridPoint v1 = subtract(tri.p2, tri.p1);
    gridPoint v2 = subtract(tri.p3, tri.p1);
    return (v1.x * v2.y - v1.y * v2.x == 0);
}

point planeNormalUp(triangle plane) {
    point norm = crossProduct( subtract(plane.p2, plane.p1), subtract(plane.p3, plane.p1));
    if(norm.z >= 0) {
        return norm;
    } else {
        return subtract(point(0,0,0), norm);
    }
}

//Above or on
bool pointAbovePlane(gridPoint p, triangle plane) {
    if(p == plane.p1 || p == plane.p2 || p == plane.p3) {
        return true;
    }
    
    point normal = planeNormalUp(plane);
    
    point v = subtract(p, plane.p1);
    
    double dot = dotProduct(normal, v);
    
    return dot >= 0;
}

bool triangleNotVisible(gridPoint p, triangle tri) {
    point normal = planeNormalUp(tri);
    point vectorToPlane = subtract(tri.p1, p);
    double dot = dotProduct(normal, vectorToPlane);
    return dot <= 0;
}

bool pointAboveHull(gridPoint point, vector<triangle> hull) {
    //assuming the hull is convex, if the point is below any plane it is below the hull:
    for(auto h = hull.begin(); h != hull.end(); h++) {
        if(!pointAbovePlane(point, *h)) {
            return false;
        }
    }
    return true;
}




/////Algorithm implementation/////

//removes the found triangles from the hull
vector<triangle> extractVisibleTriangles(gridPoint p, vector<triangle> &hull) {
    //partition the list so that the visible are at the end and the not visible are at the beginning
    //(it's computationally cheaper to remove from the end of the list):
    vector<triangle>::iterator bound = partition (hull.begin(), hull.end(), bind(triangleNotVisible,p,std::placeholders::_1));
    //copy out visible to be returned:
    vector<triangle> visible(bound, hull.end()); 
    //remove visible from the hull:
    hull.erase(bound, hull.end()); 
    
    return visible;
}

vector<edge> removeDuplicateEdges(vector<edge> edges) {
    //To find the duplicates we sort the list of edges so that 
    //equal edges will be next to each other in the list.
        
    //sort them
    sort(edges.begin(), edges.end(), sortEdges);
    
    //copy out non-duplicates:
    vector<edge> edgesOut;
    unsigned int i = 0;
    while(i + 1 < edges.size()) { //i < size - 1, but in case size is zero (it's unsigned) we use i + 1 < size
        if(edges[i] == edges[i + 1]) {  //it's a duplicate
            i += 2;
        } else { //it's not a duplicate, copy it out:
            edgesOut.push_back(edges[i]); 
            i++;
        }
    }
    
    //if there's still one entry left, it isn't a duplicate:
    if(i + 1 == edges.size()) { 
        edgesOut.push_back(edges[i]);
    }
    
    return edgesOut;
}

vector<edge> getEdgesOfTriangles(vector<triangle> group) {
    //We only want the edges forming the boundray of the triangle group.
    //All interior edges will show up twice, so if we take all edges and 
    //remove the duplicates we will be left with the boundary. 
    
    //find all the edges in the group:
    vector<edge> edges;
    for(auto tri = group.begin(); tri != group.end(); tri++) {
        edges.push_back(edge((*tri).p1, (*tri).p2));
        edges.push_back(edge((*tri).p2, (*tri).p3));
        edges.push_back(edge((*tri).p3, (*tri).p1));
    }
    
    return removeDuplicateEdges(edges);
}

void iterate(vector<gridPoint> &points, vector<triangle> &hull) {
    //select a point at random
    int pointToRemoveIndex = rand() % points.size();
    gridPoint thePoint = points[pointToRemoveIndex];
    
    //remove it from points:
    points[pointToRemoveIndex] = points.back();
    points.pop_back();
    
    if(pointAboveHull(thePoint, hull)) {
        //if it's above the hull, do nothing
        return;
    } else {
        //otherwise we must extend the hull to include the point:
    
        //find and remove all triangles in the hull visible to the point
        vector<triangle> visibleTriangles = extractVisibleTriangles(thePoint, hull);
    
        //find the edges of the visible triangles
        vector<edge> edges = getEdgesOfTriangles(visibleTriangles);
    
        //add new triangles from the point to the edges
        for(auto edge = edges.begin(); edge != edges.end(); ++edge) {
            triangle newTri(thePoint, *edge);
            //remove denerate (vertically oriented) triangles:
            if(!isDegenerate2d(newTri)) {
                hull.push_back(newTri);
            }
        }
    }
}

vector<triangle> getInitialHull(vector<gridPoint> &points, gridPoint topLeft, gridPoint topRight, gridPoint bottomLeft, gridPoint bottomRight) {
    //we take the four corners plus the minimum point in the set to generate the initial convex hull:
    
    gridPoint minPoint = points[0];
    for(auto p = points.begin(); p != points.end(); p++) {
        if((*p).z < minPoint.z) {
            minPoint = *p;
        }
    }
    
    //remove the corners from the points list:
    points.erase(std::remove(points.begin(), points.end(), topLeft), points.end());
    points.erase(std::remove(points.begin(), points.end(), topRight), points.end());
    points.erase(std::remove(points.begin(), points.end(), bottomLeft), points.end());
    points.erase(std::remove(points.begin(), points.end(), bottomRight), points.end());
    
    //remove the minPoint
    points.erase(std::remove(points.begin(), points.end(), minPoint), points.end());
    
    triangle tri1 = triangle(topLeft, topRight, minPoint);
    triangle tri2 = triangle(topRight, bottomRight, minPoint);
    triangle tri3 = triangle(bottomRight, bottomLeft, minPoint);
    triangle tri4 = triangle(bottomLeft, topLeft, minPoint);
    
    vector<triangle> hull;
    //remove denerate (vertically oriented) triangles:
    if(!isDegenerate2d(tri1)){ hull.push_back(tri1); };
    if(!isDegenerate2d(tri2)){ hull.push_back(tri2); };
    if(!isDegenerate2d(tri3)){ hull.push_back(tri3); };
    if(!isDegenerate2d(tri4)){ hull.push_back(tri4); };
    
    return hull;
}

vector<triangle> generateConvexHull(vector<gridPoint> points, gridPoint topLeft, gridPoint topRight, gridPoint bottomLeft, gridPoint bottomRight) {
    vector<triangle> hull = getInitialHull(points, topLeft, topRight, bottomLeft, bottomRight);
    
    while(points.size() > 0) {
        cout << "points remaining " << points.size() << "\n";
        iterate(points, hull);
    }
    
    return hull;
}

vector<triangle> generateConvexHullFromData(int width, int height, vector<double> data) {
    vector<gridPoint> points;
    int dataIndex = 0;
    
    //column major order
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            points.push_back(gridPoint(i, j, data[dataIndex]));
            dataIndex++;
        }
    }
    
    vector<triangle> hull = generateConvexHull(points, points[0], points[height - 1], points[height * (width - 1)], points[width * height - 1]);
    return hull;
}




/////Using the convex hull to get output data/////

//Sets the boolean array isOnHull to true or false for each point.
//A point is on the hull if it is the endpoint of one of the hull triangles.
void setIsOnHull(bool* isOnHull, int width, int height, vector<triangle> hull) {
    for(int i = 0; i < width * height; i++) {
        isOnHull[i] = false;
    }
    
    for(auto tri = hull.begin(); tri != hull.end(); tri++) {
        triangle t = *tri;
        isOnHull[t.p1.y * width + t.p1.x] = true;
        isOnHull[t.p2.y * width + t.p2.x] = true;
        isOnHull[t.p3.y * width + t.p3.x] = true;
    }
}

//Determines if a point is on the edge of the domain
bool notOnEdge(gridPoint p, int width, int height) {
    if(p.x == 0 || p.y == 0 || p.x == width - 1 || p.y == height - 1) {
        return false;
    } else {
        return true;
    }
}

//Removes edges that lie entirely on the edge of the domain
vector<edge> removeAreaEdges(vector<edge> edges, int width, int height) {
    vector<edge> edgesOut;
    for(auto e = edges.begin(); e != edges.end(); e++) {
        if(notOnEdge((*e).p1, width, height) || notOnEdge((*e).p2, width, height)) {
            edgesOut.push_back(*e);
        }
    }
    return edgesOut;
}

//"large" triangles are those of area greater than 1/2. The convex areas of the input
//values will be tiled by "small" triangles of area 1/2. Thus the "large" triangles
//outline the areas of coexistence.
vector<triangle> getLargeTriangles(vector<triangle> tris) {
    vector<triangle> outTris;
    for(auto t = tris.begin(); t != tris.end(); t++) {
        triangle tri = *t;
        gridPoint u = subtract(tri.p3, tri.p1);
        gridPoint v = subtract(tri.p2, tri.p1);
        int twiceArea2d = abs(u.x * v.y - u.y * v.x);
        if(twiceArea2d != 1) {
            outTris.push_back(tri);
        }
    }
    return outTris;
}

//Coexistence lines are determined by assuming the triangles are relatively long. 
//We take the shortest edge and find its midpoint, and then return an edge from that 
//midpoint to the point opposite the short side.
tuple<double, double, double, double> getCoexistLine(triangle tri) {
    //vectors for each edge:
    gridPoint u = subtract(tri.p3, tri.p1);
    gridPoint v = subtract(tri.p2, tri.p1);
    gridPoint w = subtract(tri.p3, tri.p2);
    
    //norm^2 for each edge:
    int uu = dotProduct2d(u, u);
    int vv = dotProduct2d(v, v);
    int ww = dotProduct2d(w, w);
    
    //find the edge:
    if(uu <= vv && uu <= ww) { //u is the short edge
        return tuple<double, double, double, double>((tri.p3.x + tri.p1.x) / 2.0, (tri.p3.y + tri.p1.y) / 2.0, tri.p2.x, tri.p2.y);
    }
    if(vv <= uu && vv <= ww) { //v is the short edge
        return tuple<double, double, double, double>((tri.p2.x + tri.p1.x) / 2.0, (tri.p2.y + tri.p1.y) / 2.0, tri.p3.x, tri.p3.y);
    }
    if(ww <= uu && ww <= vv) { //w is the short edge
        return tuple<double, double, double, double>((tri.p3.x + tri.p2.x) / 2.0, (tri.p3.y + tri.p2.y) / 2.0, tri.p1.x, tri.p1.y);
    }
    
    //should never get here
    throw "No shortest edge";
}

//Connects a set of edges into paths. Paths are sequences of points, an edge is a path of length 2.
//Sequences of edges sharing common endpoints are converted into equivalent paths.
vector<vector<gridPoint>> connectEdges(vector<edge> edges) {
    deque<vector<gridPoint>> paths;
    
    //Push all edges into the list of paths:
    for(auto e = edges.begin(); e != edges.end(); e++) {
        vector<gridPoint> path;
        path.push_back((*e).p1);
        path.push_back((*e).p2);
        paths.push_back(path);
    }
    
    //It will take at most this many iterations to connect all paths:
    int iterations = paths.size();
    
    while(iterations >= 0) {
        iterations--;
        
        //We take the first path, then search for any paths that connect to it at the front or back.
        //If we find one we merge the paths into one.
        //The resulting path (whether or not it was merged with another one) is pushed to the back of
        //the queue.
        
        vector<gridPoint> path = paths[0];
        paths.pop_front();
        
        //search for connections and merge:
        for(int i = 0; i < (int)(paths.size()); i++) {
            gridPoint first = path[0];
            gridPoint last  = path.back();
            
            //front to back
            if(paths[i][0] == last) {
                auto foundPath = paths[i];
                
                //copy found path to end of path
                path.insert(path.end(), foundPath.begin() + 1, foundPath.end());
                //remove found path
                paths.erase(paths.begin() + i);
                i--;
                continue;
            }
            
            //back to front
            if(paths[i].back() == first) {
                auto foundPath = paths[i];
                
                //copy path to end of found path:
                foundPath.insert(foundPath.end(), path.begin() + 1, path.end());
                //this is now our path:
                path = foundPath;
                //remove the found path from the list:
                paths.erase(paths.begin() + i);
                i--;
                continue;
            }
            
            //front to front
            if(paths[i][0] == first) {
                auto foundPath = paths[i];

                //reverse the path:
                reverse(path.begin(), path.end());
                //now append the found path to it:
                path.insert(path.end(), foundPath.begin() + 1, foundPath.end());
                //remove found path
                paths.erase(paths.begin() + i);
                i--;
                continue;
            }

            //back to back
            if(paths[i].back() == last) {
                auto foundPath = paths[i];
                
                //reverse the found path:
                reverse(foundPath.begin(), foundPath.end());
                //now append the found path to path:
                path.insert(path.end(), foundPath.begin() + 1, foundPath.end());
                //remove found path
                paths.erase(paths.begin() + i);
                i--;
                continue;
            }
        }
        
        paths.push_back(path);
    }
    
    return vector<vector<gridPoint>>(paths.begin(), paths.end());
}



/////writing out data for the hull/////

//Writes out a grid of 1s and 0s indicating which points in the input were on the resulting
//convect hull. (On hull is "1", not on hull is "0").
void writeOnHull(string filename, vector<triangle> hull, int width, int height) {
    ofstream outfile( filename.c_str(), ofstream::out );
    
    bool onHull[width * height];
    
    setIsOnHull(onHull, width, height, hull);
        
    int idx = 0;
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            outfile << (onHull[idx] ? "1" : "0");
            idx++;
        }
        outfile << "\n";
    }
    
}

//Outputs all triangles.
void writeTris(string filename, vector<triangle> tris, double xMin, double yMin, double xStep, double yStep) {
    ofstream outfile( filename.c_str(), ofstream::out );
    
    for(unsigned int i = 0; i < tris.size(); i++) {
        auto tri = tris[i];
        outfile << "{{"  << tri.p1.x * xStep + xMin << "," << tri.p1.y * yStep + yMin << "," << tri.p1.z
                << "},{" << tri.p2.x * xStep + xMin << "," << tri.p2.y * yStep + yMin << "," << tri.p2.z
                << "},{" << tri.p3.x * xStep + xMin << "," << tri.p3.y * yStep + yMin << "," << tri.p3.z << "}}";
        
        if(i + 1 < tris.size()) {
            outfile << ",\n";
        }
        
    }
}

//Outputs the coexisting lines (as determined by getCoexistLine) of the given convex hull.
void writeCoexistLines(string filename, vector<triangle> hull, double xMin, double yMin, double xStep, double yStep) {
    auto largeTris = getLargeTriangles(hull);
    
    vector<tuple<double, double, double, double>> edges;
    for(auto tri = largeTris.begin(); tri != largeTris.end(); tri++) {
        edges.push_back(getCoexistLine(*tri));
    }
    
    ofstream outfile( filename.c_str(), ofstream::out );
    
    for(unsigned int i = 0; i < edges.size(); i++) {
        outfile << "{{" 
            << get<0>(edges[i]) * xStep + xMin << "," 
            << get<1>(edges[i]) * yStep + yMin << "},{" 
            << get<2>(edges[i]) * xStep + xMin << "," 
            << get<3>(edges[i]) * yStep + yMin << "}}";
        
        if(i + 1 < edges.size()) {
            outfile << ",\n";
        }
    }
    
}

//Writes the paths outlining the coexistence region
void writeCoexistOutline(string filename, vector<triangle> hull, int width, int height, double xMin, double yMin, double xStep, double yStep) {
    ofstream outfile( filename.c_str(), ofstream::out );
    
    
    //Use setIsOnHull to find the areas of coexistence
    bool onHull[width * height];
    
    setIsOnHull(onHull, width, height, hull);
    
    vector<edge> edges;
    
    //For each point not on the hull, add the four edges corresponding to a square surrounding that point:
    int idx = 0;
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            if(!onHull[idx]) {
                //0 for the z values, it won't be printed anyway
                edges.push_back(edge(gridPoint(i    , j    , 0), gridPoint(i + 1, j    , 0)));
                edges.push_back(edge(gridPoint(i    , j    , 0), gridPoint(i    , j + 1, 0)));
                edges.push_back(edge(gridPoint(i + 1, j    , 0), gridPoint(i + 1, j + 1, 0)));
                edges.push_back(edge(gridPoint(i    , j + 1, 0), gridPoint(i + 1, j + 1, 0)));
            }
            idx++;
        }
    }
    
    //Remove all duplicate edges and edges on the perimeter of the domain, then connect them into paths:
    auto paths = connectEdges(removeAreaEdges(removeDuplicateEdges(edges), width, height));
    
    //Print the paths:
    for(unsigned int i = 0; i < paths.size(); i++) {
        outfile << "{";
        for(unsigned int j = 0; j < paths[i].size(); j++) {
            //shift the point by -.5, -.5 to account for placment of squares at (i,j) to (i+1,j+1) above
            double x = (paths[i][j].x - .5) * xStep + xMin; 
            double y = (paths[i][j].y - .5) * yStep + yMin;
            outfile << "{" << x << "," << y << "}";
            if(j + 1 < paths[i].size()) {
                outfile << ",";
            }
        }
        outfile << "}";
        if(i + 1 < paths.size()) {
            outfile << ",";
        }
        outfile << "\n";
    }
}



/////tests/////

void test_connectEdges() {
    vector<edge> edges;

    edges.push_back(edge(gridPoint(1,2,3.0), gridPoint(1,3,3.0)));
    edges.push_back(edge(gridPoint(7,6,3.0), gridPoint(1,3,3.0)));

    edges.push_back(edge(gridPoint(5,6,3.0), gridPoint(2,9,3.0)));
    edges.push_back(edge(gridPoint(7,9,3.0), gridPoint(5,6,3.0)));

    edges.push_back(edge(gridPoint(7,3,3.0), gridPoint(8,3,3.0)));
    edges.push_back(edge(gridPoint(7,3,3.0), gridPoint(2,6,3.0)));


    edges.push_back(edge(gridPoint(3,8,3.0), gridPoint(7,7,3.0)));
    edges.push_back(edge(gridPoint(7,7,3.0), gridPoint(7,9,3.0)));


    edges.push_back(edge(gridPoint(1,1,3.0), gridPoint(5,5,3.0)));
    edges.push_back(edge(gridPoint(5,5,3.0), gridPoint(4,3,3.0)));
    edges.push_back(edge(gridPoint(4,3,3.0), gridPoint(1,1,3.0)));
    
    edges.push_back(edge(gridPoint(23,22,0.0), gridPoint(24,22,0.0)));
    edges.push_back(edge(gridPoint(24,22,0.0), gridPoint(25,22,0.0)));
    edges.push_back(edge(gridPoint(25,21,0.0), gridPoint(25,22,0.0)));
    edges.push_back(edge(gridPoint(25,21,0.0), gridPoint(26,21,0.0)));
    edges.push_back(edge(gridPoint(26,21,0.0), gridPoint(27,21,0.0)));
    
    auto paths = connectEdges(edges);
    
    for(unsigned int i = 0; i < paths.size(); i++) {
        cout << "{";
        for(unsigned int j = 0; j < paths[i].size(); j++) {
            cout << printPointNoZ(paths[i][j]);
        }
        cout << "}\n";
    }
    
}

void test_getEdgesOfTriangles() {
    cout << "test_getEdgesOfTriangles\n";
    gridPoint a(0,0,1), b(0,1,2), c(1,0,3), d(1,1,4), e(2,0,5);
    vector<triangle> tris1, tris2, tris3, tris4, tris5;
    tris1.push_back(triangle(a,b,c));
    tris1.push_back(triangle(b,c,d));
    
    //different ordering of shared edge:
    tris2.push_back(triangle(a,b,c));
    tris2.push_back(triangle(c,b,d));

    tris3.push_back(triangle(a,b,c));
    
    //tris4 is empty
    
    tris5.push_back(triangle(a,b,c));
    tris5.push_back(triangle(c,b,d));
    tris5.push_back(triangle(b,d,e));
    
    auto edges1 = getEdgesOfTriangles(tris1);
    auto edges2 = getEdgesOfTriangles(tris2);
    auto edges3 = getEdgesOfTriangles(tris3);
    auto edges4 = getEdgesOfTriangles(tris4);
    auto edges5 = getEdgesOfTriangles(tris5);
    
    if(edges1.size() != 4 || edges2.size() != 4 || edges3.size() != 3 || edges4.size() != 0 || edges5.size() != 5) {
        cout << edges1.size() << "\n";
        cout << edges2.size() << "\n";
        cout << edges3.size() << "\n";
        cout << edges4.size() << "\n";
        cout << edges5.size() << "\n";
        cout << "Unexpected number of edges\n";
    }
}

void test_isDegenerate2d() {
    cout << "test_isDegenerate2d\n";
    for(int i = 0; i < 100; i++) {
        int rx = rand() % 5;
        int ry = rand() % 5;
        int a = rand() % 5;
        int b = rand() % 5;
        int c = rand() % 5;
        gridPoint p1(a * rx, a * ry, rand()%100);
        gridPoint p2(b * rx, b * ry, rand()%100);
        gridPoint p3(c * rx, c * ry, rand()%100);
        triangle tri = triangle(p1, p2, p3);
        if(!isDegenerate2d(tri)) {
            cout << "test " << i;
            cout << " Unexpected: " << printTri(tri);
            cout << "\n";
        }
    }
}



/////input/////
vector<string> readFile(string filename) {
    vector <string> data;
    ifstream infile( filename.c_str() );

    while (infile)
    {
        string line;
        if (!getline( infile, line )) 
            break;

        istringstream ss( line );

        while (ss)
        {
            string s;
            if (!getline( ss, s, ',' )) 
                break;
            data.push_back( s );
        }
    }
    
    if (!infile.eof()) {
        cerr << "??\n";
    }
    
    return data;
}

std::string trim(const std::string& str, const std::string& whitespace = " \t")
{
    const unsigned int strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const int strEnd = str.find_last_not_of(whitespace);
    const int strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

template <class T> T string_convert(const std::string& s) {
	std::istringstream i(s);
	T x;
	if (!(i >> x)) {
		//return NAN;
	}
   return x;
}

vector<double> convertStringData(vector<string> dataIn) {
    vector<double> dataOut;
    for(auto datum = dataIn.begin(); datum != dataIn.end(); datum++) {
        double dataNum = string_convert<double>(trim(*datum));
        dataOut.push_back(dataNum);
    }
    return dataOut;
}


void run(vector<string> args) {
    
    //Check that the user had provided the correct number of arguments:
    int numArgs = args.size();
    
    if(numArgs != 8) {
        cout << "Arguments: [size x] [size y] [file for data values] [output prefix] [xmin] [ymin] [xstep] [ystep]\n";
        return;
    }
    
    //Convert input arguments:
    unsigned int width  = string_convert<unsigned int>(trim(args[0]));
    unsigned int height = string_convert<unsigned int>(trim(args[1]));
    string inputFileName = args[2];
    string outPrefix = args[3];
    double xmin  = string_convert<double>(trim(args[4]));
    double ymin  = string_convert<double>(trim(args[5]));
    double xstep = string_convert<double>(trim(args[6]));
    double ystep = string_convert<double>(trim(args[7]));
    
    //Get the input data in string format:
    vector<string> stringData = readFile(inputFileName);
    
    //Check that the data size matches the expected size from the width and height:
    if(stringData.size() != width * height) {
        cout << "Data size (" << stringData.size() << ") doesn't match expected size (" << width << " * " << height << ")\n";
        return;
    }
    
    //Convert the string data into doubles:
    vector<double> data = convertStringData(stringData);
    
    //Generate the convex hull:
    auto hull = generateConvexHullFromData(width, height, data);
    
    //Output data:
    writeTris(outPrefix + "_hull_tris.txt", hull, xmin, ymin, xstep, ystep);
    
    writeOnHull(outPrefix + "_on_hull.txt", hull, width, height);
    
    writeCoexistOutline(outPrefix + "_coexist_region.txt", hull, width, height, xmin, ymin, xstep, ystep);
    
    writeCoexistLines(outPrefix + "_coexist_lines.txt", hull, xmin, ymin, xstep, ystep);
}

int main(int argc, char *argv[]) {
    srand (time(NULL));
    
    int numArgs = argc - 1;
    vector<string> args;
    if (numArgs > 0) {
        args.assign(argv + 1, argv + argc);
    }

    run(args);
}