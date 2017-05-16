#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <functional>
#include <vector>
#include <set>
#include "stringConvert.h"
#include <stdlib.h>
#include <time.h>
#include <tuple>


using namespace std;


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
    edge(gridPoint p1, gridPoint p2) : p1(p1), p2(p2) {};
    
};





class triangle {
public:
    gridPoint p1, p2, p3;
    triangle(gridPoint p1, gridPoint p2, gridPoint p3) : p1(p1), p2(p2), p3(p3) {};
    triangle(gridPoint p, edge e) : p1(p), p2(e.p1), p3(e.p2) {};
};




//Implement so that edges can be sorted/tested for equality. We want edges to be considered equal regardless of the order
//the endpoints are specified, hence getSortedXY, which orders the endpoints into a tuple. z-values are ignored. 

//Sort endpoints into a tuple
tuple<int,int,int,int> getEdgeSortedXY(edge e) {
    tuple<int,int> t1(e.p1.x, e.p1.y), t2(e.p2.x, e.p2.y);
    
    //sort the tuples by which is "lower"
    tuple<int,int> tLower, tUpper;
    if(t1 > t2) {
        tLower = t2;
        tUpper = t1;
    } else {
        tLower = t1;
        tUpper = t2;            
    }
    
    return tuple<int,int,int,int>(get<0>(tLower), get<1>(tLower), get<0>(tUpper), get<1>(tUpper));
}

//For use with std::sort. Equivalent edges will be sorted next to each other
bool sortEdges (edge lhs, edge rhs) {
    auto leftTuple  = getEdgeSortedXY(lhs);
    auto rightTuple = getEdgeSortedXY(rhs);
    
    return leftTuple < rightTuple;
}

//Equality check
bool operator==(const edge &lhs, const edge &rhs) {
    return (!sortEdges(lhs, rhs)) && (!sortEdges(rhs, lhs));
}

//get a triangle's edges
void addTriEdges(vector<edge> &edges, triangle tri) {
    edges.push_back(edge(tri.p1, tri.p2));
    edges.push_back(edge(tri.p2, tri.p3));
    edges.push_back(edge(tri.p3, tri.p1));
}






string printPoint(point p) {
    return "(" + to_string(p.x) + "," + to_string(p.y) + "," + to_string(p.z) + ")";
}

string printPoint(gridPoint p) {
    //return "(" + to_string(p.x) + "," + to_string(p.y) + ",-)";
    return "{" + to_string(p.x) + "," + to_string(p.y) + "," + to_string(p.z) + "}";
}

string printEdge(edge e) {
    return "[" + printPoint(e.p1) + "," + printPoint(e.p2) + "]";
}

string printTri(triangle tri) {
    return "{" + printPoint(tri.p1) + "," + printPoint(tri.p2) + "," + printPoint(tri.p3) + "}";
}

void printTris(vector<triangle> tris) {
    cout << "Size " << tris.size() << ":\n";
    for(int i = 0; i < tris.size(); i++) {
        cout << "  " << printTri(tris[i]) << ",\n";
    }
}








double dotProduct(point v1, point v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
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

void check_convex(vector<triangle> hull) {
    for(int i = 0; i < hull.size(); i++) {
        triangle tri = hull[i];
        for(int j = 0; j < hull.size(); j++) {
            if(j != i) {
                if(!pointAbovePlane(hull[j].p1, tri) || !pointAbovePlane(hull[j].p2, tri) || !pointAbovePlane(hull[j].p3, tri)) {
                    cout << "not convex!\n";
                    return;
                }   
            }
        }
    }
}




bool triangleNotVisible(gridPoint p, triangle tri) {
    point normal = planeNormalUp(tri);
    point vectorToPlane = subtract(tri.p1, p);
    double dot = dotProduct(normal, vectorToPlane);
    return dot <= 0;
}

void check_visible(gridPoint p, vector<triangle> visible, vector<triangle> notVisible) {
    for(int i = 0; i < visible.size(); i++) {
        if(triangleNotVisible(p, visible[i])) {
            cout << "'visible' triangle is not visible\n";
        }
    }
    for(int i = 0; i < notVisible.size(); i++) {
        if(!triangleNotVisible(p, notVisible[i])) {
            cout << "'not visible' triangle is visible\n";
        }
    }
}

bool pointAboveHull(gridPoint point, vector<triangle> hull) {
    for(int i = 0; i < hull.size(); i++) {
        if(!pointAbovePlane(point, hull[i])) {
            return false;
        }
    }
    return true;
}

//removes the found triangles from the hull
vector<triangle> extractVisibleTriangles(gridPoint p, vector<triangle> &hull) {
    using namespace std::placeholders;
    vector<triangle>::iterator bound = partition (hull.begin(), hull.end(), bind(triangleNotVisible,p,_1));
    vector<triangle> visible(bound, hull.end()); //visible are at the end of the partition
    hull.erase(bound, hull.end()); //remove visible
    return visible;
}

vector<edge> getEdgesOfTriangles(vector<triangle> group) {
    vector<edge> allEdges;
    for(int i = 0; i < group.size(); i++) {
        addTriEdges(allEdges, group[i]);
    }
    sort(allEdges.begin(), allEdges.end(), sortEdges);
    
    vector<edge> edgesOut;
    int i = 0;
    while(i + 1 < allEdges.size()) { //in case size() is zero (it's unsigned)
        if(allEdges[i] == allEdges[i + 1]) {
            i += 2;
        } else {
            edgesOut.push_back(allEdges[i]);
            i++;
        }
    }
    
    if(i == allEdges.size() - 1) {
        edgesOut.push_back(allEdges[i]);
    }
    
    return edgesOut;
}


void iterate(vector<gridPoint> &points, vector<triangle> &hull) {
    //select a point at random
    int pointToRemoveIndex = rand() % points.size();
    gridPoint thePoint = points[pointToRemoveIndex];
    //remove it from points:
    points[pointToRemoveIndex] = points.back();
    points.pop_back();
    
    //if it's above the hull, do nothing
    if(pointAboveHull(thePoint, hull)) {
        return;
    }
    
    //otherwise we must extend the hull to include the point:
    
    //find and remove all visible triangles in the hull
    vector<triangle> visibleTriangles = extractVisibleTriangles(thePoint, hull);
    
    check_visible(thePoint, visibleTriangles, hull);
    
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

vector<triangle> getInitialHull(vector<gridPoint> &points, gridPoint topLeft, gridPoint topRight, gridPoint bottomLeft, gridPoint bottomRight) {
    gridPoint minPoint = points[0];
    for(int i = 0; i < points.size(); i++) {
        if(points[i].z < minPoint.z) {
            minPoint = points[i];
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
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            points.push_back(gridPoint(i, j, data[dataIndex]));
            cout << data[dataIndex] << ",";
            dataIndex++;
        }
    }
    cout << "\n";
    cout << printPoint((point)points[0]) << printPoint((point)points[width - 1]) << printPoint((point)points[width * (height - 1)]) << printPoint((point)points[width * height - 1]) << "\n";
    
    vector<triangle> hull = generateConvexHull(points, points[0], points[width - 1], points[width * (height - 1)], points[width * height - 1]);
    
    return hull;
}

void setIsOnHull(bool* isOnHull, int width, int height, vector<triangle> hull) {
    for(int i = 0; i < width * height; i++) {
        isOnHull[i] = false;
    }
    
    for(int i = 0; i < hull.size(); i++) {
        triangle t = hull[i];
        isOnHull[t.p1.y * width + t.p1.x] = true;
        isOnHull[t.p2.y * width + t.p2.x] = true;
        isOnHull[t.p3.y * width + t.p3.x] = true;
    }
}




void test_addTriEdges() {
    cout << "test_addTriEdges\n";
    
    gridPoint a(0,0,1), b(0,1,2), c(1,0,3), d(1,1,4);
    triangle tri1(a,b,c);
    triangle tri2(b,c,d);
    
    vector<edge> edges;
    addTriEdges(edges, tri1);
    
    if(edges.size() != 3) {
        cout << "Unexpected " << edges.size() << "\n";
    }
    
    addTriEdges(edges, tri2);
    if(edges.size() != 6) {
        cout << "Unexpected " << edges.size() << "\n";
        for(auto edge = edges.begin(); edge != edges.end(); ++edge) {
            cout << printEdge(*edge);
        }
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

vector<double> convertStringData(vector<string> dataIn) {
    vector<double> dataOut;
    for(int j = 0; j < dataIn.size(); j++) {
        double dataNum = string_convert<double>(trim(dataIn[j]));
        dataOut.push_back(dataNum);
    }
    return dataOut;
}

void run(vector<string> args) {
    int numArgs = args.size();
    
    if(numArgs != 3) {
        cout << "Arguments: [size x] [size y] [file for data values] [file for on_hull]\n";
        return;
    }
    
    int width  = string_convert<double>(trim(args[0]));
    int height = string_convert<double>(trim(args[1]));
    string inputFileName = args[2];
    
    vector<string> stringData = readFile(inputFileName);
    
    if(stringData.size() != width * height) {
        cout << "Data size (" << stringData.size() << ") doesn't match expected size (" << width << " * " << height << ")\n";
        return;
    }
    
    vector<double> data = convertStringData(stringData);
    
    auto hull = generateConvexHullFromData(width, height, data);
    
    printTris(hull);
    
    bool onHull[width * height];
    
    setIsOnHull(onHull, width, height, hull);
        
    int idx = 0;
    for(int j = 0; j < height; j++) {
        for(int i = 0; i < width; i++) {
            cout << onHull[idx] ? "1" : "0";
            idx++;
        }
        cout << "\n";
    }
    
    check_convex(hull);
}

int main(int argc, char *argv[]) {
    srand (time(NULL));
    
    int numArgs = argc - 1;
    vector<string> args;
    if (numArgs > 0) {
        args.assign(argv + 1, argv + argc);
    }

    test_isDegenerate2d();        
    test_getEdgesOfTriangles();
    test_addTriEdges();
    
    run(args);
}