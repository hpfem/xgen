# include <stdlib.h>
# include <stdio.h>
# include <math.h>
#include <ios>      // Required for streamsize
#include <iostream>
#include <istream>
#include <limits>   // Required for numeric_limits

struct Point {
  double x, y;

  Point(double ix, double iy): x(ix), y(iy) {}
  Point() {}

  double abs() {return sqrt(x*x + y*y);}
  int operator == (Point &B) {return x==B.x && y==B.y;}

  Point &operator += (Point B) {x += B.x; y += B.y; return *this;}
  Point &operator -= (Point B) {x -= B.x; y -= B.y; return *this;}
  Point &operator *= (double k) {x *= k; y *= k; return *this;}
  Point &operator /= (double k) {x /= k; y /= k; return *this;}
  Point &operator *= (int k) {x *= k; y *= k; return *this;}
  Point &operator /= (int k) {x /= k; y /= k; return *this;}
  Point operator + (Point B) {return Point(x + B.x, y + B.y);}
  Point operator - (Point B) {return Point(x - B.x, y - B.y);}
  Point operator * (double k) {return Point(k*x, k*y);}
  Point operator / (double k) {return Point(x/k, y/k);}
  Point operator * (int k) {return Point(k*x, k*y);}
  Point operator / (int k) {return Point(x/k, y/k);}
  double operator * (Point B) {return x*B.x + y*B.y;}
};

struct Element {
  int a, b, c;                              // integer indices in the array of points
  double ab_angle, ac_angle, bc_angle;      // central angles of edges
                                            // this angle is zero when an edge is straight
  Element() {}
  Element(int in1, int in2, int in3) {
    a = in1; b = in2; b = in3;
    ab_angle = ac_angle = bc_angle = 0;
  } 
  Element(int in1, int in2, int in3, double ab, double ac, double bc) {
    a = in1; b = in2; c = in3;
    ab_angle = ab; ac_angle = ac; bc_angle = bc;
  } 
};

struct BoundarySegment {
  Point P_start;       // Starting point.
  Point P_end;         // End point.
  int subdiv;          // Equidistant subdivision.
  double alpha;        // 0 for straight edges, nonzero for circular arcs.
  int marker;          // Boundary marker.
  BoundarySegment *next;
  BoundarySegment(int m, double x, double y, int subdivision, double alp) 
  {
    marker = m;
    P_start.x = x;
    P_start.y = y;
    P_end.x = -9999;
    P_end.y = -9999;
    subdiv = subdivision;
    alpha = alp;
    next = NULL;
  }
  double CalculateLength();
};

struct BoundaryEdge {
  int A, B;            // Vertex indices.
  double alpha;        // 0 for straight edges, nonzero for circular arcs.
  int component_index; // Index of boundary component (closed loop).
  int marker;          // Boundary marker.
  BoundaryEdge *next;
  BoundaryEdge() {};
  BoundaryEdge(int a, int b, int ci, int m, double alp) 
  {
    A = a; B = b; 
    component_index = ci;
    marker = m;
    alpha = alp;
    next = NULL;
  }
};

// Pair of interest indices representing a (possibly curvilinear) 
// boundary edge.
struct BoundaryPair {
  int A, B;
  double alpha;
  BoundaryPair *next;
  BoundaryPair(int a, int b, double alp) {A = a; B = b; alpha = alp, next = NULL;}
};

// Represents the boundary (algebraically) via pairs of integer indices
// of grid points. Changes during meshing.
struct BoundaryPairsList {
  BoundaryPair *First, *Last, *Ptr;
  BoundaryPairsList() {First = Last = Ptr = NULL;}
  void Init() {Ptr = First;}
  bool GetNext(int &a, int &b, double &alpha);
  bool Contains(int C);
  void Append(int a, int b, double alp);
  bool Delete(int a, int b, double &alp);
  void DeleteLast();
  void DeleteAll();
  void GetLast(int &a, int &b, double &alp) {a = Last->A; b = Last->B; alp = Last->alpha;}
  void Change(int a, int b, int c, int d, double angle);
  void InsertAfter(int a, int b, int c, int d, double angle);
  bool IsEmpty() {return (First == NULL) ? true : false;}
};

// Pair of points representing a boundary line. Note that multiple
// boundary lines are used to represent a curvilinear boundary edge.
struct BoundaryLine {
  Point A, B;
  BoundaryLine *next;
  BoundaryLine(Point a, Point b) {A = a; B = b; next = NULL;}
};

// Represents boundary (geometrially) via pairs of grid points
// (their physical coordinates). Used for domain area calculation, 
// determination whether a point is inside or outside of the domain,
// etc. Does not change during meshing.
struct BoundaryLinesList {
  BoundaryLine *First, *Last, *Ptr;
  BoundaryLinesList() {First = Last = Ptr = NULL;}
  void Init() {Ptr = First;}
  bool GetNext(Point &a, Point &b);
  void DeleteLast();
  void AddStraightLine(Point a, Point b);
  void AddCircularArc(Point a, Point b, double angle, int subdiv);
  void DeleteAll();
};

// Serves for output to mesh file and for checking whether 
// an edge liea on the boundary. Does not change during 
// meshing.
struct BoundaryEdgeList {
  BoundaryEdge *First, *Last, *Ptr;
  BoundaryEdgeList() {First = Last = Ptr = NULL;}
  void Init() {Ptr = First;}
  bool GetNext(BoundaryEdge &I);
  void Append(int a, int b, int component_index, int marker, double alpha);
  void DeleteLast();
  void DeleteAll();
};

struct ElemBox {
  Element E;
  ElemBox *next;
  ElemBox(int a, int b, int c) {
    E.a = a; E.b = b; E.c = c; 
    E.ab_angle = 0; E.ac_angle = 0; E.bc_angle = 0;
    next = NULL;
  }
  ElemBox(int a, int b, int c, double ab, double ac, double bc) {
    E.a = a; E.b = b; E.c = c; 
    E.ab_angle = ab; E.ac_angle = ac; E.bc_angle = bc;
    next = NULL;
  }
};

struct ElemList {
  ElemBox *First, *Ptr, *Last;
  ElemList() {First = Last = Ptr = NULL;}
  void Init() {Ptr = First;}
  bool GetNext(int &a, int &b, int &c, double &ab, double &ac, double &bc);
  void Append(int a, int b, int c, double ab, double ac, double bc);
  void DeleteAll();
  ~ElemList();
};

struct PointBox {
  Point P;
  PointBox *next;
  PointBox(double pos_x, double pos_y) : P(pos_x, pos_y) {next = NULL;}
};

struct PointList {
  PointBox *First, *Last, *Ptr;
  PointList() {First = Last = Ptr = NULL;}
  void Init() {Ptr = First;}
  bool GetNext(Point &T);
  bool GetNext(int i, double &pos_x, double &pos_y);
  void Append(Point P);
  void Append(double pos_x, double pos_y);
  Point DeleteAll();
};  

struct BoundaryComponent {
  BoundarySegment *First_segment, *Last_segment;  
  BoundaryComponent *next;

  BoundaryComponent() {
    First_segment = Last_segment = NULL;
    next = NULL;
  }
};

// Stores all segments of the domain's boundary as they come from the user.
// Can create a list of all boundary vertices (points), list of boundary 
// lines (for plotting, area calculation, determination whether a point is 
// inside or outside, etc), list of boundary pairs (for meshing), and list 
// of boundary edges for output to mesh file.
struct BoundaryType {
  BoundaryComponent *First_bdy_component, *Last_bdy_component;
  BoundaryType() {First_bdy_component = Last_bdy_component = NULL; Length = -1;}  
  void CreateNewBoundaryComponent();
  void AddBoundarySegment(int marker, double x, double y, int subdiv, double alpha);
  BoundaryPairsList *CreateBoundaryPairsList();
  BoundaryLinesList *CreateBoundaryLinesList();
  PointList *CreateBoundaryPointList();
  BoundaryEdgeList *CreateBoundaryEdgeList();
  void CalculateLength();
  double GiveLength() {return Length;};
  ~BoundaryType();
  double Length;
};

// Macro sqr(..):
# define sqr(a) ((a)*(a))

// Errors & warnings:

void XgError(char *who, char *what);
void XgError(const char *who, const char *what);
void XgError(char *what);
void XgError(const char *what);
void XgWarning(char *who, char *what);
void XgMessage(char *what);

class Xgen {
  public:
  Xgen(bool nogui, int steps_to_take, bool overlay);
  char* XgGiveName();
  double XgGiveTimestep();
  int XgGiveDimension();
  long XgGiveNpoin();
  long XgGiveInteriorPtsNum();
  long XgGiveBoundaryPtsNum();
  long XgGiveNelem();
  void XgGiveLimits(Point &min, Point &max);
  void XgInitPointList();
  bool XgGiveNextPoint(Point &p);
  void XgInitInteriorPointList();
  bool XgGiveNextInteriorPoint(Point &p);
  void XgInitBoundaryPairsList();
  void XgInitBoundaryLinesList();
  bool XgGiveNextBoundaryPair(Point &p, Point &q, double &alpha);
  bool XgGiveNextBoundaryLine(Point &p, Point &q);
  void XgInitBoundaryEdgeList();
  bool XgIsBoundaryEdge(int a, int b);
  double XgGetBdyEdgeAngle(int A, int B); // returns zero if AB is not a boundary edge
  bool XgGiveNextBoundaryEdge(BoundaryEdge &edge_ptr);
  void XgInitElementList();
  bool XgGiveNextElement(Point &p, Point &q, Point &r, 
                         double &ab_angle, double &ac_angle, double &bc_angle);
  bool XgGiveNextElement(Element &element);
  void XgForgetGrid();
  void XgSetPointsRandom();
  void XgSetPointsOverlay();
  void XgInputPoints(FILE *f, int *error, int *mem);
  void XgOutputPoints(FILE *f);
  void XgOutput(FILE *f);
  void XgRemoveBoundaryPoints();
  bool XgCreateNextTriangle(int &A, int &B, int &C, double &AB_angle, 
                            double &AC_angle, double &BC_angle, bool &finished);
  void XgNextShift(Point *p_old, Point *p_new);
  bool XgMouseAdd(Point p);
  bool XgMouseRemove(Point p_from, Point *p_where);
  bool XgRemoveInteriorPoint(Point p);
  bool XgIsEmpty();
  void XgTimeInc();
  void XgTimeDec();
  void XgSetTimestepConstant(double constant);
  Point XgGivePoint(long pos_in_list);
  Element XgGiveElement(long pos_in_list);
  bool CreateMeshBatchMode();
  Point* XgGivePoints() {return Points;}
  ~Xgen();
 

  protected:
  void XgInit(char *cfg_filename);
  void XgSetTimestep(double init_timestep);
  void XgCreateNewBoundaryComponent();
  void XgAddBoundarySegment(int marker, double x, double y, int subdiv, double alpha);
  virtual void XgUserOutput(FILE *f);
  virtual double XgCriterion(Point a, Point b, Point c);
  virtual void XgReadData(FILE *f) = 0;
  bool nogui;        // if true, operates in batch mode.
  int steps_to_take; // used with -nogui only: number of relaxation steps to be done.
  bool overlay;      // if true, initial equidistant overlay point 
                     // distribution will be used, otherwise random 
  Point* Points; 

  private:
  int BdyComponentsNum;
  bool Removal_of_boundary_points_needed;
  double H, DeltaT, TimestepConst; int Nstore;
  char *Name; int Dimension; int InteriorPtsNum, 
  BoundaryPtsNum, Npoin, Nelem; BoundaryType Boundary;
  BoundaryEdgeList* BEL; 
  BoundaryPairsList* BPL; 
  BoundaryLinesList* BLL; 
  PointList* PL; 
  ElemList* EL; 
  double Xmax, Xmin, Ymax, 
  Ymin, Area, First_DeltaT; int GoThroughPointsPtr, 
  RedrawInteriorPtr, IterationPtr, First_Npoin, First_Nstore;
  void Shift(int i, Point Impuls);
  bool FindAdmissibleThirdVertex(int A, int B, int &wanted);
  Point GetImpuls(int i);
  Point GiveImpuls(int i, int j);
  bool IsInside(Point &P);
  bool BoundaryIntersectionCheck(int a, int b);
  bool EdgesIntersect(Point a, Point b, Point c, Point d);
  bool IsLeft(Point A, Point B, Point C);  // checks whether C is on the left of AB
  bool IsRight(Point A, Point B, Point C); // checks whether C is on the right of AB
};

// Calling Motif:

void XgMainLoop(Xgen *User, int argc, char *argv[]);























