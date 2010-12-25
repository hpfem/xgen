# include <stdlib.h>
# include <stdio.h>
# include <math.h>

// Public types:

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
  int n1, n2, n3;
  Element() {}
  Element(int in1, int in2, int in3) {
    n1 = in1; n2 = in2; n3 = in3;
  } 
};

struct BoundaryEdge {
  int A, B;            // Vertex indices.
  double alpha;        // 0 for straight edges, nonzero for circular arcs.
  int component_index; // Index of boundary component (closed loop).
  int marker;          // Boundary marker.
  BoundaryEdge *next;

  BoundaryEdge() {}
  BoundaryEdge(int a, int b, int ci, int m, double alp) 
  {
    A = a; B = b; 
    component_index = ci;
    marker = m;
    alpha = alp;
    next = NULL;
  }
};

// Private types:

struct Box {
  int A, B;
  Box *next;
  Box(
   int a, int b
  ) {A = a; B = b; next = NULL;}
};

struct LineList {
  Box *First, *Last, *Ptr;
  LineList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  bool Get_ptr(int &a, int &b);
  bool Is_there(int C);
  void Add(int a, int b);
  bool Delete(int a, int b);
  void Delete_last();
  void Delete_list_of_lines();
  void Get_last(int &a, int &b) {a = Last->A; b = Last->B;}
  void Change_last(int a, int b);
  bool Is_empty() {return (First == NULL) ? 1 : 0;}
};

struct BoundaryEdgeList {
  BoundaryEdge *First, *Last, *Ptr;
  BoundaryEdgeList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  bool Get_ptr(BoundaryEdge &I);
  void Add(int a, int b, int component_index, int marker, double alpha);
  void Delete_last();
  void Remove();
};

struct ElemBox {
  Element E;
  ElemBox *next;
  ElemBox(int a, int b, int c) {
    E.n1 = a; E.n2 = b; E.n3 = c; 
    next = NULL;
  }
};

struct ElemList {
  ElemBox *First, *Ptr, *Last;
  ElemList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  bool Get(int &a, int &b, int &c);
  void Add(int a, int b, int c);
  void Remove();
  ~ElemList();
};

struct PointBox {
  Point P;
  PointBox *next;
  PointBox(double pos_x, double pos_y) : P(pos_x, pos_y) {
    next = NULL;
  }
};

struct PointList {
  PointBox *First, *Last, *Ptr;
  PointList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  bool Get_ptr(Point &T);
  bool Get(int i, double &pos_x, double &pos_y);
  void Add(Point P);
  void Add(double pos_x, double pos_y);
  Point Delete();
};  

struct Edge {
  Point P;
  int marker, subdiv;
  double alpha;
  Edge *next;

  Edge(int m, double x, double y, int n, double alp) {
    P.x = x;
    P.y = y;
    marker = m;
    subdiv = n;
    alpha = alp;
    next = NULL;
  }
};

struct BdyComponent {
  Edge *First_edge, *Last_edge;  
  BdyComponent *next;

  BdyComponent() {
    First_edge = Last_edge = NULL;
    next = NULL;
  }
};

struct Boundary {
  BdyComponent *First_bdy_component, *Last_bdy_component;
  Boundary() {First_bdy_component = Last_bdy_component = NULL;}  
  void Add_bdy_component();
  void Add_bdy_segment(int marker, double x, double y, int subdiv, double alpha);
  LineList *Create_list_of_lines();
  PointList *Create_list_of_points();
  BoundaryEdgeList *CreateBoundaryEdges();
  double GiveBoundaryLength();
  ~Boundary();
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

// class Xgen:

class Xgen {
  public:
  Xgen(bool nogui, int steps_to_take, bool overlay);
  char* XgGiveName();
  double XgGiveTimestep();
  int XgGiveDimension();
  long XgGiveNpoin();
  long XgGiveNfree();
  long XgGiveNbound();
  long XgGiveNelem();
  void XgGiveLimits(Point &min, Point &max);
  void XgInitPointList();
  bool XgGiveNextPoint(Point &p);
  void XgInitFreePointList();
  bool XgGiveNextFreePoint(Point &p);
  void XgInitBoundaryLineList();
  bool XgGiveNextBoundaryLine(Point &p, Point &q);
  void XgInitBoundaryEdgeList();
  bool XgGiveNextBoundaryEdge(BoundaryEdge &edge_ptr);
  void XgInitElementList();
  bool XgGiveNextElement(Point &p, Point &q, Point &r);
  bool XgGiveNextElement(Element &element);
  void XgForgetGrid();
  void XgSetPointsRandom();
  void XgSetPointsOverlay();
  void XgInputPoints(FILE *f, int *error, int *mem);
  void XgOutputPoints(FILE *f);
  void XgOutput(FILE *f);
  bool XgNextTriangle(Point &p, Point &q, Point &r, bool &ready);
  void XgNextShift(Point *p_old, Point *p_new);
  bool XgMouseAdd(Point p);
  bool XgMouseRemove(Point p_from, Point *p_where);
  bool XgRemoveFreePoint(Point p);
  bool XgIsEmpty();
  void XgTimeInc();
  void XgTimeDec();
  void XgSetTimestepConstant(double constant);
  Point XgGivePoint(long pos_in_list);
  Element XgGiveElement(long pos_in_list);
  bool CreateMeshBatchMode();
  ~Xgen();

  protected:
  void XgInit(char *cfg_filename);
  void XgSetTimestep(double init_timestep);
  void XgAddBdyComponent();
  void XgAddBdySegment(int marker, double x, double y, int subdiv, double alpha);
  virtual void XgUserOutput(FILE *f);
  virtual double XgCriterion(Point a, Point b, Point c);
  virtual void XgReadData(FILE *f) = 0;
  bool nogui;        // if true, operates in batch mode.
  int steps_to_take; // used with -nogui only: number of relaxation steps to be done.
  bool overlay;      // if true, initial equidistant overlay point 
                     // distribution will be used, otherwise random 

  private:
  double H, DeltaT, TimestepConst; int Nstore;
  char *Name; int Dimension; int Nfree, 
  Nbound, Npoin, Nelem; Point *E; Boundary Info; 
  LineList *L; PointList *EL; ElemList *ElemL; 
  BoundaryEdgeList *BIL; double Xmax, Xmin, Ymax, 
  Ymin, Area, First_DeltaT; int GoThroughPointsPtr, 
  RedrawFreePtr, IterationPtr, First_Npoin, First_Nstore;
  void Shift(int i, Point Impuls);
  bool Find_nearest_left(int A, int B, int &wanted);
  Point Get_impuls(int i);
  Point Give_impuls(int i, int j);
  bool Is_inside(Point &P);
  bool Boundary_intact(Point a, Point b, Point c);
  bool Intact(Point a, Point b, Point c, Point d);
  bool Is_left(Point A, Point B, Point C);
  bool Is_right(Point A, Point B, Point C);
};

// Calling Motif:

void XgMainLoop(Xgen *User, int argc, char *argv[]);























