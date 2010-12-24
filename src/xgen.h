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

struct BoundaryInfo {
  int A, B;
  int ComponentIndex;
  int EdgeIndex;
  BoundaryInfo *Next;

  BoundaryInfo() {}
  BoundaryInfo(
   int a, 
   int b, 
   int CompInd, 
   int EdInd
  ) {
    A = a; B = b; 
    ComponentIndex = CompInd;
    EdgeIndex = EdInd;
    Next = NULL;
  }
};

// Private types:

struct Box {
  int A, B;
  Box *Next;
  Box(
   int a, int b
  ) {A = a; B = b; Next = NULL;}
};

struct LineList {
  Box
   *First,
   *Last,
   *Ptr;
  LineList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  int Get_ptr(int *a, int *b);
  int Is_there(int C);
  void Add(int a, int b);
  int Delete(int a, int b);
  void Delete_last();
  void Delete_list_of_lines();
  void Get_last(
   int *a, int *b
  ) {*a = Last->A; *b = Last->B;}
  void Change_last(int a, int b);
  int Is_empty() {return (First == NULL) ? 1 : 0;}
};

struct BoundaryInfoList {
  BoundaryInfo
   *First,
   *Last,
   *Ptr;
  BoundaryInfoList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  int Get_ptr(BoundaryInfo *I);
  void Add(int a, int b, 
   int CompInd, int EdInd
  );
  void Delete_last();
  void Remove();
};

struct ElemBox {
  Element E;
  ElemBox *Next;
  ElemBox(int a, 
   int b, int c
  ) {
    E.n1 = a; E.n2 = b; E.n3 = c; 
    Next = NULL;
  }
};

struct ElemList {
  ElemBox
   *First,
   *Ptr,
   *Last;
  ElemList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  int Get(int *a, 
   int *b, int *c
  );
  void Add(int a, 
   int b, int c
  );
  void Remove();
  ~ElemList();
};

struct PointBox {
  Point P;
  PointBox *Next;
  PointBox(double pos_x, double pos_y) : P(pos_x, pos_y) {
    Next = NULL;
  }
};

struct PointList {
  PointBox 
   *First,
   *Last,
   *Ptr;
  PointList() {First = Last = Ptr = NULL;}
  void Init_ptr() {Ptr = First;}
  int Get_ptr(Point *T);
  int Get(int i, double *pos_x, double *pos_y);
  void Add(Point P);
  void Add(double pos_x, double pos_y);
  Point Delete();
};  

struct Edge {
  Point P;
  int 
   Kind, Number_of_lines;
  Edge *Next;

  Edge(int k, double x, double y, int n) {
    P.x = x;
    P.y = y;
    Kind = k;
    Number_of_lines = n;
    Next = NULL;
  }
};

struct Component {
  Edge 
   *First_edge,
   *Last_edge;  
  Component *Next;

  Component() {
    First_edge = Last_edge = NULL;
    Next = NULL;
  }
};

struct Information {
  Component 
   *First_component, 
   *Last_component;
  Information() {
    First_component = Last_component = NULL; 
  }  
  void New_component();
  void Add_edge(int Kind, double x, double y, int Number_of_lines);
  LineList *Create_list_of_lines();
  PointList *Create_list_of_points();
  BoundaryInfoList *CreateBoundaryInfo();
  double GiveBoundaryLength();
  ~Information();
};

// Additional i/o functions:

int Get(FILE *f, Point *what);
int Get(FILE *f, Element *what);
int Get(FILE *f, BoundaryInfo *what);
void Put(FILE *f, Point what);
void Put(FILE *f, Element what);
void Put(FILE *f, BoundaryInfo what);
void Put(Point what);
void Put(Element what);
void Put(BoundaryInfo what);
void PutNl(FILE *f, Point what);
void PutNl(FILE *f, Element what);
void PutNl(FILE *f, BoundaryInfo what);
void PutNl(Point what);
void PutNl(Element what);
void PutNl(BoundaryInfo what);

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
  void XgGiveLimits(Point *min, Point *max);
  void XgInitPointList();
  int XgGiveNextPoint(Point *p);
  void XgInitFreePointList();
  int XgGiveNextFreePoint(Point *p);
  void XgInitBoundaryLineList();
  int XgGiveNextBoundaryLine(Point *p, Point *q);
  void XgInitBoundaryInfoList();
  int XgGiveNextBoundaryInfo(BoundaryInfo *info);
  void XgInitElementList();
  int XgGiveNextElement(Point *p, Point *q, Point *r);
  int XgGiveNextElement(Element *element);
  void XgForgetGrid();
  void XgSetPointsRandom();
  void XgSetPointsOverlay();
  void XgInputPoints(FILE *f, int *error, int *mem);
  void XgOutputPoints(FILE *f);
  void XgOutput(FILE *f);
  int XgNextTriangle(Point *p, Point *q, Point *r, int *ready);
  void XgNextShift(Point *p_old, Point *p_new);
  int XgMouseAdd(Point p);
  int XgMouseRemove(Point p_from, Point *p_where);
  int XgRemoveFreePoint(Point p);
  int XgIsEmpty();
  void XgTimeInc();
  void XgTimeDec();
  void XgSetTimestepConstant(double constant);
  Point XgGivePoint(long pos_in_list);
  Element XgGiveElement(long pos_in_list);
  ~Xgen();

  protected:
  void XgInit(char *cfg_filename);
  void XgSetTimestep(double init_timestep);
  void XgNewComponent();
  void XgAddEdge(int index, double x, double y);
  void XgAddEdge(int index, double x, double y, int lines);
  void XgAddEdge(int index, Point P);
  void XgAddEdge(int index, Point P, int lines);
  virtual void XgUserOutput(FILE *f);
  virtual double XgCriterion(Point a, Point b, Point c);
  virtual void XgReadData(FILE *f) = 0;
  bool nogui;  // if trye, operates in batch mode.
  int steps_to_take; // used with -nogui only: number of relaxation steps to be done.
  bool overlay; // if true, initial equidistant overlay point 
                // distribution will be used, otherwise random 

  private:
  double H, DeltaT, TimestepConst; int Nstore;
  char *Name; int Dimension; int Nfree, 
  Nbound, Npoin, Nelem; Point *E; Information Info; 
  LineList *L; PointList *EL; ElemList *ElemL; 
  BoundaryInfoList *BIL; double Xmax, Xmin, Ymax, 
  Ymin, Area, First_DeltaT; int GoThroughPointsPtr, 
  RedrawFreePtr, IterationPtr, First_Npoin, First_Nstore;
  void Shift(int i, Point Impuls);
  int Find_nearest_left(int A, int B, 
  int *wanted);
  Point Get_impuls(int i);
  Point Give_impuls(int i, int j);
  int Is_inside(Point *P);
  int Boundary_intact(Point a, Point b, Point c);
  int Intact(Point a, Point b, Point c, Point d);
  int Is_left(Point A, Point B, Point C);
  int Is_right(Point A, Point B, Point C);
};

// Calling Motif:

void XgMainLoop(Xgen *User, int argc, char *argv[]);























