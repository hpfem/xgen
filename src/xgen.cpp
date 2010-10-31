/**
 **  2D interactive grid generator
 **            class  xgen 
 **             1.3.1998
 **           (Pavel Solin)
 **
 **        All rights reserved.
 **      Not for commercial use!
 **/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <strings.h>
# include <unistd.h>
# include <time.h>

# include "disc.h"
# include "xgen.h"

double XG_ZERO = 1e-12;

/*
 *   Additional i/o functions:
 */

int REMOVE_BDY_PTS_ACTIVE = 1;

int Get(FILE *f, Point *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  what->x = (double)atof(str);
  if(!Get(f, str)) return 0;
  what->y = (double)atof(str);
  return 1;
}

int Get(FILE *f, Element *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  what->n1 = (int)atoi(str);
  if(!Get(f, str)) return 0;
  what->n2 = (int)atoi(str);
  if(!Get(f, str)) return 0;
  what->n3 = (int)atoi(str);
  return 1;
}

int Get(FILE *f, BoundaryInfo *what) {
  char str[255];
  if(!Get(f, str)) return 0;
  what->A = (int)atoi(str);
  if(!Get(f, str)) return 0;
  what->B = (int)atoi(str);
  if(!Get(f, str)) return 0;
  what->EdgeIndex = (int)atol(str);
  return 1;
}

void Put(FILE *f, Point what) {
  fprintf(f, "%g %g\n", what.x, what.y);
}

void Put(FILE *f, Element what) {
  fprintf(f, "%d %d %d\n", what.n1 + 1, what.n2 + 1, what.n3 + 1);
}

void Put(FILE *f, BoundaryInfo what) {
  fprintf(f, "%d %d %d\n", what.A + 1, what.B + 1, what.EdgeIndex);
}

void Put(Point what) {
  fprintf(stderr, "%g %g\n", what.x, what.y);
}

void Put(Element what) {
  fprintf(stderr, "%d %d %d\n", what.n1 + 1, what.n2 + 1, what.n3 + 1);
}

void Put(BoundaryInfo what) {
  fprintf(stderr, "%d %d %d\n", what.A + 1, what.B + 1, what.EdgeIndex);
}

/*
 *   Errors & Warnings:
 */

void XgError(char *who, char *what) {
  fprintf(stderr, "%s Error: %s\n", who, what);
  exit(0);
}

void XgError(const char *who, const char *what) {
  fprintf(stderr, "%s Error: %s\n", who, what);
  exit(0);
}

void XgError(char *what) {
  fprintf(stderr, "Error: %s\n", what);
  exit(0);
}

void XgError(const char *what) {
  fprintf(stderr, "Error: %s\n", what);
  exit(0);
}

void XgWarning(char *who, char *what) {
  fprintf(stderr, "%s Warning: %s\n", who, what);
}

void XgWarning(char *what) {
  fprintf(stderr, "Warning: %s\n", what);
}

void XgMessage(char *what) {
  fprintf(stderr, "%s\n", what);
}

/**
 **   struct LineList:
 **/

void LineList::Change_last(
 int a, int b
) {
  {Last->A = a; Last->B = b;}
}

int LineList::Get_ptr(
 int *a, int *b
) {
  if(Ptr != NULL) {
    *a = Ptr->A;
    *b = Ptr->B;
    Ptr = Ptr->Next;
    return 1;
  }
  else return 0;
}

void LineList::Delete_list_of_lines() {
  Ptr = First;
  Box *p;
  while(Ptr != NULL) {
    p = Ptr;
    Ptr = Ptr -> Next;
    delete p;
  }
  First = Last = NULL;
}

int LineList::Is_there(int C) {
  Box *P = First;
  while(P != NULL) {
    if((P->A == C) || (P->B == C)) return 1;
    P = P->Next;
  }
  return 0;
}

void LineList::Add(
 int a, int b
) {
  if (Last == NULL) {
    First = Last = new Box(a, b);
  }
  else {
    Last->Next = new Box(a, b);
    Last = Last->Next;
  }
}

void LineList::Delete_last() {
  Box *p = First;
  if (p->Next == NULL) {
    delete p;
    First = Last = NULL;
  }
  else {
    if (p->Next->Next == NULL) {
      delete Last;
      Last = First;
    }
    else {
      while (p->Next->Next != NULL) p = p->Next;
      Last = p;
      delete p->Next;
    }
    Last->Next = NULL;
  }
}

int LineList::Delete(
 int a, int b
) {
  Box *p = First, *l = NULL;
  while (p != NULL) {
    if ((p->A == a && p->B == b) || (p->A == b && p->B == a)) {
      if (l == NULL) {
        First = p->Next;
	if (p->Next == NULL) Last = NULL;
      }
      else {
	l->Next = p->Next;
	if (p->Next == NULL) {
	  Last = l;
	  Last->Next = NULL;
        }
      }
      delete p;
      return 1;
    }
    else {
      l = p;
      p = p->Next;
    }
  }
  return 0;
}

/**
 **   struct BoundaryInfoList:
 **/

int BoundaryInfoList::Get_ptr(BoundaryInfo *I) {
  if(Ptr != NULL) {
    *I = *Ptr;
    Ptr = Ptr->Next;
    return 1;
  }
  else return 0;
}

void BoundaryInfoList::Add(int a, 
 int b, int CompInd, int EdInd
) {
  if (Last == NULL) {
    First = Last = new BoundaryInfo(a, b, CompInd, EdInd);
  }
  else {
    Last->Next = new BoundaryInfo(a, b, CompInd, EdInd);
    Last = Last->Next;
  }
}

void BoundaryInfoList::Delete_last() {
  BoundaryInfo *p = First;
  if (p->Next == NULL) {
    delete p;
    First = Last = NULL;
  }
  else {
    if (p->Next->Next == NULL) {
      delete Last;
      Last = First;
    }
    else {
      while (p->Next->Next != NULL) p = p->Next;
      Last = p;
      delete p->Next;
    }
    Last->Next = NULL;
  }
}

void BoundaryInfoList::Remove() {
  BoundaryInfo *p = Ptr = First;
  while(p != NULL) {
    Ptr = Ptr->Next;
    delete p;
    p = Ptr;
  }
  Last = First = NULL;
}

/**
 **    struct ElemList:
 **/

void ElemList::Add(int a, 
 int b, int c
) {
  if (Last == NULL) {
    First = Last = Ptr = new ElemBox(a, b, c);
  }
  else {
    Last->Next = new ElemBox(a, b, c);
    Last = Last->Next;
    if(Last == NULL) {
      XgError(
       "ElemList",
       "Not enough memory for triangles."
      ); 
    }
  }
}

int ElemList::Get(int *a, 
 int *b, int *c
) {
  if(Ptr != NULL) {
    *a = Ptr->E.n1;
    *b = Ptr->E.n2;
    *c = Ptr->E.n3;
    Ptr = Ptr->Next;
    return 1;
  }
  else return 0;
}

void ElemList::Remove() {
  ElemBox *p;
  Ptr = First;
  while(Ptr != NULL) { 
    p = Ptr;
    Ptr = Ptr->Next;
    delete p;
  }
  First = Last = NULL; 
}

/**
 **   struct PointList:
 **/

int PointList::Get(int i, 
 double *pos_x, double *pos_y
) {
  PointBox *P = First;
  for(int j=0; j<i; j++) {
    if(P == NULL) return 0;
    P = P->Next;
  }
  *pos_x = P->P.x;
  *pos_y = P->P.y;
  return 1;
}

int PointList::Get_ptr(Point *T) {
  if(Ptr != NULL) {
    T->x = Ptr->P.x;
    T->y = Ptr->P.y;
    Ptr = Ptr->Next; 
    return 1;
  }
  else return 0;
}

void PointList::Add(Point P) {
  Add(P.x, P.y);
}

void PointList::Add(double pos_x, double pos_y) {
//  Put("PointList: adding point "); 
//  Put(pos_x); Put(" "); 
//  PutNl(pos_y);
  if (Last == NULL) {
    First = Last = new PointBox(pos_x, pos_y);
  }
  else {
    Last->Next = new PointBox(pos_x, pos_y);
    Last = Last->Next;
  }
}

Point PointList::Delete() {
  Point P(0, 0);

  if(First == NULL) return P;
  else {
    PointBox *p = First;
    P = p->P;
    First = First->Next;
    delete p;
    return P;
  }
}

/**
 **   struct Information:
 **/

void Information::New_component() {
//  PutNl("Info: adding new component.");
  if(Last_component == NULL) {
    First_component = Last_component = new Component();
    if(First_component == NULL) {
      XgError(
       "struct Information",
       "Not enough memory for new component." 
      );
    }
  }
  else {
    Component *help = Last_component;
    Last_component = new Component();
    help -> Next = Last_component;
    if(help -> Next == NULL) {
      XgError(
       "struct Information",
       "Not enough memory for new component." 
      );
    }
  }
}

void Information::Add_edge(
  int Kind, double x, double y, int Number_of_lines
) {
  if(Last_component == NULL) {
    First_component = Last_component = new Component();
    if(First_component == NULL) {
      XgError(
       "struct Information",
       "Not enough memory for new edge." 
      );
    }    
  }
  if(Last_component -> Last_edge == NULL) {
    Last_component -> Last_edge = 
    Last_component -> First_edge = 
    new Edge(Kind, x, y, Number_of_lines);
    if(Last_component -> Last_edge == NULL) {
      XgError(
       "struct Information",
       "Not enough memory for new edge." 
      );
    } 
//    Put("Info: adding edge "); Put(Kind); Put(" ");
//    PutNl(Number_of_lines);
  } 
  else {
    Edge *help = Last_component -> Last_edge;
    Last_component -> Last_edge = 
    new Edge(Kind, x, y, Number_of_lines);
    help -> Next = Last_component -> Last_edge;
    if(help -> Next == NULL) {
      XgError(
       "struct Information",
       "Not enough memory for new edge." 
      );
    }
//    Put("Info: adding edge "); Put(Kind); Put(" ");
//    PutNl(Number_of_lines);
  }
}

LineList *Information::Create_list_of_lines() {
  LineList *l = new LineList;
  int Count = 0, First_point_of_component = 0;
  int Number_of_lines = 0;

//  PutNl("\nCreating LineList:");
  Component *Component_pointer = First_component;
  while(Component_pointer != NULL) {
    Edge *Edge_pointer = Component_pointer -> First_edge;  
    First_point_of_component = Count;
    while(Edge_pointer != NULL) {
      Number_of_lines = Edge_pointer -> Number_of_lines;
//      Put("Number_of_lines: "); 
//      PutNl(Number_of_lines);
      for(int i=0; i<Number_of_lines; i++) {
//        Put("LineList: adding "); Put(Count); Put(" "); 
//        PutNl(Count+1);
        l->Add(Count, Count + 1);
        Count++;
      }
      Edge_pointer = Edge_pointer -> Next;
    }
    l->Delete_last();
//    PutNl("LineList: deleting last item.");
    l->Add(Count-1, First_point_of_component);
//    Put("LineList: adding "); Put(Count-1); Put(" ");
//    PutNl(First_point_of_component);
    Component_pointer = Component_pointer -> Next;
  }
  return l;
}       

PointList *Information::Create_list_of_points() {
  PointList *l = new PointList;
  Point 
   First_point_of_component(0, 0),
   First_point(0, 0), 
   Next_point(0, 0), 
   P(0, 0);
  int Number_of_lines;

//  Put("\nCreating PointList:\n");
  Component *Component_pointer = First_component;
  while(Component_pointer != NULL) {
    Edge *Edge_pointer = Component_pointer -> First_edge;  
    First_point_of_component = First_point = Edge_pointer->P;
    while(Edge_pointer->Next != NULL) {
      Number_of_lines = Edge_pointer -> Number_of_lines;
//      Put("Number_of_lines: "); PutNl(Number_of_lines);
      Edge_pointer = Edge_pointer->Next;
      Next_point = Edge_pointer->P;
      P = (Next_point - First_point)/Number_of_lines;
      for(int i=0; i<Number_of_lines; i++) {
        l->Add(First_point + P*i);
      }
      First_point = Next_point;
    }
    Number_of_lines = Edge_pointer -> Number_of_lines;
//    Put("Number_of_lines: "); PutNl(Number_of_lines);
    P = (First_point_of_component - First_point)/Number_of_lines;
    for(int i=0; i<Number_of_lines; i++) l->Add(First_point + P*i);
    Component_pointer = Component_pointer -> Next;
  }
  return l;
}       

double Information::GiveBoundaryLength() {
  double length = 0;
  Point 
   First_point_of_component(0, 0),
   First_point(0, 0), 
   Next_point(0, 0);

//  Put("\nComputing boundary length:\n");
  Component *Component_pointer = First_component;
  while(Component_pointer != NULL) {
    Edge *Edge_pointer = Component_pointer -> First_edge;  
    First_point_of_component = First_point = Edge_pointer->P;
    while(Edge_pointer->Next != NULL) {
      Edge_pointer = Edge_pointer->Next;
      Next_point = Edge_pointer->P;
      length += (Next_point - First_point).abs();
      First_point = Next_point;
    }
    length += (First_point_of_component - First_point).abs();
    Component_pointer = Component_pointer -> Next;
  }
  return length;
}       

BoundaryInfoList *Information::CreateBoundaryInfo() {
  BoundaryInfoList *l = new BoundaryInfoList;
  int Count = 0,
   First_point_of_component = 0;
  int Number_of_lines = 0, ComponentIndex = 0;
  int 
   Last_kind = 0;

//  Put("\nCreating BoundaryInfoList:\n");
  Component *Component_pointer = First_component;
  while(Component_pointer != NULL) {
    Edge *Edge_pointer = Component_pointer -> First_edge;
    First_point_of_component = Count;
    while(Edge_pointer != NULL) {
      Number_of_lines = Edge_pointer -> Number_of_lines;
      for(int i=0; i<Number_of_lines; i++) {
        l->Add(
            Count, Count + 1, 
            ComponentIndex,
            Edge_pointer -> Kind
           );
//        Put("BoundaryInfoList: adding "); 
//        Put(Count); Put(" "); 
//        Put(Count+1); Put(" ");
//        Put(ComponentIndex); Put(" ");
//        Put(Edge_pointer -> Kind);
//        Put("\n");
        Count++;
      }
      Last_kind = Edge_pointer -> Kind;
      Edge_pointer = Edge_pointer -> Next;
    }
    l->Delete_last();
//    Put("deleting last item, ");    
    l->Add(
        Count - 1, First_point_of_component,
        ComponentIndex,
        Last_kind
       );
//    Put("BoundaryInfoList: adding "); 
//    Put(Count-1); Put(" ");
//    Put(First_point_of_component); Put(" ");
//    Put(ComponentIndex); Put(" ");
//    Put(Last_kind);
//    Put("\n");
    Component_pointer = Component_pointer -> Next;;
    ComponentIndex++;
  }
  return l;
}       

Information::~Information() {
  Component *comp_ptr = First_component;
  while (comp_ptr != NULL) {
    Edge *edge_ptr = comp_ptr -> First_edge;
    while(edge_ptr != NULL) {
      Edge *help_edge = edge_ptr;
      edge_ptr = edge_ptr -> Next;
      delete help_edge;
    }
    Component *help_comp = comp_ptr;
    comp_ptr = comp_ptr -> Next;
    delete help_comp;
  }
  First_component = Last_component = NULL;
}

/**
 **   class Xgen: public methods
 **/

char* Xgen::XgGiveName() {
  char *c = (char*)malloc(255);
  if(c == NULL) XgError("No free memory to allocate.");
  c[0] = '\0';
  strcpy(c, Name);
  return c;
}

double Xgen::XgGiveTimestep() {
  return DeltaT;
}

int Xgen::XgGiveDimension() {
  return Dimension;
}

long Xgen::XgGiveNpoin() {
  return Npoin;
}

long Xgen::XgGiveNfree() {
  return Nfree;
}

long Xgen::XgGiveNbound() {
  return Nbound;
}

long Xgen::XgGiveNelem() {
  return Nelem;
}

void Xgen::XgGiveLimits(Point *min, Point *max) {
  min->x = Xmin;
  min->y = Ymin;
  max->x = Xmax;
  max->y = Ymax;
}

void Xgen::XgInitPointList() {
  GoThroughPointsPtr = 0;
}

int Xgen::XgGiveNextPoint(Point *p) {
  if(GoThroughPointsPtr < Npoin) {
    *p = E[GoThroughPointsPtr++];
    return 1;
  }
  else return 0;
} 

void Xgen::XgInitFreePointList() {
  RedrawFreePtr = Nbound;
}

int Xgen::XgGiveNextFreePoint(Point *p) {
  if(RedrawFreePtr < Npoin) {
    *p = E[RedrawFreePtr++];
    return 1;
  }
  else return 0;
}

void Xgen::XgInitBoundaryLineList() {
  L->Init_ptr();
}

int Xgen::XgGiveNextBoundaryLine(Point *p, Point *q) {
  int A, B;
  if(L->Get_ptr(&A, &B)) {
    *p = E[A];
    *q = E[B];
    return 1;
  }
  else return 0;
}

void Xgen::XgInitBoundaryInfoList() {
  BIL->Init_ptr();
}

int Xgen::XgGiveNextBoundaryInfo(BoundaryInfo *info) {
  if(BIL->Get_ptr(info)) return 1;
  else return 0;
}

void Xgen::XgInitElementList() {
  ElemL->Init_ptr();
}

int Xgen::XgGiveNextElement(Point *p, Point *q, Point *r) {
  int A, B, C;
  if(ElemL->Get(&A, &B, &C)) {
    *p = E[A];
    *q = E[B];
    *r = E[C];
    return 1;
  }
  else return 0;
}

int Xgen::XgGiveNextElement(Element *element) {
  if(ElemL->Get(&(element->n1), 
   &(element->n2), &(element->n3))) return 1;
  else return 0;
}

void Xgen::XgForgetGrid() {
  ElemL->Remove();
  L->Delete_list_of_lines();
  L = Info.Create_list_of_lines();
  Nelem = 0;
  REMOVE_BDY_PTS_ACTIVE = 1;
}

//points are randomly set
void Xgen::XgStartAgain() {
  Npoin = First_Npoin;
  Nstore = First_Nstore;
  if(Nbound > Npoin) XgError("Internal (1) in XgStartAgain().");
  Nfree = Npoin - Nbound;
#ifdef RAND_MAX
  double Rand_max = RAND_MAX;
#else
  double Rand_max = pow(2.0, 15.0);
#endif
  srand(1);
  for(int i=Nbound; i<Npoin; i++) {
    do {
      E[i].x =
        rand()*(Xmax - Xmin)/Rand_max + Xmin;
      E[i].y =
	rand()*(Ymax - Ymin)/Rand_max + Ymin;
    } while(!Is_inside(E + i));
  }
  printf("Domain filled with random points, Npoin = %d, Nfree = %d.\n", Npoin, Nfree);
}

int odd(int i) {
  if(i & 1) return 1;
  else return 0;
}

//points are equidistantly placed
void Xgen::XgStartAgain2() {
  //Npoin = First_Npoin;
  //Nstore = First_Nstore;
  //Nfree = Npoin - Nbound;
  //#ifdef RAND_MAX
  //  double Rand_max = RAND_MAX;
  //#else
  //  double Rand_max = pow(2.0, 15.0);
  //#endif
  //srand(1);
  //for(int i=Nbound; i<Npoin; i++) {
  //  do {
  //    E[i].x =
  //      rand()*(Xmax - Xmin)/Rand_max + Xmin;
  //    E[i].y =
  //      rand()*(Ymax - Ymin)/Rand_max + Ymin;
  //  } while(!Is_inside(E + i));
  //}

  //I can place Nfree + Nstore new points!
  int Nx = (int)((Xmax-Xmin)/H);
  int Ny = (int)((Ymax-Ymin)/(H*sqrt(3)/2.0));
  Npoin = Nbound;

#ifdef RAND_MAX
  double Rand_max = RAND_MAX;
#else
  double Rand_max = pow(2.0, 15.0);
#endif
  srand(1);
  for(int i=0; i<Ny; i++) {
    for(int j=0; j<Nx; j++) {
      Point P;
      double fluct_x = rand()*0.1*H/Rand_max;
      double fluct_y = rand()*0.1*H/Rand_max;
      P.x = Xmin + H/2 + H*j + 0.5*H*odd(i) + fluct_x;
      P.y = Ymin + (i+1)*H*sqrt(3)/2.0 + fluct_y;
      if(Is_inside(&P)) {
        if(Npoin >= First_Npoin + First_Nstore) 
          XgError("not enough space for new points, increase the Nstore parameter."); 
        E[Npoin++] = P;
      }
    }
  }      

  if(Nbound > Npoin) XgError("Internal (1) in XgStartAgain2().");
  Nfree = Npoin - Nbound;
  Nstore = First_Nstore + First_Npoin - Npoin;
  IterationPtr = Nbound;
  RedrawFreePtr = Nbound;
  GoThroughPointsPtr = Nbound;

  printf("Domain filled with overlay points, Npoin = %d, Nfree = %d.\n", Npoin, Nfree);
}

void Xgen::XgOutput(FILE *f) {
  BIL = Info.CreateBoundaryInfo();
  XgUserOutput(f);
  BIL -> Remove();
} 

void Xgen::XgInputPoints(FILE *f, int *error, int *mem) {
  *error = 0;
  *mem = 0;
  
  Nstore = Nstore + Npoin - Nbound;
  Npoin = Nbound;
  IterationPtr = Nbound;

  Point P, *E0;
  E0 = E + Nbound;
  while(Get(f, &P)) if(Is_inside(&P)) {
    Npoin++; Nstore--;
    if(Nstore == 0) break;
    *E0 = P;
    E0++;
  }

  if(Nbound > Npoin) XgError("Internal (1) in XgInputPoints().");
  Nfree = Npoin - Nbound;
  return;
} 

void Xgen::XgOutputPoints(FILE *f) {
  fprintf(f, "# List of free points:\n");

  Point P;
  GoThroughPointsPtr = Nbound;
  while(XgGiveNextPoint(&P)) fprintf(f, "%g %g\n", P.x, P.y);
} 

void Undraw_point1(Point );

int Xgen::XgNextTriangle(Point *p, Point *q, Point *r, int *ready) {

  *ready = 0;
  int A, B, C;

  //removing points lying near the boundary
  if(REMOVE_BDY_PTS_ACTIVE) {
    double coeff = 4.0;  ///THIS COEFFICIENT HAS TO BE SET HERE!

    int was_deleted = 0;
    Point a, b;
    XgInitBoundaryLineList();                                           
    while(XgGiveNextBoundaryLine(&a, &b)) {
      double h = sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
      double dh = h/coeff;
      XgInitFreePointList();
      Point c;
      while(XgGiveNextFreePoint(&c)) {
        Point ab, ab2, ab3;
        ab.x = 0.5*(a.x + b.x);
        ab.y = 0.5*(a.y + b.y);
        ab2.x = 0.5*(ab.x + b.x);
        ab2.y = 0.5*(ab.y + b.y);
        ab3.x = 0.5*(a.x + ab.x);
        ab3.y = 0.5*(a.y + ab.y);
        double dist_a_c = 
          sqrt((a.x - c.x)*(a.x - c.x) + (a.y - c.y)*(a.y - c.y));
        double dist_b_c = 
          sqrt((b.x - c.x)*(b.x - c.x) + (b.y - c.y)*(b.y - c.y));
        double dist_ab_c = 
          sqrt((ab.x - c.x)*(ab.x - c.x) + (ab.y - c.y)*(ab.y - c.y));
        double dist_ab2_c = 
          sqrt((ab2.x - c.x)*(ab2.x - c.x) + (ab2.y - c.y)*(ab2.y - c.y));
        double dist_ab3_c = 
          sqrt((ab3.x - c.x)*(ab3.x - c.x) + (ab3.y - c.y)*(ab3.y - c.y));

        if(dist_a_c < dh || dist_b_c < dh || dist_ab_c < dh ||
           dist_ab2_c < dh || dist_ab3_c < dh) {
          printf("Point [%g, %g] closer than h/%g to the boundary -> deleting.\n",
            c.x, c.y, coeff);
          was_deleted++;
          if(!XgRemoveFreePoint(c)) 
            XgError("internal while removing boundary close points.");
          //break;
        }
      }                
    }                 
    REMOVE_BDY_PTS_ACTIVE = 0;
    if(was_deleted > 0) printf("Deleted %d points, Npoin = %d, Nfree = %d.\n",
      was_deleted, Npoin, Nfree);
  }

  L->Get_last(&A, &B);

  if(!Find_nearest_left(A, B, &C)) return 0;

  if (L->Delete(C, A)) {
    if (L->Delete(B, C)) L->Delete_last();
    else L->Change_last(C, B);
  }
  else {
  L->Change_last(A, C);
    if (!L->Delete(B, C)) L->Add(C, B);
  }

  if(L->Is_empty()) *ready = 1; 

  // they are well oriented!
  ElemL -> Add(A, B, C);

  Nelem++;

  *p = E[A];
  *q = E[B];
  *r = E[C];
  
  return 1;
}

void Xgen::XgNextShift(Point *p_old, Point *p_new) {
  if(Npoin > Nbound) {
    if(IterationPtr >= Npoin) {
      IterationPtr = Nbound;
//      { static int count = 0; 
//        PutNl(count++);
//      }
    }
    Point impuls = Get_impuls(IterationPtr);
    *p_old = E[IterationPtr];
    Shift(IterationPtr, impuls);
    *p_new = E[IterationPtr++];
  }
}

int Xgen::XgMouseAdd(Point p) {
  if(Nstore == 0) return 0;
  if(!Is_inside(&p)) return 0;
  E[Npoin] = p;
  Nfree++;
  Npoin++;
  Nstore--;
  IterationPtr=Nbound;
  printf("Point at [%g, %g] added to the end of the list, new Npoin = %d, Nfree = %d.\n", 
    p.x, p.y, Npoin, Nfree);
  return 1;
}

int Xgen::XgMouseRemove(Point p_from, Point *p_where) {
  if(Nfree == 0) return 0;
  int wanted = Nbound;  //case it is the last free one
  double min = 1e100;
  for(int j=Nbound; j < Npoin; j++) {
    Point Dr = E[j] - p_from;
    if(Dr.abs() < min) {
      min = Dr.abs();
      wanted = j;
    }
  }
  *p_where = E[wanted];
  for(int j = wanted + 1; j < Npoin; j++) E[j - 1] = E[j];
  Npoin--;
  Nfree--;
  Nstore++;
  IterationPtr = Nbound;
  printf("Point no. %d at [%g, %g] removed by mouse, Npoin = %d, Nfree = %d.\n", 
    wanted, E[wanted].x, E[wanted].y, Npoin, Nfree);
  return 1;
}

int Xgen::XgRemoveFreePoint(Point p) {
  if(Nfree <= 0) return 0;
  int wanted = -1;  //number of the point to be found
  double min = 1e100;
  //going through free points
  for(int j=Nbound; j < Npoin; j++) {
    double dist = sqrt((E[j].x - p.x)*(E[j].x - p.x) + 
                       (E[j].y - p.y)*(E[j].y - p.y));
    if(dist < min) {
      min = dist;
      wanted = j;
    }
  }
  if(wanted == -1) XgError("internal in XgRemoveFreePoint.");
  for(int j = wanted; j< Npoin-1; j++) E[j] = E[j+1];
  Npoin--;
  Nfree--;
  Nstore++;
  RedrawFreePtr--;

  Undraw_point1(p);
  IterationPtr = Nbound;

  return 1;
}

int Xgen::XgIsEmpty() {
  return (Npoin == Nbound);
}

void Xgen::XgTimeInc() {
    DeltaT /= TimestepConst;
}

void Xgen::XgTimeDec() {
    DeltaT *= TimestepConst;
}

void Xgen::XgSetTimestepConstant(double constant) {
  TimestepConst = constant;
  if(TimestepConst <= 0) XgError("Invalid TimestepConst.");
  if(TimestepConst >= 1) XgError("Invalid TimestepConst.");
}

Point Xgen::XgGivePoint(long pos_in_list) {
  if(pos_in_list < 0 || pos_in_list >= Npoin)
    XgError("Invalid point index used.");
  return E[pos_in_list];
}

Element Xgen::XgGiveElement(long pos_in_list) {
  if(pos_in_list < 0 || pos_in_list >= Nelem)
    XgError("Invalid element index used.");
  ElemBox *ptr = ElemL->First;
  for(long i=0; i<pos_in_list; i++) ptr = ptr->Next;
  return ptr->E;
}

/**
 **   class Xgen: protected methods
 **/

void Xgen::XgInit(char *cfg_filename) {
  if(cfg_filename == NULL) XgError("Empty filename detected in XgInit().\nVerify number of command line parameters.");

  printf("\n");
  printf("-----------------------------------------\n");
  printf("  XGEN, a 2D interactive mesh generator  \n");
  printf("            Linux version 6.0            \n");
  printf("        Last modified Dec 5, 2003        \n");
  printf("    Revised and updated November 2010    \n");
  printf("   Copyright (C) 1995-2010 Pavel Solin   \n");
  printf("          e-mail: solin@unr.edu          \n");
  printf("    Distributed under the BSD license    \n");
  printf("-----------------------------------------\n");

  //initializing variables:
  Nfree = Npoin = First_Npoin = Nelem = Nbound = 0;
  Xmax = Xmin = Ymax = Ymin = 0;
  H = Area = 0;
  E = NULL;
  L = new LineList;
  EL = new PointList;
  BIL = new BoundaryInfoList;
  ElemL = new ElemList;
  IterationPtr = 0;
  RedrawFreePtr = 0;
  GoThroughPointsPtr = 0;
  TimestepConst = 0.8;

  //sorry, but:
  Dimension = 2;

  //getting the application name from the cfg filename
  char Help_str[50];
  strcpy(Help_str, cfg_filename);
  for(int i = 0; i<50; i++) {
    if((Help_str[i] == '\0') || (Help_str[i] == '.')) {
      if((Name = (char*)malloc((i+1)*sizeof(char))) == NULL) {
        XgError(
         "class Xgen",
         "Not enough memory for Name string." 
        );
      }
      Help_str[i] = '\0';
      strcpy(Name, Help_str);
      break;
    }
  }

  //case user wouldn't redefine:
  First_Nstore = Nstore = 10000;
  First_DeltaT = DeltaT = 1;
  H = 1;

  //executing user-redefined function:
  FILE *f; 

  f = fopen(cfg_filename, "rb");
  if(f == NULL) {
    XgError(
     "class Xgen",
     "Couldn't open configuration file." 
    );
  }        

  XgReadData(f);
  fclose(f); 

  //controlling the user:
  if(Nbound == 0) {
    XgError(
     "class Xgen",
     "Bad configuration file: Nbound = 0." 
    );
  }

  if(DeltaT <= 0) {
    XgError(
     "class Xgen",
     "DeltaT not initialized in XgReadData()." 
    );
  }
 
  //storing DeltaT case user would want to refresh it:
  First_DeltaT = DeltaT;

  //storing Nstore case user would want to refresh it:
  First_Nstore = Nstore;

  //computing H:
  H = Info.GiveBoundaryLength()/Nbound;

  //creating initial point list:
  EL = Info.Create_list_of_points();

  //getting extrems of boundary coordinates:
  Point T;
  EL->Init_ptr();
  Xmin = EL->First->P.x;            
  Xmax = EL->First->P.x;            
  Ymin = EL->First->P.y;            
  Ymax = EL->First->P.y;            
  while(EL->Get_ptr(&T)) {                      
    if(T.x < Xmin) Xmin = T.x;
    if(T.x > Xmax) Xmax = T.x;
    if(T.y < Ymin) Ymin = T.y;
    if(T.y > Ymax) Ymax = T.y;
  }
  EL->Init_ptr();

  //getting the region area:
  L = Info.Create_list_of_lines();
  int A, B;
  double xa, xb, ya, yb;
  Area = 0;                                                
  L->Init_ptr();                                           
  while(L->Get_ptr(&A, &B)) {                                    
    EL->Get(A, &xa, &ya);                                  
    EL->Get(B, &xb, &yb);                                  
    Area += (ya + yb - 2*Ymin)*(xa - xb)/2;                
  }                          
  L->Init_ptr();                                           

  //controlling the user:
  if(Area <= 0) {
    XgError(
     "class Xgen",
     "Bad boundary (or its orientation) in XgReadData()." 
    );
  }

  //getting optimal number of electrons that will be added:
  int help;
  help = (int)
  ((Area/(H*H*sqrt(3)/4) - 0.9*Nbound)/2 + 0.5);
  if(help <= 0) Nfree = 0;
  else Nfree = help;
  Npoin = Nfree + Nbound;

  //initializing iterations loop:
  IterationPtr = Nbound;

  //storing Npoin case user would want to refresh it:
  First_Npoin = Npoin;

  //copying user-added boundary electrons to the
  //beginning of electron list:
  E = (Point*)malloc((Npoin + Nstore)*sizeof(Point)); 
  if(E == NULL) {                                          
    XgError(
     "class Xgen",
     "Not enough memory for points." 
    );                                             
  }                                                      
  else {                                                 
    for(int i = 0; i<Nbound; i++) E[i] = EL->Delete();       
  }

  //free electrons are stochastic set:
  XgStartAgain();

  //at this moment class ready.
}

void Xgen::XgSetTimestep(double init_timestep) {
  DeltaT = init_timestep;
}

void Xgen::XgNewComponent() {
  Info.New_component();
}

void Xgen::XgAddEdge(int index, double x, double y, int lines) {
  Info.Add_edge(index, x, y, lines);
  Nbound += lines;
}

void Xgen::XgAddEdge(int index, double x, double y) {
  Info.Add_edge(index, x, y, 1);
  Nbound++;
}

void Xgen::XgAddEdge(int index, Point P, int lines) {
  Info.Add_edge(index, P.x, P.y, lines);
  Nbound += lines;
}

void Xgen::XgAddEdge(int index, Point P) {
  Info.Add_edge(index, P.x, P.y, 1);
  Nbound++;
}

/*
double Xgen::XgCriterion(Point p, Point q, Point r) {
  Point s = (p + q)/2;
  return (r - s).abs();
}
*/

double Xgen::XgCriterion(Point p, Point q, Point r) {
  return ((r-p)*(r-q))/((r-p).abs()*(r-q).abs());
}

int XgVertexInElem(int v, Element *e) {
  if(v == e->n1 || v == e->n2 || v == e->n3) return 1;
  else return 0;
}

void Xgen::XgUserOutput(FILE *f) {
  fprintf(f, "# XGEN mesh in Hermes2D format\n");
  fprintf(f, "# Project: "); fprintf(f, "%s\n", XgGiveName());
  fprintf(f, "# Edges are positively oriented\n");

  //PutNl(f, "* Number of points:");
  //PutNl(f, XgGiveNpoin());
  //PutNl(f, "* Number of elements:");
  //PutNl(f, XgGiveNelem());
  //PutNl(f, "* Boundary data number:");
  //PutNl(f, XgGiveNbound());

  //writing list of mesh vertices
  //(two coordinates per line)
  fprintf(f, "\n# Vertices:\nvertices =\n{\n");
  Point P;
  XgInitPointList();
  int counter = 0;
  while(XgGiveNextPoint(&P)) {
    counter++;
    if (counter < XgGiveNpoin()) fprintf(f, "  { %g, %g },\n", P.x, P.y);
    else fprintf(f, "  { %g, %g }\n", P.x, P.y);
  }
  fprintf(f, "}\n");

  //writing list of elements
  //(every element is defined via three
  //indices of corresponding grid points,
  // elements are positively oriented)
  fprintf(f, "\n# Elements:\nelements =\n{\n");
  Element E;
  XgInitElementList();
  counter = 0;
  while(XgGiveNextElement(&E)) {
    counter++;
    if (counter < XgGiveNelem()) fprintf(f, "  { %d, %d , %d, 0 },\n", E.n1, E.n2, E.n3);
    else fprintf(f, "  { %d, %d , %d, 0 }\n", E.n1, E.n2, E.n3);
  }
  fprintf(f, "}\n");

  //writing indices for all boundary edges
  //(their vertices are ordered according
  //to the orientation of the boundary)
  fprintf(f, "\n# Boundary data:\n");  
  fprintf(f, "# (bdy_vrt_1 bdy_vrt_2 edge_index)\nboundaries =\n{\n");  
  BoundaryInfo I;
  XgInitBoundaryInfoList();
  counter = 0;
  while(XgGiveNextBoundaryInfo(&I)) {
    counter++;
    if (counter < XgGiveNbound()) fprintf(f, "  { %d, %d , %d},\n", I.A, I.B, I.EdgeIndex);
    else fprintf(f, "  { %d, %d , %d}\n", I.A, I.B, I.EdgeIndex);
  }
  fprintf(f, "}\n");

  /*
  //writing boundary connectivity
  //(to every boundary edge the index of the adjacent
  //element and the index of the adjacent ghostcell)
  PutNl(f, "\n* Boundary connectivity:");  
  PutNl(f, "* (bdy_edge_nr adjacent_element_nr adjac_ghostcell_nr)");  
  int nelem = XgGiveNelem();
  BoundaryInfo Bin;
  XgInitBoundaryInfoList();
  int edge_count = 0;
  while(XgGiveNextBoundaryInfo(&Bin)) {
    edge_count++;
    Element Elem;
    XgInitElementList();
    int elem_count = 0;
    while(XgGiveNextElement(&Elem)) {
      elem_count++;
      int vrt_A = Bin.A;
      int vrt_B = Bin.B;
      if(XgVertexInElem(vrt_A, &Elem) && XgVertexInElem(vrt_B, &Elem)) {
        fprintf(f, "%ld %ld %ld\n", edge_count, elem_count, edge_count + nelem);
        break;
      }
    }
  }
  */
}

/**
 **   class Xgen: private methods
 **/

Point Xgen::Give_impuls(int j, int i) {
  double r, I;      //j to i
  //Point P;

  Point dr = E[i] - E[j];
  if(fabs(dr.x) > H || fabs(dr.y) > H) return Point(0, 0);
  r = dr.abs();
  if(r < XG_ZERO) return Point(1, 1);
  if(r >= H) return Point(0, 0);   //now r>XG_ZERO !
  I = H*H/(10*r);
  return Point(I*dr.x/r, I*dr.y/r);
}

Point Xgen::Get_impuls(int i) {
  double I, I0;
  Point Impuls(0, 0);
  for(int j=0; j<Npoin; j++) {
    if(j != i) Impuls += Give_impuls(j, i);
  }
  I = Impuls.abs();
  I0 = H/(sqr(DeltaT));
  if(I > I0) Impuls *= I0/I;
  return Impuls;
}

void Xgen::Shift(int i, Point Impuls) {
  Point Velocity(0, 0);
  Point Old_pos = E[i];

  Velocity = Impuls*DeltaT;
  E[i] += Velocity*(DeltaT*0.5);

  if(!Is_inside(E + i)) E[i] = Old_pos;
}

int Xgen::Is_inside(Point *P) {
  int Flag = 0;
  int a, b;

  L->Init_ptr();
  while((L->Get_ptr(&a, &b)) == 1) {
    Point A = E[a];
    Point B = E[b];
    if((A.x <= P->x) && (P->x < B.x)) {
      if(Is_right(A, B, *P)) Flag +=1;
    }
    if((A.x >= P->x) && (P->x > B.x)) {
      if(Is_left(A, B, *P)) Flag -=1;
    }
  }
  if(Flag == 0) return 0;
  else return 1;
}

int Xgen::Intact(
 Point a, Point b, Point c, Point d
) {
  double p, t, s;
  Point f = b - a, e = d - c;

  if(fabs(p = f.y*e.x - f.x*e.y) < XG_ZERO) return 0;

  t = (e.y*a.x - e.x*a.y + c.y*e.x - c.x*e.y)/p;
  if(fabs(e.x) > XG_ZERO) s = (a.x + t*f.x - c.x)/e.x;
  else s = (a.y + t*f.y - c.y)/e.y;
  if(t > XG_ZERO && t < 1 - XG_ZERO && s > XG_ZERO && s < 1 - XG_ZERO) 
   return 1;
  else return 0;
}

int Xgen::Boundary_intact(
 Point a, Point b, Point c
) {
  int d, e;
  L -> Init_ptr();
  while(L->Get_ptr(&d, &e)) {
    if(Intact(a, c, E[d], E[e]) ||
     Intact(b, c, E[d], E[e])) return 1;
  }
  return 0;
}

double AreaSize(Point a, Point b) {
  Point ab;
  ab = b - a;
  return ab.abs();
}

double Distance(Point P, Point a, Point b) {
  Point ab;
  double A, B, C, Nsize;
  ab = b - a;
  A = -ab.y;
  B = ab.x;
  C = -(A*a.x + B*a.y);
  Nsize = sqrt(A*A + B*B);
  return fabs((A*P.x + B*P.y + C)/Nsize);
}

double TriaVolume(Point a, Point b, Point c) {
  return Distance(a, b, c)*AreaSize(b, c)/2;
}

int Xgen::Find_nearest_left(
 int A, int B, 
 int *wanted
) {
  double m, Min;
  long int wanted0 = -1;
  Min = 1e50;
  for (int C = 0; C < Npoin; C++) {
    if (C != A && C != B) {
      if (Is_left(E[A], E[B], E[C])) {
        m = XgCriterion(E[A], E[B], E[C]);
        if(m < Min) {
          if(Is_inside(E + C) || L->Is_there(C)) {
            if(!Boundary_intact(E[A], E[B], E[C])) {
              if(TriaVolume(E[A], E[B], E[C]) > XG_ZERO) {
  	        Min = m;
	        wanted0 = C;
              }
            }
          }
        }
      }
    }
  }
  if(wanted0 < 0) return 0;
  *wanted = wanted0;
  return 1;
}

int Xgen::Is_left(Point A, Point B, Point C) {
  double x1 = B.x - A.x,
	 y1 = B.y - A.y,
	 x2 = C.x - A.x,
	 y2 = C.y - A.y;

  if((x1*y2 - y1*x2) > 0) return 1;
  else return 0;
}

int Xgen::Is_right(Point A, Point B, Point C) {
  double x1 = B.x - A.x,
	 y1 = B.y - A.y,
	 x2 = C.x - A.x,
	 y2 = C.y - A.y;

  if((x1*y2 - y1*x2) < 0) return 1;
  else return 0;
}

Xgen::~Xgen() {
  free(E);
  free(Name);
}

/**
 **    Motif interface:
 **/

# include <X11/Intrinsic.h>
# include <Xm/Xm.h>
# include <Xm/PushB.h>
# include <Xm/Form.h>

# include <Xm/Frame.h>
# include <Xm/ScrolledW.h>
# include <Xm/ScrollBar.h>

# include <Xm/DrawingA.h>
# include <Xm/Label.h>
# include <Xm/Text.h> 
# include <Xm/FileSB.h>
# include <Xm/RowColumn.h>
# include <Xm/ToggleB.h>
# include <Xm/CascadeB.h>

static struct {
  Xgen *RPtr;
  GC 
   draw_gc, 
   undraw_gc;
  Widget 
   topLevel,
   form,
   form1,  form2, 
   scrolled, horiz_bar,
   vert_bar, frame,
   menu_button, menu,
   about,
   timeinc, timedec,
   zoominc, zoomdec,
   savepoints,
   initforget, gridsave,
   quit, menucancel,
   bailout, info,
   blackboard;
  Widget 
   quitdialog,
   quitd_yes, quitd_no, quitd_label;
  Widget
   savepointsdialog;
  Widget 
   initdialog,
   initd_load, 
   initd_set, initd_cancel, initd_label,
   initd_load_d, initd_load_text,
   initd_load_ok, initd_load_cancel,
   initd_load_label;
  Widget 
   savedialog,
   saved_text, saved_ok,
   saved_cancel, saved_label;
  Widget 
   infodialog, info_l1,
   info_l2, info_l3,
   info_l4, info_l5,
   info_l6;
  XmString 
   info_l1_label1,
   info_l1_label2,
   info_l1_label3,
   info_l1_label4;
  char 
   info_npoin_str[50],
   info_nfree_str[50],
   info_nelem_str[50],
   info_timestep_str[50];
  XmString 
   menu_name, about_name,
   timeinc_name, timedec_name,
   zoominc_name, zoomdec_name,
   savepoints_name,
   forget_name, init_name, 
   grid_name, save_name, 
   quit_name, bailout_name,
   info_name;
  XmString 
   yes_label, no_label,
   ok_label,
   cancel_label;
  Pixel 
   black_color,
   white_color;
  int 
   Optimal_height,
   Optimal_width;
  int 
   Grid_working,
   Grid_failed,
   Grid_ready;
  Point 
   Ratio;
  int 
   MAX_X, 
   MAX_Y;
  int 
   Init_pos_x,
   Init_pos_y;
  int 
   SubLoopsNumber,
   InitSizeLimit,
   BoundaryRedrawInterval;
  char 
   PointPath[255], 
   GridPath[255];
  double 
   ZoomConst, Measure;
  Point 
   AreaSize, min, max;
  double 
   size_sliderx, size_slidery,
   pos_sliderx, pos_slidery; 
  int 
   max_sliderx, max_slidery;
  int 
   IsZoominc, IsInfo, IsAbout; 
} RS; 

void Draw_point(int blackboard_x, int blackboard_y) {
  XDrawArc(
   XtDisplay(RS.blackboard),
   XtWindow(RS.blackboard),
   RS.draw_gc,
   blackboard_x,
   blackboard_y,
   1, 1,
   0, 23040
 );
}

void Draw_point1(Point P) {
  int blackboard_x = -RS.Init_pos_x + (int)(RS.Ratio.x*P.x/RS.Measure),
   blackboard_y = RS.Init_pos_y + RS.MAX_Y - (int)(RS.Ratio.y*P.y/RS.Measure);
  Draw_point(blackboard_x, blackboard_y);
}

void Undraw_point(int blackboard_x, int blackboard_y) {
  XDrawArc(
   XtDisplay(RS.blackboard),
   XtWindow(RS.blackboard),
   RS.undraw_gc,
   blackboard_x,
   blackboard_y,
   1, 1,
   0, 23040
 );
}

void Undraw_point1(Point P) {
  int blackboard_x = - RS.Init_pos_x + (int)(RS.Ratio.x*P.x/RS.Measure),
   blackboard_y = RS.Init_pos_y + RS.MAX_Y - (int)(RS.Ratio.y*P.y/RS.Measure);
  Undraw_point(blackboard_x, blackboard_y);
}

void Draw_point2(Point P) {
  XDrawArc(
   XtDisplay(RS.blackboard),
   XtWindow(RS.blackboard),
   RS.draw_gc,
   -RS.Init_pos_x + (int)(RS.Ratio.x*P.x/RS.Measure + 0.5) - 1,
    RS.Init_pos_y + RS.MAX_Y - (int)(RS.Ratio.y*P.y/RS.Measure + 0.5) - 1,
    2, 2,
    0, 23040
  );
}

void Draw_line(Point A, Point B) {
  XDrawLine(
    XtDisplay(RS.blackboard),
    XtWindow(RS.blackboard),
    RS.draw_gc,
    -RS.Init_pos_x + (int)(RS.Ratio.x*A.x/RS.Measure + 0.5),
    RS.Init_pos_y + RS.MAX_Y - (int)(RS.Ratio.y*A.y/RS.Measure + 0.5),
    -RS.Init_pos_x + (int)(RS.Ratio.x*B.x/RS.Measure + 0.5),
    RS.Init_pos_y + RS.MAX_Y - (int)(RS.Ratio.y*B.y/RS.Measure + 0.5)
  );
}                             

void Graphics_const_init() {
  RS.Init_pos_x = 
   (int)(RS.pos_sliderx*RS.AreaSize.x*RS.Ratio.x/RS.max_sliderx/RS.Measure + 
   RS.min.x*RS.Ratio.x/RS.Measure + 0.5);
  RS.Init_pos_y = 
   (int)((RS.max_slidery - RS.size_slidery - RS.pos_slidery)
   *RS.AreaSize.y*RS.Ratio.y/RS.max_slidery/RS.Measure + 
   RS.min.y*RS.Ratio.y/RS.Measure + 0.5);
}

void Graphics_const_init(int iwidth, int iheight) {
  RS.MAX_X = iwidth - 1;
  RS.MAX_Y = iheight - 1;
  RS.Ratio.x = RS.MAX_X / RS.AreaSize.x;  
  RS.Ratio.y = RS.MAX_Y / RS.AreaSize.y;
  Graphics_const_init();
}

void Show_normal_situation() {
  Point a, b;             
           
  RS.RPtr->XgInitBoundaryLineList();                                           
  while(RS.RPtr->XgGiveNextBoundaryLine(&a, &b)) {
    Draw_point2(a);
    Draw_line(a, b);                       
  }                 
  RS.RPtr->XgInitFreePointList();
  while(RS.RPtr->XgGiveNextFreePoint(&a)) Draw_point1(a);
}

void Redraw_boundary() {
  Point a, b;             
           
  RS.RPtr->XgInitBoundaryLineList();                                           
  while(RS.RPtr->XgGiveNextBoundaryLine(&a, &b)) {
    Draw_point2(a);
    Draw_line(a, b);                       
  }                 
}

void Redraw_grid() {                                           
  Point a, b, c;                                               
  RS.RPtr->XgInitElementList();   
  while(RS.RPtr->XgGiveNextElement(&a, &b, &c)) {
    Draw_line(a, b);
    Draw_line(a, c);
    Draw_line(b, c);   
  }                 
}

void Redraw_blackboard() {
  if(RS.Grid_working || RS.Grid_ready) Redraw_grid();
  Show_normal_situation();
} 

void Get_optimal_size(int *width, int *height) {
  if(RS.AreaSize.x > RS.AreaSize.y) {
    *width = RS.InitSizeLimit;
    *height = (int)
     (RS.InitSizeLimit*RS.AreaSize.y/RS.AreaSize.x + 0.5);
  }
  else {
    *height = RS.InitSizeLimit;
    *width = (int)
     (RS.InitSizeLimit*RS.AreaSize.x/RS.AreaSize.y + 0.5);
  }
}   

void DeactivateWindow() {
  XtVaSetValues(
   RS.initforget,
   XmNsensitive, False,
   NULL
  );
  XtVaSetValues(
   RS.gridsave,
   XmNsensitive, False,
   NULL
  );     
  XtVaSetValues(
   RS.quit,
   XmNsensitive, False,
   NULL
  );     
  XtVaSetValues(
   RS.blackboard,
   XmNsensitive, False,
   NULL
  );
}

void ActivateWindow() {
  XtVaSetValues(
   RS.initforget,
   XmNsensitive, True,
   NULL
  );     
  XtVaSetValues(
   RS.gridsave,
   XmNsensitive, True,
   NULL
  );     
  XtVaSetValues(
   RS.quit,
   XmNsensitive, True,
   NULL
  );     
  XtVaSetValues(
   RS.blackboard,
   XmNsensitive, True,
   NULL
  );     
}

/*ARGSUSED*/
void QuitDialogYES(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  exit(0);
}

/*ARGSUSED*/
void QuitDialogNO(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtUnmanageChild(RS.quitdialog);
  ActivateWindow();
}

/*ARGSUSED*/
void InitDialogCANCEL(
 Widget w, XtPointer client_data, XtPointer call_data
) {  
  XtUnmanageChild(RS.initdialog);
  ActivateWindow();
}

/*ARGSUSED*/
void SaveDialogCANCEL(
 Widget w, XtPointer client_data, XtPointer call_data
) {  
  //w, client_data, call_data,
  XtUnmanageChild(RS.savedialog);
  ActivateWindow();
}

/*ARGSUSED*/
void SavePointsDialogCANCEL(
 Widget w, XtPointer client_data, XtPointer call_data
) {  
  XtUnmanageChild(RS.savepointsdialog);
  XtVaSetValues(
   RS.savepoints,
   XmNsensitive, True,
   NULL
  );  
}

/*ARGSUSED*/
void InitLoadDialogCANCEL(
 Widget w, XtPointer client_data, XtPointer call_data
) {  
  XtUnmanageChild(RS.initd_load_d);
  XtVaSetValues(
   RS.initd_load,
   XmNsensitive, True,
   NULL
  );
  XtVaSetValues(
   RS.initd_cancel,
   XmNsensitive, True,
   NULL
  );
  XtVaSetValues(
   RS.initd_set,
   XmNsensitive, True,
   NULL
  );
  XtVaSetValues(
   RS.initd_label,
   XmNsensitive, True,
   NULL
  );
}

void SavePointsErrorDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtUnmanageChild((Widget)client_data);

/* DO NOT CLEAR THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savepointsdialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, True,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savepointsdialog,
   XmNsensitive, True,
   NULL
  );
}

void InitDiskErrorDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  //w, client_data, call_data,
  XtUnmanageChild((Widget)client_data);

/* DO NOT CLEAR THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.initd_load_d, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, True,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.initd_load_d,
   XmNsensitive, True,
   NULL
  );
}

void SaveDiskErrorDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  //w, client_data, call_data,
  XtUnmanageChild((Widget)client_data);

/* DO NOT CLEAR THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, True,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savedialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, True,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savedialog,
   XmNsensitive, True,
   NULL
  );
}

void TrianglesErrorDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  //w, client_data, call_data,
  ActivateWindow();
  RS.Grid_failed = 0;
  char elem_info[50];
  strcpy(elem_info, RS.info_nelem_str);
  strcat(elem_info, "0");
  XmString info_l5_label = 
  XmStringCreate(elem_info, XmFONTLIST_DEFAULT_TAG);
  if(!RS.IsInfo) {
    XtVaSetValues(
     RS.info_l5,
     XmNlabelString, info_l5_label,
     NULL
    );
    XtVaSetValues(
     RS.info_l1,
     XmNlabelString, RS.info_l1_label1,
     NULL
    );
  }
  XtUnmanageChild((Widget)client_data);
  RS.RPtr->XgForgetGrid();
  XtVaSetValues(
   RS.initforget,
   XmNlabelString, RS.init_name,
   NULL
  );
  XtVaSetValues(
   RS.gridsave,
   XmNlabelString, RS.grid_name,
   NULL
  );   
  XClearWindow(
   XtDisplay(RS.blackboard),  
   XtWindow(RS.blackboard)    
  );    
  Show_normal_situation();
}

void SavePointsErrorDialog(char *event, char *where) {
/*
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savepointsdialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, False,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savepointsdialog,
   XmNsensitive, False,
   NULL
  );

  Widget errdialog, label1, label2, ok;

  int i = 0;
  Arg args[20];
  
  XtSetArg(args[i], XtNtitle, "error"); i++;
  XtSetArg(args[i], XmNminWidth, 350); i++;
  XtSetArg(args[i], XmNminHeight, 200); i++;
  XtSetArg(args[i], XmNwidth, 350); i++;
  XtSetArg(args[i], XmNheight, 200); i++;

  char msg[255];
  strcpy((char*)msg, "errdialog");
  errdialog = XmCreateFormDialog(
   RS.savepointsdialog, msg, args, i
  );

  XmString label1_label = 
   XmStringCreate(event, XmFONTLIST_DEFAULT_TAG);

  char str[255];
  strcpy(str, "'");
  strcat(str, where);
  strcat(str, "'");

  XmString label2_label = 
   XmStringCreate(str, XmFONTLIST_DEFAULT_TAG);

  label1 = XtVaCreateManagedWidget(
   "label1",
   xmLabelWidgetClass,
   errdialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 15,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 30,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label1_label,
   NULL
  ); 

  XtManageChild(label1);

  label2 = XtVaCreateManagedWidget(
   "label2",
   xmLabelWidgetClass,
   errdialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 30,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 45,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label2_label,
   NULL
  ); 

  XtManageChild(label2);

  ok = XtVaCreateManagedWidget(
   "Ok",
   xmPushButtonWidgetClass,
   errdialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 68,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 85,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 38,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 62,
   XmNlabelString, RS.ok_label,
   NULL
  );  
 
  XtManageChild(ok);

  XtAddCallback(ok, 
   XmNactivateCallback, SavePointsErrorDialogOK, (XtPointer)errdialog);

  XtManageChild(errdialog); 
}

void InitDiskErrorDialog(char *event, char *where) {
  Widget dialog, label1, label2, ok;

/* DO NOT CLEAR THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.initd_load_d, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, False,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.initd_load_d,
   XmNsensitive, False,
   NULL
  );

  Arg args[20];
  int i;

  i=0;
  XtSetArg(args[i], XtNtitle, "error"); i++;
  XtSetArg(args[i], XmNminWidth, 350); i++;
  XtSetArg(args[i], XmNminHeight, 200); i++;
  XtSetArg(args[i], XmNwidth, 350); i++;
  XtSetArg(args[i], XmNheight, 200); i++;

  char msg[255];
  strcpy((char*)msg, "Dialog");
  dialog = XmCreateFormDialog(
   RS.initd_load_d,
   msg,
   args,
   i
  );

  XmString label1_label = 
   XmStringCreate(event, XmFONTLIST_DEFAULT_TAG);

  char str[255];
  strcpy(str, "'");
  strcat(str, where);
  strcat(str, "'");

  XmString label2_label = 
   XmStringCreate(str, XmFONTLIST_DEFAULT_TAG);

  label1 = XtVaCreateManagedWidget(
   "label1",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 15,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 30,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label1_label,
   NULL
  ); 

  XtManageChild(label1);

  label2 = XtVaCreateManagedWidget(
   "label2",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 30,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 45,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label2_label,
   NULL
  ); 

  XtManageChild(label2);

  ok = XtVaCreateManagedWidget(
   "Ok",
   xmPushButtonWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 68,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 85,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 38,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 62,
   XmNlabelString, RS.ok_label,
   NULL
  );  
 
  XtManageChild(ok);

  XtAddCallback(ok, 
   XmNactivateCallback, InitDiskErrorDialogOK, (XtPointer)dialog);

  XtManageChild(dialog); 
}

void TrianglesErrorDialog() {
  Widget dialog, label1, label2, label3, ok;

  Arg args[20];
  int i;  

  i=0;
  XtSetArg(args[i], XtNtitle, "Message"); i++;
  XtSetArg(args[i], XmNminWidth, 350); i++;
  XtSetArg(args[i], XmNminHeight, 200); i++;
  XtSetArg(args[i], XmNwidth, 350); i++;
  XtSetArg(args[i], XmNheight, 200); i++;

  char msg[255];
  strcpy((char*)msg, "Dialog");
  dialog = XmCreateFormDialog(
   RS.topLevel,
   msg,
   args,
   i
  );

  strcpy((char*)msg, "Sorry, the algorithm failed.");
  XmString label1_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "Please try to modify");
  XmString label2_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "the parameter XG_ZERO in the source.");
  XmString label3_label = 
   XmStringCreate(msg, 
     XmFONTLIST_DEFAULT_TAG);

  label1 = XtVaCreateManagedWidget(
   "label1",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 10,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 25,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label1_label,
   NULL
  ); 

  XtManageChild(label1);

  label2 = XtVaCreateManagedWidget(
   "label2",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 35,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 44,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label2_label,
   NULL
  ); 

  XtManageChild(label2);

  label3 = XtVaCreateManagedWidget(
   "label3",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 45,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 54,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label3_label,
   NULL
  ); 

  XtManageChild(label3);

  ok = XtVaCreateManagedWidget(
   "Ok",
   xmPushButtonWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 68,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 85,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 38,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 62,
   XmNlabelString, RS.ok_label,
   NULL
  );  
 
  XtManageChild(ok);

  XtAddCallback(ok, 
   XmNactivateCallback, TrianglesErrorDialogOK, (XtPointer)dialog);

  XtManageChild(dialog); 
}

void SaveDiskErrorDialog(char *event, char *where) {
/*
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, False,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savedialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, False,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savedialog,
   XmNsensitive, False,
   NULL
  );

  Widget dialog, label1, label2, ok;

  Arg args[20];
  int i=0;  

  i=0;
  XtSetArg(args[i], XtNtitle, "error"); i++;
  XtSetArg(args[i], XmNminWidth, 350); i++;
  XtSetArg(args[i], XmNminHeight, 200); i++;
  XtSetArg(args[i], XmNwidth, 350); i++;
  XtSetArg(args[i], XmNheight, 200); i++;

  char msg[255];
  strcpy((char*)msg, "Dialog");
  dialog = XmCreateFormDialog(
   RS.savedialog,
   msg,
   args,
   i
  );

  XmString label1_label = 
   XmStringCreate(event, XmFONTLIST_DEFAULT_TAG);

  char str[255];
  strcpy(str, "'");
  strcat(str, where);
  strcat(str, "'");

  XmString label2_label = 
   XmStringCreate(str, XmFONTLIST_DEFAULT_TAG);

  label1 = XtVaCreateManagedWidget(
   "label1",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 15,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 30,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label1_label,
   NULL
  ); 

  XtManageChild(label1);

  label2 = XtVaCreateManagedWidget(
   "label2",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 30,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 45,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label2_label,
   NULL
  ); 

  XtManageChild(label2);

  ok = XtVaCreateManagedWidget(
   "Ok",
   xmPushButtonWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 68,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 85,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 38,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 62,
   XmNlabelString, RS.ok_label,
   NULL
  );  
 
  XtManageChild(ok);

  XtAddCallback(ok, 
   XmNactivateCallback, SaveDiskErrorDialogOK, (XtPointer)dialog);

  XtManageChild(dialog); 
}

/*ARGSUSED*/
void SavePointsDialogOK(Widget w, XtPointer client_data, XtPointer call_data) {
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_TEXT
  );  
  char *name = XmTextGetString(child);

  char msg[255];
  strcpy((char*)msg, "Couldn't create file");
  FILE *f = fopen(name, "wb");
  if(f == NULL) {
    SavePointsErrorDialog(msg, name);
  }
  else {
    RS.RPtr->XgOutputPoints(f);
    fclose(f);
    XtUnmanageChild(RS.savepointsdialog);
    XtVaSetValues(
     RS.savepoints,
     XmNsensitive, True,
     NULL
    );   
  }
}

/*ARGSUSED*/
void InitLoadDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_TEXT
  );
  
  char msg[255];
  strcpy((char*)msg, "Unable to open file");
  char *name = XmTextGetString(child);
  FILE *f = fopen(name, "rb");
  if(f == NULL) {
    InitDiskErrorDialog(msg, name);
  }
  else {
    int Out, Error, NoMatch, Mem;
    Out = 0; NoMatch = 0;
    RS.RPtr->XgInputPoints(f, &Error, &Mem);
    if(Out) {  
      strcpy((char*)msg, "A point out of the area found in file");
      InitDiskErrorDialog(
       msg, name);
    }
    else {
      if(Error) {
        strcpy((char*)msg, "Error while reading in file");
        InitDiskErrorDialog(
         msg, name);
      }
      else {  
        if(NoMatch) {
          strcpy((char*)msg, "A value doesn't match in file");
          InitDiskErrorDialog(
           msg, name);
        }
        else {
          if(Mem) {
            strcpy((char*)msg, "Not enough memory for points from");
            InitDiskErrorDialog(
             msg, name);
          }
          else {
            ActivateWindow();   
            int npoin = RS.RPtr->XgGiveNpoin();
            char points_info[255];
            strcpy(points_info, RS.info_npoin_str);
            char number[20];
            sprintf(number, "%d", npoin);
            strcat(points_info, number);
            XmString info_l2_label = 
             XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
            int nfree = RS.RPtr->XgGiveNfree();
            strcpy(points_info, RS.info_nfree_str);
            sprintf(number, "%d", nfree);
            strcat(points_info, number);
            XmString info_l3_label = 
             XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
            if(!RS.IsInfo) {
              XtVaSetValues(
               RS.info_l2,
               XmNlabelString, info_l2_label,
               NULL
              );
              XtVaSetValues(
               RS.info_l3,
               XmNlabelString, info_l3_label,
               NULL
              );
            }
            XtUnmanageChild(RS.initd_load_d);
            XtUnmanageChild(RS.initdialog); 
            XClearWindow(
             XtDisplay(RS.blackboard),  
             XtWindow(RS.blackboard)    
            );    
            Show_normal_situation();
          }
        } 
      }
    }
    fclose(f); 
  }
}

/*ARGSUSED*/
void SaveDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {;
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_TEXT
  );
  
  char *name = XmTextGetString(child);
  FILE *f = fopen(name, "wb");
  if(f == NULL) {
    char msg[255];
    strcpy((char*)msg, "Unable to create file");
    SaveDiskErrorDialog(msg, name);
  }
  else {
    ActivateWindow();   
    RS.RPtr->XgOutput(f);
    fclose(f);
    RS.RPtr->XgForgetGrid();
    RS.Grid_ready = 0;
    char elem_info[50];
    strcpy(elem_info, RS.info_nelem_str);
    strcat(elem_info, "0");
    XmString info_l5_label = 
    XmStringCreate(elem_info, XmFONTLIST_DEFAULT_TAG);
    if(!RS.IsInfo) {
      XtVaSetValues(
       RS.info_l5,
       XmNlabelString, info_l5_label,
       NULL
      );
      XtVaSetValues(
       RS.info_l1,
       XmNlabelString, RS.info_l1_label1,
       NULL
      );
    }
    XtVaSetValues(
     RS.initforget,
     XmNlabelString, RS.init_name,
     NULL
    );
    XtVaSetValues(
     RS.gridsave,
     XmNlabelString, RS.grid_name,
     NULL
    );
    XtUnmanageChild(RS.savedialog);   
    XClearWindow(
     XtDisplay(RS.blackboard),  
     XtWindow(RS.blackboard)    
    );    
    Show_normal_situation();
  }
}

/*ARGSUSED*/
void SavePoints(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  Arg args[20]; int i=0;  
  XtSetArg(args[i], XtNtitle, "Request"); i++;

  XtSetArg(args[i], XmNfileTypeMask, XmFILE_REGULAR); i++;  
  XtSetArg(args[i], XmNminWidth, 310); i++;  
  XtSetArg(args[i], XmNminHeight, 430); i++;  
  XtSetArg(args[i], XmNwidth, 310); i++;  
  XtSetArg(args[i], XmNheight, 430); i++;  

  XmString label;

  //is there the points directory?
  char test_filename[100];
  strcpy(test_filename, RS.PointPath);
  strcat(test_filename, "test.test");
  FILE *f = fopen(test_filename, "wb");
  if(f != NULL) {
    label = 
     XmStringCreate(RS.PointPath, XmFONTLIST_DEFAULT_TAG);
    XtSetArg(args[i], XmNdirectory, label); i++;
    fclose(f);
    unlink(test_filename);
  }

  char msg[255];
  strcpy((char*)msg, "The points file name:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);  
  XtSetArg(args[i], XmNselectionLabelString, label); i++;
  strcpy((char*)msg, "Filter:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNfilterLabelString, label); i++;
  strcpy((char*)msg, "Directory:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNdirListLabelString, label); i++;
  strcpy((char*)msg, "File:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNlistLabelString, label); i++;
  strcpy((char*)msg, "(empty)");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNnoMatchString, label); i++;
  strcpy((char*)msg, "*");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);   
  XtSetArg(args[i], XmNpattern, label); i++;

  strcpy((char*)msg, "savepointsdialog");
  RS.savepointsdialog = XmCreateFileSelectionDialog(
   RS.topLevel,
   msg,
   args,
   i
  );

  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_HELP_BUTTON
  );
  XtUnmanageChild(child);

  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SEPARATOR
  ); 
  XtUnmanageChild(child);

  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_OK_BUTTON
  ); 
  strcpy((char*)msg, "Ok");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_APPLY_BUTTON
  ); 
  strcpy((char*)msg, "Filter");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_CANCEL_BUTTON
  ); 
  strcpy((char*)msg, "Cancel");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );

  XtVaSetValues(
   RS.savepoints,
   XmNsensitive, False,
   NULL
  );

  XtAddCallback(RS.savepointsdialog, 
   XmNokCallback, SavePointsDialogOK, NULL);
  XtAddCallback(RS.savepointsdialog, 
   XmNcancelCallback, SavePointsDialogCANCEL, NULL); 

  XtManageChild(RS.savepointsdialog);
}

/*ARGSUSED*/
void InitDialogLOAD(
 Widget w, XtPointer client_data, XtPointer call_data
) {  
  Arg args[20]; int i=0;  
  XtSetArg(args[i], XtNtitle, "Request"); i++;

  XtSetArg(args[i], XmNfileTypeMask, XmFILE_REGULAR); i++;  
  XtSetArg(args[i], XmNminWidth, 310); i++;  
  XtSetArg(args[i], XmNminHeight, 430); i++;  
  XtSetArg(args[i], XmNwidth, 310); i++;  
  XtSetArg(args[i], XmNheight, 430); i++;  

  XmString label;

  //is there the points directory?
  char test_filename[100];
  strcpy(test_filename, RS.PointPath);
  strcat(test_filename, "test.test");
  FILE *f = fopen(test_filename, "wb");
  if(f != NULL) {
    label = 
     XmStringCreate(RS.PointPath, XmFONTLIST_DEFAULT_TAG);
    XtSetArg(args[i], XmNdirectory, label); i++;
    fclose(f);
    unlink(test_filename);
  }

  char msg[255];
  strcpy((char*)msg, "The points file name:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);  
  XtSetArg(args[i], XmNselectionLabelString, label); i++;
  strcpy((char*)msg, "Filter:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNfilterLabelString, label); i++;
  strcpy((char*)msg, "Directory:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNdirListLabelString, label); i++;
  strcpy((char*)msg, "File:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNlistLabelString, label); i++;
  strcpy((char*)msg, "(empty)");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNnoMatchString, label); i++;
  strcpy((char*)msg, "*.*");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);   
  XtSetArg(args[i], XmNpattern, label); i++;

  strcpy((char*)msg, "initd_load_d");
  RS.initd_load_d = XmCreateFileSelectionDialog(
   RS.initdialog,
   msg,
   args,
   i
  );

  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_HELP_BUTTON
  );
  XtUnmanageChild(child);

  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SEPARATOR
  ); 
  XtUnmanageChild(child);

  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_OK_BUTTON
  ); 
  strcpy((char*)msg, "Ok");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_APPLY_BUTTON
  ); 
  strcpy((char*)msg, "Filter");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_CANCEL_BUTTON
  ); 
  strcpy((char*)msg, "Cancel");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );

  XtAddCallback(RS.initd_load_d, 
   XmNokCallback, InitLoadDialogOK, NULL);
  XtAddCallback(RS.initd_load_d, 
   XmNcancelCallback, InitLoadDialogCANCEL, NULL); 

  XtManageChild(RS.initd_load_d); 

  XtVaSetValues(
   RS.initd_load,
   XmNsensitive, False,
   NULL
  );
  XtVaSetValues(
   RS.initd_set,
   XmNsensitive, False,
   NULL
  );
  XtVaSetValues(
   RS.initd_cancel,
   XmNsensitive, False,
   NULL
  );
  XtVaSetValues(
   RS.initd_label,
   XmNsensitive, False,
   NULL
  );
}

/*ARGSUSED*/
void SaveDialog() {  
  Arg args[20]; int i=0;  
  XtSetArg(args[i], XtNtitle, "Request"); i++;
  XtSetArg(args[i], XmNfileTypeMask, XmFILE_REGULAR); i++;  
  XtSetArg(args[i], XmNminWidth, 310); i++;  
  XtSetArg(args[i], XmNminHeight, 430); i++;  
  XtSetArg(args[i], XmNwidth, 310); i++;  
  XtSetArg(args[i], XmNheight, 430); i++;  

  XmString label;

  //is there the grid directory?
  char test_filename[100];
  strcpy(test_filename, RS.GridPath);
  strcat(test_filename, "test.test");
  FILE *f = fopen(test_filename, "wb");
  if(f != NULL) {
    label = 
     XmStringCreate(RS.GridPath, XmFONTLIST_DEFAULT_TAG);
    XtSetArg(args[i], XmNdirectory, label); i++;
    fclose(f);
    unlink(test_filename);
  }

  char msg[255];
  strcpy((char*)msg, "Mesh file name:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);  
  XtSetArg(args[i], XmNselectionLabelString, label); i++;
  strcpy((char*)msg, "Filter:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNfilterLabelString, label); i++;
  strcpy((char*)msg, "Directory:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNdirListLabelString, label); i++;
  strcpy((char*)msg, "File:");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNlistLabelString, label); i++;
  strcpy((char*)msg, "(empty)");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtSetArg(args[i], XmNnoMatchString, label); i++;
  strcpy((char*)msg, "*.*");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);   
  XtSetArg(args[i], XmNpattern, label); i++;

  strcpy((char*)msg, "savedialog");
  RS.savedialog = XmCreateFileSelectionDialog(
   RS.topLevel,
   msg,
   args,
   i
  );

  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_HELP_BUTTON
  );
  XtUnmanageChild(child);

  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SEPARATOR
  ); 
  XtUnmanageChild(child);

  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_OK_BUTTON
  ); 
  strcpy((char*)msg, "Ok");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_APPLY_BUTTON
  ); 
  strcpy((char*)msg, "Filter");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_CANCEL_BUTTON
  ); 
  strcpy((char*)msg, "Cancel");
  label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);    
  XtVaSetValues(
   child,
   XmNdefaultButtonShadowThickness, 0,
   XmNlabelString, label,
   NULL
  );

  XtAddCallback(RS.savedialog, 
   XmNokCallback, SaveDialogOK, NULL);
  XtAddCallback(RS.savedialog, 
   XmNcancelCallback, SaveDialogCANCEL, NULL); 

  XtManageChild(RS.savedialog); 
}

/*ARGSUSED*/
void InitDialogSET(
 Widget w, XtPointer client_data, XtPointer call_data
) {  
  ActivateWindow();
  RS.RPtr->XgStartAgain2();
  int npoin = RS.RPtr->XgGiveNpoin();
  char points_info[255];
  strcpy(points_info, RS.info_npoin_str);
  char number[20];
  sprintf(number, "%d", npoin);
  strcat(points_info, number);
  XmString info_l2_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
  int nfree = RS.RPtr->XgGiveNfree();
  strcpy(points_info, RS.info_nfree_str);
  sprintf(number, "%d", nfree);
  strcat(points_info, number);
  XmString info_l3_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
  if(!RS.IsInfo) {
    XtVaSetValues(
     RS.info_l2,
     XmNlabelString, info_l2_label,
     NULL
    );
    XtVaSetValues(
     RS.info_l3,
     XmNlabelString, info_l3_label,
     NULL
    );
  }
  XtUnmanageChild(RS.initdialog);
  XClearWindow(
    XtDisplay(RS.blackboard),  
    XtWindow(RS.blackboard)    
  );    
  Show_normal_situation();
}

void QuitDialog(XtPointer client_data) {
  //quit dialog:
  Arg args[20];
  int i;  

  i=0;
  XtSetArg(args[i], XtNtitle, "Question"); i++;
  XtSetArg(args[i], XmNminWidth, 400); i++;
  XtSetArg(args[i], XmNminHeight,100); i++;
  XtSetArg(args[i], XmNwidth, 400); i++;
  XtSetArg(args[i], XmNheight,100); i++;

  char msg[255];
  strcpy((char*)msg, "quitdialog");
  RS.quitdialog = XmCreateFormDialog(
   RS.topLevel,
   msg,
   args,
   i
  );

  Widget subdialog, top;

  top = XtVaCreateManagedWidget(
   "top",
   xmFormWidgetClass,
   RS.quitdialog,
   XmNtopAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNbottomOffset, 39,
   NULL
  ); 

  XtManageChild(top);

  subdialog = XtVaCreateManagedWidget(
   "subdialog",
   xmFormWidgetClass,
   RS.quitdialog,
   XmNtopAttachment, XmATTACH_WIDGET,
   XmNtopWidget, top,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNbottomOffset, 5,
   XmNleftAttachment, XmATTACH_FORM,
   XmNleftOffset, 0,
   XmNrightAttachment, XmATTACH_FORM,
   XmNrightOffset, 0,
   NULL
  ); 

  XtManageChild(subdialog);

  strcpy((char*)msg, "Do you really want to quit?");
  XmString label_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  RS.quitd_label = XtVaCreateManagedWidget(
   "quitd_label",
   xmLabelWidgetClass,
   top,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label_label,
   NULL
  ); 

  RS.quitd_yes = XtVaCreateManagedWidget(
   "quitd_yes",
   xmPushButtonWidgetClass,
   subdialog,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 4,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 34,
   XmNlabelString, RS.yes_label,
   NULL
  ); 

  RS.quitd_no = XtVaCreateManagedWidget(
   "quitd_no",
   xmPushButtonWidgetClass,
   subdialog,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 66,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 96,
   XmNlabelString, RS.no_label,
   NULL
  ); 

  XtManageChild(RS.quitd_label);
  XtManageChild(RS.quitd_yes);
  XtManageChild(RS.quitd_no);

  XtAddCallback(RS.quitd_yes, 
   XmNactivateCallback, QuitDialogYES, NULL);
  XtAddCallback(RS.quitd_no, 
   XmNactivateCallback, QuitDialogNO, NULL);

  XtManageChild(RS.quitdialog);
  DeactivateWindow();
}

void InitDialog() {
  //init dialog:
  Arg args[20];
  int i;  

  i=0;
  XtSetArg(args[i], XtNtitle, "question"); i++;
  XtSetArg(args[i], XmNminWidth, 400); i++;
  XtSetArg(args[i], XmNminHeight,100); i++;
  XtSetArg(args[i], XmNwidth, 400); i++;
  XtSetArg(args[i], XmNheight,100); i++;

  char msg[255];
  strcpy((char*)msg, "initdialog");
  RS.initdialog = XmCreateFormDialog(
   RS.form,
   msg,
   args,
   i
  );

  Widget subdialog, top;

  top = XtVaCreateManagedWidget(
   "top",
   xmFormWidgetClass,
   RS.initdialog,
   XmNtopAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNbottomOffset, 39,
   NULL
  ); 

  XtManageChild(top);

  subdialog = XtVaCreateManagedWidget(
   "subdialog",
   xmFormWidgetClass,
   RS.initdialog,
   XmNtopAttachment, XmATTACH_WIDGET,
   XmNtopWidget, top,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNbottomOffset, 5,
   XmNleftAttachment, XmATTACH_FORM,
   XmNleftOffset, 0,
   XmNrightAttachment, XmATTACH_FORM,
   XmNrightOffset, 0,
   NULL
  ); 

  XtManageChild(subdialog);

  strcpy((char*)msg, "Load points or overlay them?");
  XmString label_label = 
   XmStringCreate(
    msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Load");
  XmString load_str = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Overlay");
  XmString set_str = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "initd_label");
  RS.initd_label = XtVaCreateManagedWidget(
   msg,
   xmLabelWidgetClass,
   top,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label_label,
   NULL
  ); 

  strcpy((char*)msg, "initd_load");
  RS.initd_load = XtVaCreateManagedWidget(
   msg,
   xmPushButtonWidgetClass,
   subdialog,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 4,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 24,
   XmNlabelString, load_str,
   NULL
  ); 

  strcpy((char*)msg, "initd_set");
  RS.initd_set = XtVaCreateManagedWidget(
   msg,
   xmPushButtonWidgetClass,
   subdialog,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 25,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 45,
   XmNlabelString, set_str,
   NULL
  ); 
 
  strcpy((char*)msg, "initd_cancel");
  RS.initd_cancel = XtVaCreateManagedWidget(
   msg,
   xmPushButtonWidgetClass,
   subdialog,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 66,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 96,
   XmNlabelString, RS.cancel_label,
   NULL
  ); 

  XtManageChild(RS.initd_label);
  XtManageChild(RS.initd_load);
  XtManageChild(RS.initd_set);
  XtManageChild(RS.initd_cancel);

  XtAddCallback(RS.initd_load, 
   XmNactivateCallback, InitDialogLOAD, NULL);
  XtAddCallback(RS.initd_set, 
   XmNactivateCallback, InitDialogSET, NULL);
  XtAddCallback(RS.initd_cancel, 
   XmNactivateCallback, InitDialogCANCEL, NULL);

  XtManageChild(RS.initdialog);
  DeactivateWindow();
}

/*
 *   Callback functions:
 */

void AboutOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  RS.IsAbout = 1;
  //w, client_data, call_data,
  XtUnmanageChild((Widget)client_data);

  XtVaSetValues(
   RS.about,
   XmNsensitive, True,
   NULL
  );
}

/*ARGSUSED*/
void About(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  RS.IsAbout = 0;

  XtVaSetValues(
   RS.about,
   XmNsensitive, False,
   NULL
  );

  Arg args[20];
  int i;

  i=0;
  XtSetArg(args[i], XtNtitle, "Message"); i++;
  XtSetArg(args[i], XmNminWidth, 310); i++;
  XtSetArg(args[i], XmNminHeight, 250); i++;
  XtSetArg(args[i], XmNwidth, 310); i++;
  XtSetArg(args[i], XmNheight, 250); i++;

  char msg[255];
  strcpy((char*)msg, "Dialog");
  Widget dialog = XmCreateFormDialog(
   RS.topLevel,
   msg,
   args,
   i
  );

  strcpy((char*)msg, "XGEN");
  XmString labela_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "2D interactive mesh generator");
  XmString labelb_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "Version 6.0");
  XmString label1_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "Last changes: December 1998");
  XmString label2_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "Revised and updated: November 2010");
  XmString label3_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "Contact email:");
  XmString label4_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  strcpy((char*)msg, "solin@unr.edu");
  XmString label5_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  Widget labela = XtVaCreateManagedWidget(
   "labela",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 5,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 13,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, labela_label,
   NULL
  ); 

  XtManageChild(labela);

  Widget labelb = XtVaCreateManagedWidget(
   "labelb",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 12,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 20,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, labelb_label,
   NULL
  ); 

  XtManageChild(labelb);

  Widget label1 = XtVaCreateManagedWidget(
   "label1",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 25,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 33,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label1_label,
   NULL
  ); 

  XtManageChild(label1);

  Widget label2 = XtVaCreateManagedWidget(
   "label2",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 32,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 40,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label2_label,
   NULL
  ); 

  XtManageChild(label2);

  Widget label3 = XtVaCreateManagedWidget(
   "label3",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 45,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 53,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label3_label,
   NULL
  ); 

  XtManageChild(label3);

  Widget label4 = XtVaCreateManagedWidget(
   "label4",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 58,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 66,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label4_label,
   NULL
  ); 

  XtManageChild(label4);

  Widget label5 = XtVaCreateManagedWidget(
   "label5",
   xmLabelWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 65,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 73,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, label5_label,
   NULL
  ); 

  XtManageChild(label5);

  Widget ok = XtVaCreateManagedWidget(
   "Ok",
   xmPushButtonWidgetClass,
   dialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 80,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 92,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 40,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 60,
   XmNlabelString, RS.ok_label,
   NULL
  );  
 
  XtManageChild(ok);

  XtAddCallback(ok, 
   XmNactivateCallback, AboutOK, (XtPointer)dialog);

  XtManageChild(dialog); 
}

/*ARGSUSED*/
void TimeInc(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  RS.RPtr -> XgTimeInc();
  char time_info[50], value[50];
  strcpy(time_info, RS.info_timestep_str);
  sprintf(value, "%g", RS.RPtr->XgGiveTimestep());
  strcat(time_info, value);
  XmString info_l6_label = 
   XmStringCreate(time_info, XmFONTLIST_DEFAULT_TAG);
  if(!RS.IsInfo) {
    XtVaSetValues(
     RS.info_l6,
     XmNlabelString, info_l6_label,
     NULL
    );
  }
}

/*ARGSUSED*/
void TimeDec(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  RS.RPtr -> XgTimeDec();
  char time_info[50], value[50];
  strcpy(time_info, RS.info_timestep_str);
  sprintf(value, "%g", RS.RPtr->XgGiveTimestep());
  strcat(time_info, value);
  XmString info_l6_label = 
   XmStringCreate(time_info, XmFONTLIST_DEFAULT_TAG);
  if(!RS.IsInfo) {
    XtVaSetValues(
     RS.info_l6,
     XmNlabelString, info_l6_label,
     NULL
    );
  }
}

/*ARGSUSED*/
void ShiftX(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  int pos;
  XtVaGetValues(
   RS.horiz_bar,
   XmNvalue, &pos,
   NULL
  );
  RS.pos_sliderx = pos;
  RS.Init_pos_x = 
   (int)(RS.pos_sliderx*RS.AreaSize.x*RS.Ratio.x/RS.max_sliderx/RS.Measure + 
   RS.min.x*RS.Ratio.x/RS.Measure + 0.5);
  XClearWindow(
   XtDisplay(RS.blackboard),  
   XtWindow(RS.blackboard)    
  );    
  Redraw_blackboard();
}

/*ARGSUSED*/
void ShiftY(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  int pos;
  XtVaGetValues(
   RS.vert_bar,
   XmNvalue, &pos,
   NULL
  );
  RS.pos_slidery = pos;
  RS.Init_pos_y = 
   (int)((RS.max_slidery - RS.size_slidery - RS.pos_slidery)
   *RS.AreaSize.y*RS.Ratio.y/RS.max_slidery/RS.Measure + 
   RS.min.y*RS.Ratio.y/RS.Measure + 0.5);
  XClearWindow(
   XtDisplay(RS.blackboard),  
   XtWindow(RS.blackboard)    
  );    
  Redraw_blackboard();
}

/*ARGSUSED*/
void ZoomInc(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  if(!RS.IsZoominc) return;
  int size;
  RS.Measure *= RS.ZoomConst;
  if(RS.Ratio.x/RS.Measure*RS.AreaSize.x > 50000 ||
   RS.Ratio.y/RS.Measure*RS.AreaSize.y > 50000) {
    RS.IsZoominc = 0;
    XtVaSetValues(
     RS.zoominc,
     XmNsensitive, False,
     NULL
    );
  }
  double old_size = RS.size_sliderx;  
  RS.size_sliderx = RS.Measure*RS.max_sliderx;
  if((int)(RS.size_sliderx+0.5) < RS.max_sliderx) {
    size = (int)(RS.size_sliderx + 0.5);
    if(size < 1) size = 1;
    XtVaSetValues(
     RS.horiz_bar,
     XmNsliderSize, size,
     NULL
    );
  }
  RS.pos_sliderx += (old_size - RS.size_sliderx)/2;
  if((int)(RS.pos_sliderx+0.5) > 0) {
    XtVaSetValues(
     RS.horiz_bar,
     XmNvalue, (int)(RS.pos_sliderx+0.5),
     NULL
    );
  }
  old_size = RS.size_slidery;  
  RS.size_slidery = RS.Measure*RS.max_slidery;
  if((int)(RS.size_slidery+0.5) < RS.max_slidery) {
    size = (int)(RS.size_slidery + 0.5);
    if(size < 1) size = 1;
    XtVaSetValues(
     RS.vert_bar,
     XmNsliderSize, size,
     NULL
    );
  }
  RS.pos_slidery += (old_size - RS.size_slidery)/2;
  if((int)(RS.pos_slidery+0.5) > 0) {
    XtVaSetValues(
     RS.vert_bar,
     XmNvalue, (int)(RS.pos_slidery+0.5),
     NULL
    );
  }
  Graphics_const_init();
  XClearWindow(
   XtDisplay(RS.blackboard),  
   XtWindow(RS.blackboard)    
  );    
  Redraw_blackboard();
}

/*ARGSUSED*/
void ZoomDec(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  RS.Measure /= RS.ZoomConst;
  if(!RS.IsZoominc) {
    RS.IsZoominc = 1;
    XtVaSetValues(
     RS.zoominc,
     XmNsensitive, True,
     NULL
    );
  }
  int size;
  double old_size = RS.size_sliderx;
  RS.size_sliderx = RS.Measure*RS.max_sliderx;
  double old_pos = RS.pos_sliderx;
  RS.pos_sliderx -= (RS.size_sliderx - old_size)/2;
  if((int)(old_size+0.5) < RS.max_sliderx) {
    if((int)(RS.size_sliderx+0.5) > RS.max_sliderx) {
      size = RS.max_sliderx;
      RS.pos_sliderx = 0;
    }
    else {
      size = (int)(RS.size_sliderx+0.5);
      if((int)(RS.pos_sliderx+0.5) < 0) RS.pos_sliderx = 0;
      else {
        if((int)(RS.pos_sliderx+0.5) + 
         (int)(RS.size_sliderx+0.5) > RS.max_sliderx)
        RS.pos_sliderx = RS.max_sliderx - (int)(RS.size_sliderx + 0.5);
      }  
    }
    XtVaSetValues(
     RS.horiz_bar,
     XmNvalue, (int)(RS.pos_sliderx+0.5),
     NULL
    ); 
    if(size < 1) size = 1;
    XtVaSetValues(
     RS.horiz_bar,
     XmNsliderSize, size,
     NULL
    );
  }
  old_size = RS.size_slidery;
  RS.size_slidery = RS.Measure*RS.max_slidery;
  old_pos = RS.pos_slidery;
  RS.pos_slidery -= (RS.size_slidery - old_size)/2;
  if((int)(old_size+0.5) < RS.max_slidery) {
    if((int)(RS.size_slidery+0.5) > RS.max_slidery) {
      size = RS.max_slidery;
      RS.pos_slidery = 0;
    }
    else {
      size = (int)(RS.size_slidery+0.5);
      if((int)(RS.pos_slidery+0.5) < 0) RS.pos_slidery = 0;
      else {
        if((int)(RS.pos_slidery+0.5) + 
         (int)(RS.size_slidery+0.5) > RS.max_slidery)
        RS.pos_slidery = RS.max_slidery - (int)(RS.size_slidery + 0.5);
      }  
    }
    XtVaSetValues(
     RS.vert_bar,
     XmNvalue, (int)(RS.pos_slidery+0.5),
     NULL
    ); 
    if(size < 1) size = 1;
    XtVaSetValues(
     RS.vert_bar,
     XmNsliderSize, size,
     NULL
    );
  }
  Graphics_const_init();  
  XClearWindow(
   XtDisplay(RS.blackboard),  
   XtWindow(RS.blackboard)    
  );    
  Redraw_blackboard();
}

/*ARGSUSED*/
void InitForget(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  if(RS.Grid_ready || RS.Grid_working) {     
    RS.RPtr->XgForgetGrid();
    RS.Grid_working = 0; 
    RS.Grid_ready = 0;
    char elem_info[50];
    strcpy(elem_info, RS.info_nelem_str);
    strcat(elem_info, "0");
    XmString info_l5_label = 
     XmStringCreate(elem_info, XmFONTLIST_DEFAULT_TAG);
    if(!RS.IsInfo) {
      XtVaSetValues(
       RS.info_l5,
       XmNlabelString, info_l5_label,
       NULL
      );
      XtVaSetValues(
       RS.info_l1,
       XmNlabelString, RS.info_l1_label1,
       NULL
      );
    }        
    XtVaSetValues(
     RS.initforget,
     XmNlabelString, RS.init_name,
     NULL
    );     
    XtVaSetValues(
     RS.gridsave,
     XmNsensitive, True,
     XmNlabelString, RS.grid_name,
     NULL
    );     
    XtVaSetValues(
     RS.blackboard,
     XmNsensitive, True,
     NULL
    );
    XClearWindow(
      XtDisplay(RS.blackboard),  
      XtWindow(RS.blackboard)    
    );    
    Show_normal_situation();
  } 
  else InitDialog();     
}

/*ARGSUSED*/
void GridSave(
 Widget w, XtPointer client_data, XtPointer call_data
) { 
  if(RS.Grid_ready) {
    SaveDialog();
  }
  else {     
    if(!RS.IsInfo) {     
      XtVaSetValues(
       RS.info_l1,
       XmNlabelString, RS.info_l1_label2,
       NULL
      );
    }
    RS.Grid_working = 1;          
    XtVaSetValues(
     RS.initforget,
     XmNlabelString, RS.forget_name,
     NULL
    );          
    XtVaSetValues(
     RS.gridsave,
     XmNsensitive, False,
     NULL
    );          
    XtVaSetValues(
     RS.blackboard,
     XmNsensitive, False,
     NULL
    );          
  }
}

/*ARGSUSED*/
void Quit(
 Widget w, XtPointer client_data, XtPointer call_data
) {  
  QuitDialog(client_data);
}

/*ARGSUSED*/
static void BlackboardInput(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  int blackboard_x, blackboard_y;
  XmDrawingAreaCallbackStruct 
   *draw_data = (XmDrawingAreaCallbackStruct*)call_data;
  if(draw_data -> event -> type != ButtonRelease) return;
  blackboard_x = draw_data -> event -> xbutton.x; 
  blackboard_y = draw_data -> event -> xbutton.y;
  Point P;
  P.x =  RS.Measure*(blackboard_x + RS.Init_pos_x)/RS.Ratio.x;
  P.y = RS.Measure*(RS.MAX_Y +  RS.Init_pos_y - blackboard_y)/RS.Ratio.y;
  if(draw_data -> event -> xbutton.button == Button1) { 
    if(RS.RPtr -> XgMouseAdd(P)) {
      Draw_point1(P);
      int npoin = RS.RPtr->XgGiveNpoin();
      char points_info[255];
      strcpy(points_info, RS.info_npoin_str);
      char number[20];
      sprintf(number, "%d", npoin);
      strcat(points_info, number);
      XmString info_l2_label = 
       XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
      int nfree = RS.RPtr->XgGiveNfree();
      strcpy(points_info, RS.info_nfree_str);
      sprintf(number, "%d", nfree);
      strcat(points_info, number);
      XmString info_l3_label = 
       XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
      if(!RS.IsInfo) {
        XtVaSetValues(
         RS.info_l2,
         XmNlabelString, info_l2_label,
         NULL
        );
        XtVaSetValues(
         RS.info_l3,
         XmNlabelString, info_l3_label,
         NULL
        );
      }
    }
    return;
  }
  if(draw_data -> event -> xbutton.button == Button2) {
    Point Ret;
    if(RS.RPtr -> XgMouseRemove(P, &Ret)) {
      Undraw_point1(Ret);
      int npoin = RS.RPtr->XgGiveNpoin();
      char points_info[255];
      strcpy(points_info, RS.info_npoin_str);
      char number[20];
      sprintf(number, "%d", npoin);
      strcat(points_info, number);
      XmString info_l2_label = 
       XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
      int nfree = RS.RPtr->XgGiveNfree();
      strcpy(points_info, RS.info_nfree_str);
      sprintf(number, "%d", nfree);
      strcat(points_info, number);
      XmString info_l3_label = 
       XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
      if(!RS.IsInfo) {
        XtVaSetValues(
         RS.info_l2,
         XmNlabelString, info_l2_label,
         NULL
        );
        XtVaSetValues(
         RS.info_l3,
         XmNlabelString, info_l3_label,
         NULL
        );
      }
    }
    return;
  }
}

/*ARGSUSED*/
static void BlackboardResize(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  int size, pos;
  double division = 1;
  static int Count = 0;
  if(++Count <2) return;
  Dimension new_w, new_h;
  XtVaGetValues(       
   RS.blackboard,
   XmNheight, &new_h,   
   XmNwidth, &new_w,           
   NULL                       
  );
  if(RS.max_sliderx != 0)
    division = new_w/(double)RS.max_sliderx;
  RS.max_sliderx = new_w;
  RS.pos_sliderx = RS.pos_sliderx*division;
  pos = (int)(RS.pos_sliderx+0.5);
  if(pos < 0) pos = 0;
  size = (int)((RS.size_sliderx = new_w*RS.Measure) + 0.5);
  if(size > RS.max_sliderx) size = RS.max_sliderx;
  if(pos > new_w - size) pos = new_w - size;
  XtVaSetValues(
   RS.horiz_bar,
   XmNmaximum, new_w,
   XmNsliderSize, size,
   XmNvalue, pos,
   NULL
  );
//  Put(pos); Put("\n");
//  Put(new_w); Put("\n");
//  Put(size); Put("\n");
//  Put("\n");
  if(RS.max_slidery != 0)
    division = new_h/(double)RS.max_slidery;
  RS.max_slidery = new_h;
  RS.pos_slidery = RS.pos_slidery*division;
  pos = (int)(RS.pos_slidery+0.5);
  if(pos < 0) pos = 0;
  size = (int)((RS.size_slidery = new_h*RS.Measure) + 0.5);
  if(size > RS.max_slidery) size = RS.max_slidery;
  if(pos > new_h - size) pos = new_h - size;
  XtVaSetValues(
   RS.vert_bar,
   XmNmaximum, new_h,
   XmNsliderSize, size,
   XmNvalue, pos,
   NULL
  );
  Graphics_const_init(new_w, new_h);
  if(Count < 3) return;
  XClearWindow(
   XtDisplay(RS.blackboard),  
   XtWindow(RS.blackboard)    
  );    
  Redraw_blackboard();
}

/*ARGSUSED*/
static void BlackboardExpose(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  Redraw_blackboard();
}

/*
 *    WorkProcedure:
 */

/*ARGSUSED*/
Boolean Continue(XtPointer client_data) {
  if(RS.Grid_ready || RS.Grid_failed) return False;
  if(RS.Grid_working) {
    Point a, b, c;
    int Ready;
    if(RS.RPtr->XgNextTriangle(&a, &b, &c, &Ready)) {
      if(!RS.IsInfo) {
        char elem_info[50], number[10];
        strcpy(elem_info, RS.info_nelem_str);
        sprintf(number, "%ld", RS.RPtr->XgGiveNelem());
        strcat(elem_info, number);
        XmString info_l5_label = 
        XmStringCreate(elem_info, XmFONTLIST_DEFAULT_TAG);
        XtVaSetValues(
         RS.info_l5,
         XmNlabelString, info_l5_label,
         NULL
        );
      }
      if(Ready == 1) {
        RS.Grid_working = 0;
        RS.Grid_ready = 1;
        REMOVE_BDY_PTS_ACTIVE = 1;
        if(!RS.IsInfo) {
          XtVaSetValues(
           RS.info_l1,
           XmNlabelString, RS.info_l1_label3,
           NULL
          );
        }
        XtVaSetValues(
         RS.gridsave,
         XmNsensitive, True,
         XmNlabelString, RS.save_name,
         NULL
        );
      }
      Draw_line(a, b);
      Draw_line(a, c);
      Draw_line(b, c);
    }
    else {
      DeactivateWindow();
      RS.Grid_working = 0;
      RS.Grid_ready = 0;
      RS.Grid_failed = 1;
      if(!RS.IsInfo) { 
        XtVaSetValues(
         RS.info_l1,
         XmNlabelString, RS.info_l1_label4,
         NULL
        );
      }
      TrianglesErrorDialog();
    }   
  } 
  else {
    static int Bound_count = 0;
    if(!RS.RPtr->XgIsEmpty()) {
      for(int i=0; i<RS.SubLoopsNumber; i++) {
        Point P1, P2;
        RS.RPtr->XgNextShift(&P1, &P2);
        Undraw_point1(P1);
        Draw_point1(P2);
      }
    }
    if(++Bound_count > RS.BoundaryRedrawInterval) {
      Redraw_boundary(); 
      Bound_count = 0;
    }      
  }
  return(False); 
}

void MenuCancel(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtUnmanageChild(RS.menu);
  XtVaSetValues(
   RS.menu_button,
   XmNsensitive, True,
   NULL
  );
}

void BailOut(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  exit(0);
}

void InfoDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtUnmanageChild(RS.infodialog);
  XtVaSetValues(
   RS.info,
   XmNsensitive, True,
   NULL
  );
  RS.IsInfo = 1;
}

void Info(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtVaSetValues(
   RS.info,
   XmNsensitive, False,
   NULL
  );
  RS.IsInfo = 0;
  Widget ok;  

  int i = 0;
  Arg args[20];
  
  XtSetArg(args[i], XtNtitle, "Information"); i++;
  XtSetArg(args[i], XmNminWidth, 300); i++;
  XtSetArg(args[i], XmNminHeight, 230); i++;
  XtSetArg(args[i], XmNwidth, 300); i++;
  XtSetArg(args[i], XmNheight, 230); i++;
  XtSetArg(args[i], XmNresizable, False); i++;

  char msg[255];
  strcpy((char*)msg, "Infodialog");
  RS.infodialog = XmCreateFormDialog(
   RS.topLevel,
   msg,
   args,
   i
  );

  XmString info_l1_label;
  if(!RS.Grid_working && !RS.Grid_ready && !RS.Grid_failed) 
   info_l1_label = RS.info_l1_label1;
  if(RS.Grid_working) 
   info_l1_label = RS.info_l1_label2;
  if(RS.Grid_ready) 
   info_l1_label = RS.info_l1_label3;
  if(RS.Grid_failed) 
   info_l1_label = RS.info_l1_label4;

  int npoin = RS.RPtr->XgGiveNpoin();
  char points_info[255];
  strcpy(points_info, RS.info_npoin_str);
  char number[20];
  sprintf(number, "%d", npoin);
  strcat(points_info, number);
  XmString info_l2_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);

  int nfree = RS.RPtr->XgGiveNfree();
  strcpy(points_info, RS.info_nfree_str);
  sprintf(number, "%d", nfree);
  strcat(points_info, number);
  XmString info_l3_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);

  strcpy(points_info, "Boundary points: ");
  sprintf(number, "%d", npoin - nfree);
  strcat(points_info, number);
  XmString info_l4_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);

  points_info[0] = '\0';
  if(RS.RPtr->XgGiveNelem() != 0) {
    strcpy(points_info, RS.info_nelem_str);
    sprintf(number, "%ld", RS.RPtr->XgGiveNelem());
    strcat(points_info, number);
  }
  else {
    strcpy(points_info, RS.info_nelem_str);
    strcat(points_info, "0");
  }

  XmString info_l5_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);

  strcpy(points_info, RS.info_timestep_str);
  sprintf(number, "%g", RS.RPtr->XgGiveTimestep());
  strcat(points_info, number);
  XmString info_l6_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);

  RS.info_l1 = XtVaCreateManagedWidget(
   "info_l1",
   xmLabelWidgetClass,
   RS.infodialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 8,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 18,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, info_l1_label,
   NULL
  ); 
  XtManageChild(RS.info_l1);

  RS.info_l2 = XtVaCreateManagedWidget(
   "info_l2",
   xmLabelWidgetClass,
   RS.infodialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 22,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 32,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, info_l2_label,
   NULL
  ); 
  XtManageChild(RS.info_l2);

  RS.info_l3 = XtVaCreateManagedWidget(
   "info_l3",
   xmLabelWidgetClass,
   RS.infodialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 32,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 42,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, info_l3_label,
   NULL
  ); 
  XtManageChild(RS.info_l3);

  RS.info_l4 = XtVaCreateManagedWidget(
   "info_l4",
   xmLabelWidgetClass,
   RS.infodialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 42,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 52,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, info_l4_label,
   NULL
  ); 
  XtManageChild(RS.info_l4);

  RS.info_l5 = XtVaCreateManagedWidget(
   "info_l5",
   xmLabelWidgetClass,
   RS.infodialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 52,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 62,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, info_l5_label,
   NULL
  ); 
  XtManageChild(RS.info_l5);

  RS.info_l6 = XtVaCreateManagedWidget(
   "info_l6",
   xmLabelWidgetClass,
   RS.infodialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 62,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 72,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, info_l6_label,
   NULL
  ); 
  XtManageChild(RS.info_l6);

  ok = XtVaCreateManagedWidget(
   "Ok",
   xmPushButtonWidgetClass,
   RS.infodialog,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 80,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 94,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 38,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 62,
   XmNlabelString, RS.ok_label,
   NULL
  );  
 
  XtManageChild(ok);

  XtAddCallback(ok, 
   XmNactivateCallback, InfoDialogOK, NULL);

  XtManageChild(RS.infodialog);   
}


void Menu(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  int i;
  Arg args[20];

  i=0;  
  XtSetArg(args[i], XtNtitle, "Menu"); i++;
  XtSetArg(args[i], XmNminWidth, 90); i++;
  XtSetArg(args[i], XmNminHeight, 190); i++;
  XtSetArg(args[i], XmNwidth, 140); i++;
  XtSetArg(args[i], XmNheight, 270); i++;
  XtSetArg(args[i], XmNautoUnmanage, False); i++;
  XtSetArg(args[i], XmNfractionBase, 1000); i++;

  char msg[255];
  strcpy((char*)msg, "Menu");
  RS.menu = XmCreateFormDialog(
   RS.menu_button,
   msg,
   args,
   i
  );

  RS.about = XtVaCreateManagedWidget(
   "about",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.about_name,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 111,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNsensitive, RS.IsAbout ? True : False,
   NULL
  );  

  RS.timeinc = XtVaCreateManagedWidget(
   "timeinc",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.timeinc_name,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 111,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 222,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   NULL
  );  

  RS.timedec = XtVaCreateManagedWidget(
   "timedec",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.timedec_name,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 222,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 333,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   NULL
  ); 

  RS.zoominc = XtVaCreateManagedWidget(
   "zoominc",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.zoominc_name,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 333,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 444,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   NULL
  );  

  RS.zoomdec = XtVaCreateManagedWidget(
   "zoomdec",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.zoomdec_name,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 444,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 555,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   NULL
  );  

  RS.savepoints = XtVaCreateManagedWidget(
   "savepoints",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.savepoints_name,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 555,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 666,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   NULL
  ); 

  RS.info = XtVaCreateManagedWidget(
   "Info",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.info_name,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 666,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 777,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNsensitive, RS.IsInfo ? True : False,
   NULL
  ); 

  RS.menucancel = XtVaCreateManagedWidget(
   "menucancel",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.cancel_label,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 777,
   XmNbottomAttachment, XmATTACH_POSITION,
   XmNbottomPosition, 889,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 100,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 900,
   NULL
  ); 

  RS.bailout = XtVaCreateManagedWidget(
   "bailout",
   xmPushButtonWidgetClass,
   RS.menu,
   XmNlabelString, RS.bailout_name,
   XmNtopAttachment, XmATTACH_POSITION,
   XmNtopPosition, 889,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   NULL
  ); 

  XtManageChild(RS.about);
  XtManageChild(RS.timeinc);
  XtManageChild(RS.timedec);
  XtManageChild(RS.zoominc);
  XtManageChild(RS.zoomdec);
  XtManageChild(RS.savepoints);
  XtManageChild(RS.bailout);
  XtManageChild(RS.menucancel);
  XtManageChild(RS.info);
  XtManageChild(RS.menu);

  XtAddCallback(RS.about, 
   XmNactivateCallback, About, NULL);
  XtAddCallback(RS.timeinc, 
   XmNactivateCallback, TimeInc, NULL);
  XtAddCallback(RS.timedec, 
   XmNactivateCallback, TimeDec, NULL);
  XtAddCallback(RS.zoominc, 
   XmNactivateCallback, ZoomInc, NULL);
  XtAddCallback(RS.zoomdec, 
   XmNactivateCallback, ZoomDec, NULL);
  XtAddCallback(RS.savepoints, 
   XmNactivateCallback, SavePoints, NULL);
  XtAddCallback(RS.bailout, 
   XmNactivateCallback, BailOut, NULL);
  XtAddCallback(RS.menucancel, 
   XmNactivateCallback, MenuCancel, NULL);
  XtAddCallback(RS.info, 
   XmNactivateCallback, Info, NULL);

  XtVaSetValues(
   RS.menu_button,
   XmNsensitive, False,
   NULL
  );
}

/*
 *    Graphic setup:     
 */

static void SetUpThings() {
  XGCValues values;

  RS.black_color =
   XtScreen(RS.blackboard) -> black_pixel;
  RS.white_color =
   XtScreen(RS.blackboard) -> white_pixel;     

  //drawing white on black hardcoded:
  values.foreground = RS.white_color;
  values.background = RS.black_color;
  RS.draw_gc = XCreateGC(
    XtDisplay(RS.blackboard), 
    XtWindow(RS.blackboard),
    GCForeground | GCBackground, 
    &values
  );
  values.foreground = RS.black_color;
  values.background = RS.white_color;
  RS.undraw_gc = XCreateGC(
    XtDisplay(RS.blackboard), 
    XtWindow(RS.blackboard),
    GCForeground | GCBackground, 
    &values
  );

  //blackboards background hardcoded:
  XtVaSetValues(
   RS.blackboard,
   XmNbackground, RS.black_color,
   NULL
  );
}                                

void RSGetConfiguration(Xgen *User) {
  //at first: adding User to RS:
  RS.RPtr = User;

  //button names hardcoded:
  char msg[255];
  strcpy((char*)msg, "Menu");
  RS.menu_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "About");
  RS.about_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Timestep inc");
  RS.timeinc_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Timestep dec");
  RS.timedec_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Zoom in");
  RS.zoominc_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Zoom out");
  RS.zoomdec_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Save points");
  RS.savepoints_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Init");
  RS.init_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Remove");
  RS.forget_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Mesh");
  RS.grid_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Save");
  RS.save_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Quit");
  RS.quit_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Yes");
  RS.yes_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "No");
  RS.no_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Ok");
  RS.ok_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Cancel");
  RS.cancel_label = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Exit");
  RS.bailout_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);
  strcpy((char*)msg, "Info");
  RS.info_name = 
   XmStringCreate(msg, XmFONTLIST_DEFAULT_TAG);

  //initialization of info-board:
  strcpy((char*)msg, "Nothing is happening.");
  RS.info_l1_label1 = XmStringCreate(
    msg, XmFONTLIST_DEFAULT_TAG); 
  strcpy((char*)msg, "Mesh is being created.");
  RS.info_l1_label2 = XmStringCreate(
    msg, XmFONTLIST_DEFAULT_TAG); 
  strcpy((char*)msg, "Mesh is ready.");
  RS.info_l1_label3 = XmStringCreate(
    msg, XmFONTLIST_DEFAULT_TAG); 
  strcpy((char*)msg, "The algorithm failed.");
  RS.info_l1_label4 = XmStringCreate(
    msg, XmFONTLIST_DEFAULT_TAG);
  strcpy(RS.info_npoin_str, "Number of points: ");
  strcpy(RS.info_nfree_str, "Interior points: ");
  strcpy(RS.info_nelem_str, "Number of elements: ");
  strcpy(RS.info_timestep_str, "Timestep: ");

  //initialization of area extrems:
  Point min, max;
  User -> XgGiveLimits(&min, &max);
  RS.AreaSize = max - min;
  RS.min = min;
  RS.max = max;

  //some more app-initializations:
  RS.Grid_working = 
  RS.Grid_ready = RS.Grid_failed = 0;  
  RS.Ratio = Point(0, 0);
  RS.Init_pos_x = 0;
  RS.Init_pos_y = 0;
  RS.MAX_X = 0;
  RS.MAX_Y = 0;
  RS.Measure = 1.0;
  RS.pos_sliderx = RS.pos_slidery = 0;
  RS.max_sliderx = RS.max_slidery = 0;
  RS.size_sliderx = RS.size_slidery = 0; 
  RS.IsZoominc = RS.IsInfo = RS.IsAbout = 1;

  //case error while reading:
  strcpy(RS.PointPath, "points/");
  strcpy(RS.GridPath, "meshes/");
  RS.SubLoopsNumber = 10;
  RS.BoundaryRedrawInterval = 100;
  RS.InitSizeLimit = 600;
  RS.ZoomConst = 0.95;
  double timestep_const = 0.8;

  FILE *f= fopen(".XgenConfig", "rb");
  if(f == NULL) return;

  if(!Get(f, RS.PointPath)) return;
  if(!Get(f, RS.GridPath)) return;
  if(!Get(f, &RS.SubLoopsNumber)) return;
  if(RS.SubLoopsNumber == 0) 
   XgError("Invalid SubLoopsNumber."); 
  if(!Get(f, &RS.InitSizeLimit)) return;
  if(RS.InitSizeLimit == 0) 
   XgError("Invalid InitSizeLimit."); 
  if(!Get(f, &RS.BoundaryRedrawInterval)) return;
  if(RS.BoundaryRedrawInterval == 0) 
   XgError("Invalid BoundaryRedrawInterval."); 
  if(!Get(f, &timestep_const)) return;
  RS.RPtr->XgSetTimestepConstant(timestep_const);
  if(!Get(f, &RS.ZoomConst)) return;
  if(RS.ZoomConst <= 0 || RS.ZoomConst >= 1) 
   XgError("Invalid ZoomConst."); 

  fclose(f);
}

/* 
 *    Main program function:
 */

void XgMainLoop(Xgen *User, int argc, char *argv[]) {
  RSGetConfiguration(User);

  //getting the ApplicationsName:
  char AppName[255];
  strcpy(AppName, argv[0]);
  for(int k=1; k<argc; k++) {
    strcat(AppName, " ");
    strcat(AppName, argv[k]);
  }

  //getting optimal size of the xgen-Window:
  Get_optimal_size( 
   &(RS.Optimal_width),
   &(RS.Optimal_height)  
  );

  //installing the xgen-Window:
  XtAppContext app_context;

  RS.topLevel = XtVaAppInitialize(
   &app_context,
   "XGen",
   NULL, 0,
   &argc, argv,
   NULL,
   XtNtitle, AppName,
   NULL
  );             

  RS.form = XtVaCreateManagedWidget(
   "form",
   xmFormWidgetClass,
   RS.topLevel,
   XmNfractionBase, 1000,
   XmNheight, RS.Optimal_height,
   XmNwidth, RS.Optimal_width,
   NULL    
  );    

  RS.form2 = XtVaCreateManagedWidget(
   "form2",
   xmFormWidgetClass,
   RS.form,
   XmNtopAttachment, XmATTACH_NONE,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNfractionBase, 1000,
   NULL 
  ); 

  RS.form1 = XtVaCreateManagedWidget(
   "form1",
   xmFormWidgetClass,
   RS.form,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_WIDGET,
   XmNbottomWidget, RS.form2,
   XmNbottomOffset, 1,
   XmNleftAttachment, XmATTACH_FORM,
   XmNrightAttachment, XmATTACH_FORM,
   XmNfractionBase, 1000,
   NULL 
  ); 

  Arg args[20]; int i; 

  i=0;
  XtSetArg(args[i], XmNtopAttachment, XmATTACH_POSITION); i++; 
  XtSetArg(args[i], XmNtopPosition, 4); i++; 
  XtSetArg(args[i], XmNbottomAttachment, XmATTACH_FORM); i++; 
  XtSetArg(args[i], XmNleftAttachment, XmATTACH_POSITION); i++; 
  XtSetArg(args[i], XmNleftPosition, 4); i++; 
  XtSetArg(args[i], XmNrightAttachment, XmATTACH_FORM); i++; 
  XtSetArg(args[i], XmNscrollingPolicy, XmAPPLICATION_DEFINED); i++;   
  char msg[255];
  strcpy((char*)msg, "scrolled");
  RS.scrolled = XmCreateScrolledWindow(
   RS.form1,
   msg,
   args,
   i
  );

  strcpy((char*)msg, "blackboard");
  RS.blackboard = XtVaCreateManagedWidget(
   msg,
   xmDrawingAreaWidgetClass,
   RS.scrolled,
   NULL    
  );

  i=0; 
  strcpy((char*)msg, "horiz_bar");
  XtSetArg(args[i], XmNorientation, XmHORIZONTAL); i++;   
  RS.horiz_bar = XmCreateScrollBar(
   RS.scrolled,
   msg,
   args,
   i
  );

  i=0;
  strcpy((char*)msg, "vert_bar");
  XtSetArg(args[i], XmNorientation, XmVERTICAL); i++;   
  RS.vert_bar = XmCreateScrollBar(
   RS.scrolled,
   msg,
   args,
   i
  );

  XmScrolledWindowSetAreas(
   RS.scrolled,
   RS.horiz_bar,
   RS.vert_bar,
   RS.blackboard
  );

  XtManageChild(RS.scrolled);
  XtManageChild(RS.blackboard);
  XtManageChild(RS.horiz_bar);
  XtManageChild(RS.vert_bar);

  RS.menu_button = XtVaCreateManagedWidget(
   "menu_button",
   xmPushButtonWidgetClass,
   RS.form2,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 3,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 238,
   XmNlabelString, RS.menu_name,
   NULL
  );  
  XtManageChild(RS.menu_button);

  RS.initforget = XtVaCreateManagedWidget(
   "initforget",
   xmPushButtonWidgetClass,
   RS.form2,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 256,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 491,
   XmNlabelString, RS.init_name,
   NULL
  ); 
  
  RS.gridsave = XtVaCreateManagedWidget(
   "gridsave",
   xmPushButtonWidgetClass,
   RS.form2,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 509,
   XmNrightAttachment, XmATTACH_POSITION,
   XmNrightPosition, 746,
   XmNlabelString, RS.grid_name,
   NULL
  );

  RS.quit = XtVaCreateManagedWidget(
   "Quit",
   xmPushButtonWidgetClass,
   RS.form2,
   XmNtopAttachment, XmATTACH_FORM,
   XmNbottomAttachment, XmATTACH_FORM,
   XmNleftAttachment, XmATTACH_POSITION,
   XmNleftPosition, 764,
   XmNrightAttachment, XmATTACH_FORM,
   XmNlabelString, RS.quit_name,
   NULL
  );

  //adding callback functions:
  XtAddCallback(RS.menu_button, 
   XmNactivateCallback, Menu, NULL);
  XtAddCallback(RS.initforget,
    XmNactivateCallback, InitForget, NULL);
  XtAddCallback(RS.gridsave, 
   XmNactivateCallback, GridSave, NULL);
  XtAddCallback(RS.quit, 
   XmNactivateCallback, Quit, NULL);
  XtAddCallback(RS.blackboard, 
   XmNinputCallback, BlackboardInput, NULL);
  XtAddCallback(RS.blackboard, 
   XmNresizeCallback, BlackboardResize, NULL);
  XtAddCallback(RS.blackboard, 
   XmNexposeCallback, BlackboardExpose, NULL);
  XtAddCallback(RS.horiz_bar, 
   XmNdragCallback, ShiftX, NULL);
  XtAddCallback(RS.vert_bar, 
   XmNdragCallback, ShiftY, NULL);
  XtAddCallback(RS.horiz_bar, 
   XmNvalueChangedCallback, ShiftX, NULL);
  XtAddCallback(RS.vert_bar, 
   XmNvalueChangedCallback, ShiftY, NULL);

  //adding work-procedure:
  XtAppAddWorkProc(app_context, Continue, NULL);

  //starting the XGen-Window:
  XtRealizeWidget(RS.topLevel); 

  //graphics setup:
  SetUpThings();

  //running the XGen-Window:
  XtAppMainLoop(app_context);
}

















