/**
 **  A simple 2D interactive mesh generator Xgen 
 **  (c) Pavel Solin, 1998 - 2011
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

bool REMOVE_BDY_PTS_ACTIVE = true;

/*
 *   Error message:
 */

void XgError(char *what) {
  fprintf(stderr, "Error: %s\n", what);
  exit(0);
}

void XgError(const char *what) {
  fprintf(stderr, "Error: %s\n", what);
  exit(0);
}

/**
 **   struct BoundaryPairsList:
 **/

void BoundaryPairsList::ChangeLast(int a, int b) 
{
  {Last->A = a; Last->B = b;}
}

bool BoundaryPairsList::GetNext(int &a, int &b) 
{
  if(Ptr != NULL) {
    a = Ptr->A;
    b = Ptr->B;
    Ptr = Ptr->next;
    return true;
  }
  else return false;
}

void BoundaryPairsList::Delete() {
  Ptr = First;
  BoundaryPair *p;
  while(Ptr != NULL) {
    p = Ptr;
    Ptr = Ptr -> next;
    delete p;
  }
  First = Last = NULL;
}

bool BoundaryPairsList::Contains(int C) {
  BoundaryPair *P = First;
  while(P != NULL) {
    if((P->A == C) || (P->B == C)) return true;
    P = P->next;
  }
  return false;
}

void BoundaryPairsList::Add(int a, int b) 
{
  if (Last == NULL) {
    First = Last = new BoundaryPair(a, b);
  }
  else {
    Last->next = new BoundaryPair(a, b);
    Last = Last->next;
  }
}

void BoundaryPairsList::DeleteLast() {
  BoundaryPair *p = First;
  if (p->next == NULL) {
    delete p;
    First = Last = NULL;
  }
  else {
    if (p->next->next == NULL) {
      delete Last;
      Last = First;
    }
    else {
      while (p->next->next != NULL) p = p->next;
      Last = p;
      delete p->next;
    }
    Last->next = NULL;
  }
}

bool BoundaryPairsList::Delete(int a, int b) 
{
  BoundaryPair *p = First, *l = NULL;
  while (p != NULL) {
    if ((p->A == a && p->B == b) || (p->A == b && p->B == a)) {
      if (l == NULL) {
        First = p->next;
	if (p->next == NULL) Last = NULL;
      }
      else {
	l->next = p->next;
	if (p->next == NULL) {
	  Last = l;
	  Last->next = NULL;
        }
      }
      delete p;
      return true;
    }
    else {
      l = p;
      p = p->next;
    }
  }
  return false;
}

/**
 **   struct BoundaryLinesList:
 **/

bool BoundaryLinesList::GetNext(Point &a, Point &b) 
{
  if(Ptr != NULL) {
    a = Ptr->A;
    b = Ptr->B;
    Ptr = Ptr->next;
    return true;
  }
  else return false;
}

void BoundaryLinesList::Delete() {
  Ptr = First;
  BoundaryLine *bl;
  while(Ptr != NULL) {
    bl = Ptr;
    Ptr = Ptr -> next;
    delete bl;
  }
  First = Last = NULL;
}

void BoundaryLinesList::AddStraightLine(Point a, Point b) 
{
  if (Last == NULL) {
    First = Last = new BoundaryLine(a, b);
  }
  else {
    Last->next = new BoundaryLine(a, b);
    Last = Last->next;
  }
}

// Calculates center of circle with central angle alpha, passing through 
// points a, b.
void Calculate_circle_center(Point a, Point b, double alpha, Point &S) {
  Point c = (a + b) / 2.;
  Point ab = b - a;  
  double length_ab = ab.abs();
  double R = (length_ab / 2) / sin(alpha * M_PI / 180 / 2.);        // circle radius
  double r0 = R * cos(alpha * M_PI / 180 / 2.);                     // distance between C and circle center
  Point u = Point(-ab.y, ab.x) / length_ab;            // unit vector from C to circle center S
  S = c + u * r0;                                // circle center
  //printf("R = %g\n", R);
  //printf("r0 = %g\n", r0);
  //printf("C = %g %g\n", c.x, c.y);
  //printf("U = %g %g\n", u.x, u.y);
  //printf("S = %g %g\n", S.x, S.y);
}

// Here a, b are the arc end points and theta is an angle in degrees. With 
// theta = 0 one obtains a. With theta = alpha one obtains b.
void Calculate_point_on_arc(Point a, Point b, double alpha, double theta, Point &a_next) 
{
  // Calculate circle center.
  Point S;
  Calculate_circle_center(a, b, alpha, S);

  // Rotation matrix.
  double m_11 =  cos(theta*M_PI/180.); 
  double m_12 = -sin(theta*M_PI/180.); 
  double m_21 =  sin(theta*M_PI/180.); 
  double m_22 =  cos(theta*M_PI/180.); 
  //printf("Matrix: %g %g\n", m_11, m_12);
  //printf("        %g %g\n", m_21, m_22);

  Point Sa = a - S;           // vector from S to a
  //printf("Vector SA: %g %g\n", Sa.x, Sa.y);
  Point rotated_vec;
  rotated_vec.x = m_11 * Sa.x + m_12 * Sa.y;
  rotated_vec.y = m_21 * Sa.x + m_22 * Sa.y;
  a_next = S + rotated_vec; 
  //printf("Point A_next: %g %g\n", a_next.x, a_next.y);
}

void BoundaryLinesList::AddCircularArc(Point a, Point b, double alpha, double H) 
{
  // Number of straight lines to approximate a circle.
  int n_min = 1;
  double angle_increment = 5;   // decreasing this will produce finer plots;
  int n = (int) (alpha / angle_increment + 0.5); 
  if (n < n_min) n = n_min;

  // Angle increment correponding to one straight line.
  double theta = alpha / n;

  // Check the length of the elementary line approximating the arc since
  // we need to decrease the abgle if the length is less than this->H.
  // This is because of boundary point removal.
  Point a2;
  Calculate_point_on_arc(a, b, alpha, theta, a2); 
  Point Da2 = a2 - a;
  double h = Da2.abs();
  if (h > H) theta *= H / (2*h);

  // Add elementary lines approximating the arc with angle increment theta.
  for (int i=0; i < n; i++) {
    Point A1, A2;
    Calculate_point_on_arc(a, b, alpha, i * theta, A1);     
    Calculate_point_on_arc(a, b, alpha, (i+1) * theta, A2);     
    AddStraightLine(A1, A2); // adding the elementary line
    //printf("Adding elementary line (%g, %g), (%g, %g)\n", A1.x, A1.y, A2.x, A2.y);
  }
}

void BoundaryLinesList::DeleteLast() {
  BoundaryLine *p = First;
  if (p->next == NULL) {
    delete p;
    First = Last = NULL;
  }
  else {
    if (p->next->next == NULL) {
      delete Last;
      Last = First;
    }
    else {
      while (p->next->next != NULL) p = p->next;
      Last = p;
      delete p->next;
    }
    Last->next = NULL;
  }
}

/**
 **   struct OutputBoundaryEdgeList:
 **/

bool OutputBoundaryEdgeList::GetNext(BoundaryEdge &BE) {
  if(Ptr != NULL) {
    BE = *Ptr;
    Ptr = Ptr->next;
    return true;
  }
  else return false;
}

void OutputBoundaryEdgeList::Add(int a, int b, int component_index, int marker, double alpha) 
{
  if (Last == NULL) {
    First = Last = new BoundaryEdge(a, b, component_index, marker, alpha);
  }
  else {
    Last->next = new BoundaryEdge(a, b, component_index, marker, alpha);
    Last = Last->next;
  }
}

void OutputBoundaryEdgeList::DeleteLast() {
  BoundaryEdge *p = First;
  if (p->next == NULL) {
    delete p;
    First = Last = NULL;
  }
  else {
    if (p->next->next == NULL) {
      delete Last;
      Last = First;
    }
    else {
      while (p->next->next != NULL) p = p->next;
      Last = p;
      delete p->next;
    }
    Last->next = NULL;
  }
}

void OutputBoundaryEdgeList::Delete() {
  BoundaryEdge *p = Ptr = First;
  while(p != NULL) {
    Ptr = Ptr->next;
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
    Last->next = new ElemBox(a, b, c);
    Last = Last->next;
    if(Last == NULL) XgError("ElemList::Add(): Not enough memory."); 
  }
}

bool ElemList::GetNext(int &a, int &b, int &c) 
{
  if(Ptr != NULL) {
    a = Ptr->E.n1;
    b = Ptr->E.n2;
    c = Ptr->E.n3;
    Ptr = Ptr->next;
    return true;
  }
  else return false;
}

void ElemList::Delete() {
  ElemBox *p;
  Ptr = First;
  while(Ptr != NULL) { 
    p = Ptr;
    Ptr = Ptr->next;
    delete p;
  }
  First = Last = NULL; 
}

/**
 **   struct PointList:
 **/

bool PointList::GetNext(int i, double &pos_x, double &pos_y
) {
  PointBox *P = First;
  for(int j=0; j<i; j++) {
    if(P == NULL) return false;
    P = P->next;
  }
  pos_x = P->P.x;
  pos_y = P->P.y;
  return true;
}

bool PointList::GetNext(Point &T) {
  if(Ptr != NULL) {
    T.x = Ptr->P.x;
    T.y = Ptr->P.y;
    Ptr = Ptr->next; 
    return true;
  }
  else return false;
}

void PointList::Add(Point P) {
  Add(P.x, P.y);
}

void PointList::Add(double pos_x, double pos_y) {
  if (Last == NULL) {
    First = Last = new PointBox(pos_x, pos_y);
  }
  else {
    Last->next = new PointBox(pos_x, pos_y);
    Last = Last->next;
  }
}

Point PointList::Delete() {
  Point P(0, 0);

  if(First == NULL) return P;
  else {
    PointBox *p = First;
    P = p->P;
    First = First->next;
    delete p;
    return P;
  }
}

/**
 **   struct BoundarySegment:
 **/

double BoundarySegment::CalculateLength() 
{
  if (fabs(alpha) < 1e-3) {
    Point straight_line = P_end - P_start;
    return straight_line.abs();
  }
  else {
    Point S;
    Calculate_circle_center(P_start, P_end, alpha, S);
    Point radius_line = P_start - S;
    double R = radius_line.abs();
    double perimeter = 2*M_PI*R;
    return alpha / 360 * perimeter;   // this is exact
  }
}

/**
 **   struct Boundary:
 **/

void BoundaryType::CreateNewBoundaryComponent() {
  if(Last_bdy_component == NULL) {
    First_bdy_component = Last_bdy_component = new BoundaryComponent();
    if(First_bdy_component == NULL) 
      XgError("BoundaryType::CreateNewBoundaryComponent(): Not enough memory.");
  }
  else {
    BoundaryComponent *help = Last_bdy_component;
    Last_bdy_component = new BoundaryComponent();
    help -> next = Last_bdy_component;
    if(help -> next == NULL) 
      XgError("BoundaryType::CreateNewBoundaryComponent(): Not enough memory.");
  }
}

void BoundaryType::AddBoundarySegment(int marker, double x, double y, int subdiv, double alpha) 
{
  if(Last_bdy_component == NULL) {
    First_bdy_component = Last_bdy_component = new BoundaryComponent();
    if(First_bdy_component == NULL) 
      XgError("BoundaryType::AddBoundarySegment(): Not enough memory.");    
  }
  if(Last_bdy_component -> Last_segment == NULL) {
    Last_bdy_component -> Last_segment = 
    Last_bdy_component -> First_segment = new BoundarySegment(marker, x, y, subdiv, alpha);
    if(Last_bdy_component -> Last_segment == NULL) 
      XgError("BoundaryType::AddBoundarySegment(): Not enough memory.");    
  } 
  else {
    BoundarySegment *help = Last_bdy_component -> Last_segment;
    Last_bdy_component -> Last_segment = new BoundarySegment(marker, x, y, subdiv, alpha);
    help -> P_end = Last_bdy_component -> Last_segment -> P_start;
    Last_bdy_component -> Last_segment -> P_end = Last_bdy_component -> First_segment -> P_start;
    help -> next = Last_bdy_component -> Last_segment;
    if(help -> next == NULL) 
      XgError("BoundaryType::AddBoundarySegment(): Not enough memory."); 
  }
}

BoundaryPairsList* BoundaryType::CreateBoundaryPairsList() {
  BoundaryPairsList* bpl = new BoundaryPairsList;
  int Count = 0;

  BoundaryComponent *Bdy_component_pointer = First_bdy_component;
  while(Bdy_component_pointer != NULL) {
    BoundarySegment *Segment_pointer = Bdy_component_pointer -> First_segment;  
    int first_point_of_component = Count;
    while(Segment_pointer != NULL) {
      int subdiv = Segment_pointer -> subdiv;
      int first_point_of_segment = Count;
      for(int i=0; i<subdiv; i++) {
        //printf("Adding boundary pair %d %d\n", Count, Count + 1);
        bpl->Add(Count, Count + 1);
        Count++;
      }
      Segment_pointer = Segment_pointer -> next;
    }
    //printf("Deleting last\n");
    bpl->DeleteLast();
    //printf("Adding boundary pair %d %d\n", Count-1, first_point_of_component);
    bpl->Add(Count-1, first_point_of_component);
    Bdy_component_pointer = Bdy_component_pointer -> next;
  }
  return bpl;
}       

BoundaryLinesList* BoundaryType::CreateBoundaryLinesList(double H) {
  BoundaryLinesList* bll = new BoundaryLinesList;

  BoundaryComponent *Bdy_component_pointer = First_bdy_component;
  while(Bdy_component_pointer != NULL) {
    BoundarySegment *Segment_pointer = Bdy_component_pointer -> First_segment;  
    Point first_point_of_component = Segment_pointer -> P_start;
    while(Segment_pointer != NULL) {
      Point first_point_of_segment = Segment_pointer -> P_start;
      Point last_point_of_segment;
      if (Segment_pointer != Bdy_component_pointer -> Last_segment) 
        last_point_of_segment = Segment_pointer -> next -> P_start;
      else last_point_of_segment = first_point_of_component;
      double alpha = Segment_pointer -> alpha;
      if (fabs(alpha) < 1e-3) {
        bll->AddStraightLine(first_point_of_segment, last_point_of_segment);
      }
      else {
        bll->AddCircularArc(first_point_of_segment, last_point_of_segment, alpha, H);
      }
      Segment_pointer = Segment_pointer -> next;
    }
    Bdy_component_pointer = Bdy_component_pointer -> next;
  }
  return bll;
}       

PointList* BoundaryType::CreateBoundaryPointList() {
  PointList* bpl = new PointList;
  Point first_point_of_component, first_point, next_point, P;
  int subdiv;
  double angle;

  BoundaryComponent *Bdy_component_pointer = First_bdy_component;
  while(Bdy_component_pointer != NULL) {
    BoundarySegment *Segment_pointer = Bdy_component_pointer -> First_segment;  
    first_point_of_component = first_point = Segment_pointer->P_start;
    while(Segment_pointer->next != NULL) {
      subdiv = Segment_pointer -> subdiv;
      angle = Segment_pointer -> alpha;
      Segment_pointer = Segment_pointer->next;
      next_point = Segment_pointer->P_start;
      if (fabs(angle) < 1e-3) {
        P = (next_point - first_point)/subdiv;
        for(int i=0; i<subdiv; i++) {
          Point PP = first_point + P*i;
          bpl->Add(PP);
        }
      }
      else {
        for(int i=0; i<subdiv; i++) {
          Calculate_point_on_arc(first_point, next_point, angle, i * angle / subdiv, P); 
          bpl->Add(P);
        }
      }
      first_point = next_point;
    }
    subdiv = Segment_pointer -> subdiv;
    angle = Segment_pointer -> alpha;
    next_point = first_point_of_component;
    if (fabs(angle) < 1e-3) {
      P = (next_point - first_point)/subdiv;
      for(int i=0; i<subdiv; i++) {
        Point PP = first_point + P*i;
        bpl->Add(PP);
      }
    }
    else {
      for(int i=0; i<subdiv; i++) {
        Calculate_point_on_arc(first_point, next_point, angle, i * angle / subdiv, P); 
        bpl->Add(P);
      }
    }
    Bdy_component_pointer = Bdy_component_pointer -> next;
  }
  return bpl;
}

void BoundaryType::CalculateLength() 
{
  this->Length = 0;
  BoundaryComponent *component = First_bdy_component;
  while(component != NULL) {
    BoundarySegment *segment = component -> First_segment;  
    while(segment != NULL) {
      this->Length += segment->CalculateLength();
      segment = segment->next;
    }
    component = component -> next;
  }
}

OutputBoundaryEdgeList* BoundaryType::CreateOutputBoundaryEdgeList() {
  OutputBoundaryEdgeList* bel = new OutputBoundaryEdgeList;
  int Count = 0, First_point_of_component = 0;
  int subdiv = 0, component_index = 0;
  int last_marker, last_subdiv;
  double last_alpha;

  BoundaryComponent *Bdy_component_pointer = First_bdy_component;
  while(Bdy_component_pointer != NULL) {
    BoundarySegment *Segment_pointer = Bdy_component_pointer -> First_segment;
    First_point_of_component = Count;
    while(Segment_pointer != NULL) {
      subdiv = Segment_pointer -> subdiv;
      // Adding elementary boundary edges.
      for(int i=0; i<subdiv; i++) {
        bel->Add(
            Count, Count + 1, 
            component_index,
            Segment_pointer -> marker,
            (Segment_pointer -> alpha) / subdiv
           );
        Count++;
      }
      last_marker = Segment_pointer -> marker;
      last_alpha = Segment_pointer -> alpha;
      last_subdiv = Segment_pointer -> subdiv;
      Segment_pointer = Segment_pointer -> next;
    }
    bel->DeleteLast();
    bel->Add(
        Count - 1, First_point_of_component,
        component_index,
        last_marker,
        last_alpha / last_subdiv
       );
    Bdy_component_pointer = Bdy_component_pointer -> next;
    component_index++;
  }
  return bel;
}       

BoundaryType::~BoundaryType() {
  BoundaryComponent *component = First_bdy_component;
  while (component != NULL) {
    BoundarySegment *segment = component -> First_segment;
    while(segment != NULL) {
      BoundarySegment *help_segment = segment;
      segment = segment -> next;
      delete help_segment;
    }
    BoundaryComponent *help_component = component;
    component = component -> next;
    delete help_component;
  }
  First_bdy_component = Last_bdy_component = NULL;
}

/**
 **   class Xgen: public methods
 **/

Xgen::Xgen(bool nogui, int steps_to_take, bool overlay) {
  // Printing logo.  
  printf("\n");
  printf("-----------------------------------------\n");
  printf("  Xgen, a 2D interactive mesh generator  \n");
  printf("            Linux version 6.0            \n");
  printf("        Last modified Dec 5, 2003        \n");
  printf("    Revised and updated November 2010    \n");
  printf("   Copyright (C) 1995-2011 Pavel Solin   \n");
  printf("          e-mail: solin@unr.edu          \n");
  printf("    Distributed under the BSD license    \n");
  printf("-----------------------------------------\n");

  this->nogui = nogui;
  this->steps_to_take = steps_to_take;
  this->overlay = overlay;
  if (nogui) { 
    printf("Xgen will operate in batch mode.\n");
    printf("%d relaxation steps will be performed.\n", steps_to_take);
  }
  else printf("Xgen will operate in interactive mode.\n");
  if (overlay) printf("Xgen will use an equidistant overlay pattern for initial point positions.\n");
  else printf("Random initial point distribution will be used.\n");
}

char* Xgen::XgGiveName() {
  char *c = (char*)malloc(255);
  if(c == NULL) XgError("Xgen::XgGiveName(): Not enough memory.");
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

long Xgen::XgGiveInteriorPtsNum() {
  return InteriorPtsNum;
}

long Xgen::XgGiveBoundaryPtsNum() {
  return BoundaryPtsNum;
}

long Xgen::XgGiveNelem() {
  return Nelem;
}

void Xgen::XgGiveLimits(Point &min, Point &max) {
  min.x = Xmin;
  min.y = Ymin;
  max.x = Xmax;
  max.y = Ymax;
}

void Xgen::XgInitPointList() {
  GoThroughPointsPtr = 0;
}

bool Xgen::XgGiveNextPoint(Point &p) {
  if(GoThroughPointsPtr < Npoin) {
    p = E[GoThroughPointsPtr++];
    return true;
  }
  else return false;
} 

void Xgen::XgInitInteriorPointList() {
  RedrawInteriorPtr = BoundaryPtsNum;
}

bool Xgen::XgGiveNextInteriorPoint(Point &p) {
  if(RedrawInteriorPtr < Npoin) {
    p = E[RedrawInteriorPtr++];
    return true;
  }
  else return false;
}

void Xgen::XgInitBoundaryPairsList() {
  BPL->Init();
}

void Xgen::XgInitBoundaryLinesList() {
  BLL->Init();
}

bool Xgen::XgGiveNextBoundaryPair(Point &p, Point &q) {
  int A, B;
  if(BPL->GetNext(A, B)) {
    p = E[A];
    q = E[B];
    return true;
  }
  else return false;
}

bool Xgen::XgGiveNextBoundaryLine(Point &p, Point &q) {
  return BLL->GetNext(p, q);
}

void Xgen::XgInitElementList() {
  EL->Init();
}

bool Xgen::XgGiveNextElement(Point &p, Point &q, Point &r) {
  int A, B, C;
  if(EL->GetNext(A, B, C)) {
    p = E[A];
    q = E[B];
    r = E[C];
    return true;
  }
  else return false;
}

bool Xgen::XgGiveNextElement(Element &element) {
  if(EL->GetNext(element.n1, element.n2, element.n3)) return true;
  else return false;
}

void Xgen::XgForgetGrid() {
  EL->Delete();
  Nelem = 0;
  BPL->Delete();
  BPL = Boundary.CreateBoundaryPairsList();
  REMOVE_BDY_PTS_ACTIVE = true;
}

//points are randomly set
void Xgen::XgSetPointsRandom() {
  Npoin = First_Npoin;
  Nstore = First_Nstore;
  if(BoundaryPtsNum > Npoin) XgError("Xgen::XgSetPointsRandom(): Internal error (1).");
  InteriorPtsNum = Npoin - BoundaryPtsNum;
#ifdef RAND_MAX
  double Rand_max = RAND_MAX;
#else
  double Rand_max = pow(2.0, 15.0);
#endif
  srand(1);
  for(int i=BoundaryPtsNum; i<Npoin; i++) {
    do {
      E[i].x =
        rand()*(Xmax - Xmin)/Rand_max + Xmin;
      E[i].y =
	rand()*(Ymax - Ymin)/Rand_max + Ymin;
    } while(!IsInside(E[i]));
  }
  printf("Domain filled with random points, %d grid points, %d interior.\n", Npoin, InteriorPtsNum);
}

bool odd(int i) {
  if(i & 1) return true;
  else return false;
}

//points are equidistantly placed
void Xgen::XgSetPointsOverlay() {
  // Place InteriorPtsNum + Nstore new points!
  int Nx = (int)((Xmax-Xmin)/H);
  int Ny = (int)((Ymax-Ymin)/(H*sqrt(3)/2.0));
  Npoin = BoundaryPtsNum;

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
      if(IsInside(P)) {
        if(Npoin >= First_Npoin + First_Nstore) 
          XgError("Xgen::XgSetPointsOverlay(): Not enough space for new points, increase the Nstore parameter."); 
        E[Npoin++] = P;
      }
    }
  }      

  if(BoundaryPtsNum > Npoin) XgError("Xgen::XgSetPointsOverlay(): Internal error (1).");
  InteriorPtsNum = Npoin - BoundaryPtsNum;
  Nstore = First_Nstore + First_Npoin - Npoin;
  IterationPtr = BoundaryPtsNum;
  RedrawInteriorPtr = BoundaryPtsNum;
  GoThroughPointsPtr = BoundaryPtsNum;

  printf("Domain filled with overlay points, %d grid points, %d interior.\n", Npoin, InteriorPtsNum);
}

void Xgen::XgOutput(FILE *f) {
  OutputBoundaryEdgeList* bel = Boundary.CreateOutputBoundaryEdgeList();
  XgUserOutput(f, bel);
  bel -> Delete();
} 

void Xgen::XgInputPoints(FILE *f, int *error, int *mem) {
  *error = 0;
  *mem = 0;
  
  Nstore = Nstore + Npoin - BoundaryPtsNum;
  Npoin = BoundaryPtsNum;
  IterationPtr = BoundaryPtsNum;

  Point P, *E0;
  E0 = E + BoundaryPtsNum;
  while(Get(f, &P.x) && Get(f, &P.y)) if(IsInside(P)) {
    Npoin++; Nstore--;
    if(Nstore == 0) break;
    *E0 = P;
    E0++;
  }

  if(BoundaryPtsNum > Npoin) XgError("Xgen::XgInputPoints: Internal error (1).");
  InteriorPtsNum = Npoin - BoundaryPtsNum;
  return;
} 

void Xgen::XgOutputPoints(FILE *f) {
  fprintf(f, "# List of interior points:\n");

  Point P;
  GoThroughPointsPtr = BoundaryPtsNum;
  while(XgGiveNextPoint(P) == true) fprintf(f, "%g %g\n", P.x, P.y);
} 

void Undraw_point(Point);

bool Xgen::XgNextTriangle(Point &p, Point &q, Point &r, bool &finished) {

  finished = false;
  int A, B, C;

  // At first call, go through all boundary edges and remove points which 
  // are too close to the boundary.
  if(REMOVE_BDY_PTS_ACTIVE == true) {
    double coeff = 3.0;

    int was_deleted = 0;
    Point a, b;
    XgInitBoundaryLinesList();                                           
    while(XgGiveNextBoundaryLine(a, b) == true) {
      //printf("Boundary removal line: (%g %g) (%g %g).\n", a.x, a.y, b.x, b.y);

      //double h = sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
      double dh = H/coeff;
      XgInitInteriorPointList();
      Point c;
      while(XgGiveNextInteriorPoint(c) == true) {
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
          printf("Point [%g, %g] closer than h/%g to the boundary -> deleting.\n", c.x, c.y, coeff);
          was_deleted++;
          bool success = XgRemoveInteriorPoint(c);
          if(!success) 
            XgError("Xgen::XgNextTriangle: Internal error while removing points that lie too close to boundary.");
          //break;
        }
      }                
    }                 
    REMOVE_BDY_PTS_ACTIVE = false;
    if(was_deleted > 0) printf("Deleted %d points, %d grid points, %d interior.\n",
      was_deleted, Npoin, InteriorPtsNum);
  }

  BPL->GetLast(A, B);

  if(!FindNearestLeft(A, B, C)) {
    printf("FindNearestLeft() failed.\n");
    return false;
  }

  if (BPL->Delete(C, A)) {
    if (BPL->Delete(B, C)) BPL->DeleteLast();
    else BPL->ChangeLast(C, B);
  }
  else {
  BPL->ChangeLast(A, C);
    if (!BPL->Delete(B, C)) BPL->Add(C, B);
  }

  if(BPL->IsEmpty()) finished = true; 

  // Elements are positively oriented.
  EL -> Add(A, B, C);

  Nelem++;

  p = E[A];
  q = E[B];
  r = E[C];
  
  //printf("New element (%d %d %d) ... (%g %g) (%g %g) (%g %g)\n", 
  //       A, B, C, E[A].x, E[A].y, E[B].x, E[B].y, E[C].x, E[C].y);

  return true;
}

void Xgen::XgNextShift(Point *p_old, Point *p_new) {
  if(Npoin > BoundaryPtsNum) {
    if(IterationPtr >= Npoin) {
      IterationPtr = BoundaryPtsNum;
    }
    Point impuls = GetImpuls(IterationPtr);
    *p_old = E[IterationPtr];
    Shift(IterationPtr, impuls);
    *p_new = E[IterationPtr++];
  }
}

bool Xgen::XgMouseAdd(Point p) {
  if(Nstore == 0) return false;
  if(!IsInside(p)) return false;
  E[Npoin] = p;
  InteriorPtsNum++;
  Npoin++;
  Nstore--;
  IterationPtr=BoundaryPtsNum;
  printf("Point at [%g, %g] added to the end of the list, %d grid points, %d interior.\n", 
    p.x, p.y, Npoin, InteriorPtsNum);
  return true;
}

bool Xgen::XgMouseRemove(Point p_from, Point *p_where) {
  if(InteriorPtsNum == 0) return false;
  int wanted = BoundaryPtsNum;  // if it is the last interior one
  double min = 1e100;
  for(int j=BoundaryPtsNum; j < Npoin; j++) {
    Point Dr = E[j] - p_from;
    if(Dr.abs() < min) {
      min = Dr.abs();
      wanted = j;
    }
  }
  *p_where = E[wanted];
  for(int j = wanted + 1; j < Npoin; j++) E[j - 1] = E[j];
  Npoin--;
  InteriorPtsNum--;
  Nstore++;
  IterationPtr = BoundaryPtsNum;
  printf("Point no. %d at [%g, %g] removed by mouse, %d grid points, %d interior.\n", 
    wanted, E[wanted].x, E[wanted].y, Npoin, InteriorPtsNum);
  return true;
}

bool Xgen::XgRemoveInteriorPoint(Point p) {
  if(InteriorPtsNum <= 0) return false;
  int wanted = -1;  // Index of the wanted point.
  double min = 1e100;
  // Looping through interior points.
  for(int j = BoundaryPtsNum; j < Npoin; j++) {
    double dist = sqrt((E[j].x - p.x)*(E[j].x - p.x) + 
                       (E[j].y - p.y)*(E[j].y - p.y));
    if(dist < min) {
      min = dist;
      wanted = j;
    }
  }
  if(wanted == -1) XgError("Xgen::XgRemoveInteriorPoint(): Internal error.");
  for(int j = wanted; j< Npoin-1; j++) E[j] = E[j+1];
  Npoin--;
  InteriorPtsNum--;
  Nstore++;
  RedrawInteriorPtr--;

  if (this->nogui == false) Undraw_point(p);
  IterationPtr = BoundaryPtsNum;

  return true;
}

bool Xgen::XgIsEmpty() {
  return (Npoin == BoundaryPtsNum);
}

void Xgen::XgTimeInc() {
    DeltaT /= TimestepConst;
}

void Xgen::XgTimeDec() {
    DeltaT *= TimestepConst;
}

void Xgen::XgSetTimestepConstant(double constant) {
  TimestepConst = constant;
  if(TimestepConst <= 0) XgError("Xgen::XgSetTimestepConstant(): Invalid TimestepConst.");
  if(TimestepConst >= 1) XgError("Xgen::XgSetTimestepConstant(): Invalid TimestepConst.");
}

Point Xgen::XgGivePoint(long pos_in_list) {
  if(pos_in_list < 0 || pos_in_list >= Npoin)
    XgError("Xgen::XgGivePoint(): Invalid point index.");
  return E[pos_in_list];
}

Element Xgen::XgGiveElement(long pos_in_list) {
  if(pos_in_list < 0 || pos_in_list >= Nelem)
    XgError("Xgen::XgGiveElement(): Invalid element index.");
  ElemBox *ptr = EL->First;
  for(long i=0; i<pos_in_list; i++) ptr = ptr->next;
  return ptr->E;
}

/**
 **   class Xgen: protected methods
 **/

double min(double a, double b) {
  if (a < b) return a;
  else return b;
}

double max(double a, double b) {
  if (a > b) return a;
  else return b;
}

void Xgen::XgInit(char *cfg_filename) {
  if(cfg_filename == NULL) XgError("Xgen::XgInit(): Empty filename.\nVerify the number of command line parameters.");

  // Initializing variables.
  InteriorPtsNum = Npoin = First_Npoin = Nelem = BoundaryPtsNum = 0;
  Xmax = Xmin = Ymax = Ymin = 0;
  H = Area = 0;
  E = NULL;
  BPL = new BoundaryPairsList;
  BLL = new BoundaryLinesList;
  PL = new PointList;
  EL = new ElemList;
  IterationPtr = 0;
  RedrawInteriorPtr = 0;
  GoThroughPointsPtr = 0;
  TimestepConst = 0.8;

  // Spatial dimension.
  Dimension = 2;

  // Extracting the application name from the cfg filename.
  char Help_str[50];
  strcpy(Help_str, cfg_filename);
  for(int i = 0; i<50; i++) {
    if((Help_str[i] == '\0') || (Help_str[i] == '.')) {
      if((Name = (char*)malloc((i+1)*sizeof(char))) == NULL) 
        XgError("Xgen::XgInit(): Not enough memory."); 
      Help_str[i] = '\0';
      strcpy(Name, Help_str);
      break;
    }
  }

  // Setting defaults.
  First_Nstore = Nstore = 10000;
  First_DeltaT = DeltaT = 1;
  H = 1;

  // Reading user data.
  FILE *f; 
  f = fopen(cfg_filename, "rb");
  if(f == NULL) XgError("Couldn't open configuration file.");        
  XgReadData(f);
  fclose(f); 

  // Sanity checks.
  if(BoundaryPtsNum == 0) XgError("Bad configuration file: BoundaryPtsNum = 0.");
  if(DeltaT <= 0) XgError("DeltaT not initialized in XgReadData().");
 
  // Storing DeltaT in case user would want to refresh it.
  First_DeltaT = DeltaT;

  // Storing Nstore in case user would want to refresh it.
  First_Nstore = Nstore;

  // Creating initial point list.
  PL = Boundary.CreateBoundaryPointList();

  // Boundary calculates its length (taking into account that some segments are 
  // curvilinear). 
  Boundary.CalculateLength();

  // Average boundary edge length.
  double length = Boundary.GiveLength();
  this->H = length / BoundaryPtsNum;
  printf("boundary length = %g\n", length);

  // Constructing boundary -- list of pairs of point indices (for mesh geenration).
  BPL = Boundary.CreateBoundaryPairsList();
  BPL->Init();                                           

  // Constructing boundary -- the geometrical curve.
  BLL = Boundary.CreateBoundaryLinesList(this->H);
  BLL->Init();                                           

  // Calculating extrems of boundary coordinates.
  Xmin = min(BLL->First->A.x, BLL->First->B.x);            
  Xmax = max(BLL->First->A.x, BLL->First->B.x);            
  Ymin = min(BLL->First->A.y, BLL->First->B.y);            
  Ymax = max(BLL->First->A.y, BLL->First->B.y);            
  Point T1, T2;
  while(BLL->GetNext(T1, T2)) {       
    double xmin = min(T1.x, T2.x);            
    double xmax = max(T1.x, T2.x);            
    double ymin = min(T1.y, T2.y);            
    double ymax = max(T1.y, T2.y);            
    if (xmin < Xmin) Xmin = xmin;
    if (xmax > Xmax) Xmax = xmax;
    if (ymin < Ymin) Ymin = ymin;
    if (ymax > Ymax) Ymax = ymax;
  }
  BLL->Init();
  printf("Xmin = %g, Xmax = %g, Ymin = %g, Ymax = %g\n", Xmin, Xmax, Ymin, Ymax);

  // Calculating the domain's area.
  Point A, B;
  double xa, xb, ya, yb;
  Area = 0;                                                
  while(BLL->GetNext(A, B)) {                                    
    xa = A.x; ya = A.y;                                 
    xb = B.x; yb = B.y;                                 
    Area += (ya + yb - 2*Ymin)*(xa - xb)/2;                
  }                          
  BLL->Init();                                           

  // Sanity checks.
  if(Area <= 0)
    XgError("Bad boundary (or its orientation) in XgReadData().");

  // Calculating optimal number of interior points.
  int help;
  help = (int)
  ((Area/(H*H*sqrt(3)/4) - 0.9*BoundaryPtsNum)/2 + 0.5);
  if(help <= 0) InteriorPtsNum = 0;
  else InteriorPtsNum = help;
  Npoin = InteriorPtsNum + BoundaryPtsNum;

  // Initializing iteration loop.
  IterationPtr = BoundaryPtsNum;

  // Storing Npoin incase user would want to refresh it.
  First_Npoin = Npoin;

  // Copying boundary points to the beginning of the point list.
  E = (Point*)malloc((Npoin + Nstore)*sizeof(Point)); 
  if(E == NULL) XgError("Not enough memory for points.");                                         
  else for(int i = 0; i<BoundaryPtsNum; i++) E[i] = PL->Delete();

  // Setting initial points either randomly or 
  // using an equidistributed regular pattern.
  if(this->overlay) XgSetPointsOverlay();
  else XgSetPointsRandom();

  // In interactive mode we just return control, in batch mode
  // we loop over the points "steps_to_make" times,
  // generate mesh, save it to a file, and quit. 
  if (this->nogui == true) {
    // Run smoothing iterations.
    for (int s = 0; s < this->steps_to_take; s++) {
      this->XgInitInteriorPointList();
      for(int i=0; i<this->InteriorPtsNum; i++) {
        Point P1, P2;
        this->XgNextShift(&P1, &P2);
      }
      printf("Smoothing iteration %d finished.\n", s+1);
    }

    // Create mesh.
    printf("Starting the meshing algorithm.\n");
    REMOVE_BDY_PTS_ACTIVE = true;
    bool success = this->CreateMeshBatchMode();
    if (!success) XgError("Mesh algorithm failed.");

    // Save the mesh to a file.
    FILE *f = fopen("out.mesh", "w");
    if (f == NULL) XgError("Failed to open mesh file for writing.");
    OutputBoundaryEdgeList* bel = Boundary.CreateOutputBoundaryEdgeList();
    XgUserOutput(f, bel);
    bel -> Delete();
    fclose(f);

    // Exit the program.
    printf("Mesh saved to \"out.mesh\".\n");
    exit(0);
  }

  return;
}

void Xgen::XgSetTimestep(double init_timestep) {
  DeltaT = init_timestep;
}

void Xgen::XgCreateNewBoundaryComponent() {
  Boundary.CreateNewBoundaryComponent();
}

void Xgen::XgAddBoundarySegment(int marker, double x, double y, int subdiv, double alpha) {
  Boundary.AddBoundarySegment(marker, x, y, subdiv, alpha);
  BoundaryPtsNum += subdiv;
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

bool XgVertexInElem(int v, Element *e) {
  if(v == e->n1 || v == e->n2 || v == e->n3) return true;
  else return false;
}

void Xgen::XgUserOutput(FILE *f, OutputBoundaryEdgeList* bel) {
  fprintf(f, "# XGEN mesh in Hermes2D format\n");
  fprintf(f, "# Project: "); fprintf(f, "%s\n", XgGiveName());
  fprintf(f, "# Edges are positively oriented\n");

  //writing list of mesh vertices
  //(two coordinates per line)
  fprintf(f, "\n# Vertices:\nvertices =\n{\n");
  Point P;
  XgInitPointList();
  int counter = 0;
  while(XgGiveNextPoint(P) == true) {
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
  while(XgGiveNextElement(E) == true) {
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
  BoundaryEdge be;
  bel->Init();
  counter = 0;
  while(bel->GetNext(be) == true) {
    counter++;
    if (counter < XgGiveBoundaryPtsNum()) fprintf(f, "  { %d, %d , %d},\n", be.A, be.B, be.marker);
    else fprintf(f, "  { %d, %d , %d}\n", be.A, be.B, be.marker);
  }
  fprintf(f, "}\n");
}

/**
 **   class Xgen: private methods
 **/

Point Xgen::GiveImpuls(int j, int i) {
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

Point Xgen::GetImpuls(int i) {
  double I, I0;
  Point Impuls(0, 0);
  for(int j=0; j<Npoin; j++) {
    if(j != i) Impuls += GiveImpuls(j, i);
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

  if(!IsInside(E[i])) E[i] = Old_pos;
}

bool Xgen::IsInside(Point &P) {
  int Flag = 0;
  Point A, B;

  BLL->Init();
  while((BLL->GetNext(A, B)) == true) {
    if((A.x <= P.x) && (P.x < B.x)) {
      if(IsRight(A, B, P)) Flag +=1;
    }
    if((A.x >= P.x) && (P.x > B.x)) {
      if(IsLeft(A, B, P)) Flag -=1;
    }
  }
  if(Flag == 0) return false;
  else return true;
}

bool Xgen::EdgesIntersect(Point a, Point b, Point c, Point d) 
{
  double p, t, s;
  Point f = b - a, e = d - c;

  if(fabs(p = f.y*e.x - f.x*e.y) < XG_ZERO) return false;

  t = (e.y*a.x - e.x*a.y + c.y*e.x - c.x*e.y)/p;
  if(fabs(e.x) > XG_ZERO) s = (a.x + t*f.x - c.x)/e.x;
  else s = (a.y + t*f.y - c.y)/e.y;
  if(t > XG_ZERO && t < 1 - XG_ZERO && s > XG_ZERO && s < 1 - XG_ZERO) return true;
  else return false;
}

// Checks whether edges (a, c) or (b, c) intersect with 
// (possibly curvilinear) boundary.
bool Xgen::BoundaryIntersectionCheck(Point a, Point b, Point c) 
{
  Point d, e;
  BLL -> Init();
  while(BLL->GetNext(d, e)) {
    if(EdgesIntersect(a, c, d, e) || EdgesIntersect(b, c, d, e)) return true;
  }
  return false;
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

bool Xgen::FindNearestLeft(int A, int B, int &wanted) 
{
  double m, Min;
  long int wanted0 = -1;
  Min = 1e50;
  for (int C = 0; C < Npoin; C++) {
    if (C != A && C != B) {
      if (IsLeft(E[A], E[B], E[C])) {
        m = XgCriterion(E[A], E[B], E[C]);
        if(m < Min) {
          if(IsInside(E[C]) || BPL->Contains(C)) {
            if(!BoundaryIntersectionCheck(E[A], E[B], E[C])) {
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
  if(wanted0 < 0) return false;
  wanted = wanted0;
  return true;
}

bool Xgen::IsLeft(Point A, Point B, Point C) {
  double x1 = B.x - A.x,
	 y1 = B.y - A.y,
	 x2 = C.x - A.x,
	 y2 = C.y - A.y;

  if((x1*y2 - y1*x2) > 0) return true;
  else return false;
}

bool Xgen::IsRight(Point A, Point B, Point C) {
  double x1 = B.x - A.x,
	 y1 = B.y - A.y,
	 x2 = C.x - A.x,
	 y2 = C.y - A.y;

  if((x1*y2 - y1*x2) < 0) return true;
  else return false;
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
   info_ninterior_str[50],
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

void Draw_point(Point P) {
  int blackboard_x = -RS.Init_pos_x + (int)(RS.Ratio.x*P.x/RS.Measure),
   blackboard_y = RS.Init_pos_y + RS.MAX_Y - (int)(RS.Ratio.y*P.y/RS.Measure);
  Draw_point(blackboard_x, blackboard_y);
}

void Undraw_point(Point P) {
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

void Draw_boundary_and_points() {
  Point a, b;             
           
  // boundary (just the line)
  RS.RPtr->XgInitBoundaryLinesList();                                           
  while(RS.RPtr->XgGiveNextBoundaryLine(a, b) == true) {
    //Draw_point2(a);
    Draw_line(a, b);                       
  }          
  // boundary points       
  RS.RPtr->XgInitPointList();
  int count = 0;
  int bdy_pts_num = RS.RPtr->XgGiveBoundaryPtsNum();
  while(RS.RPtr->XgGiveNextPoint(a) == true && count < bdy_pts_num) {
    Draw_point2(a);
    count++;
  }
  // interior points
  RS.RPtr->XgInitInteriorPointList();
  while(RS.RPtr->XgGiveNextInteriorPoint(a) == true) Draw_point(a);
}

void Redraw_boundary() {
  Point a, b;             
           
  // just the line
  RS.RPtr->XgInitBoundaryLinesList();                                           
  while(RS.RPtr->XgGiveNextBoundaryLine(a, b) == true) {
    //Draw_point2(a);
    Draw_line(a, b);                       
  }            
  // boundary points       
  RS.RPtr->XgInitPointList();
  int count = 0;
  int bdy_pts_num = RS.RPtr->XgGiveBoundaryPtsNum();
  while(RS.RPtr->XgGiveNextPoint(a) == true && count < bdy_pts_num) {
    Draw_point2(a);
    count++;
  }
}

void Redraw_grid() {                                           
  Point a, b, c;                                               
  RS.RPtr->XgInitElementList();   
  while(RS.RPtr->XgGiveNextElement(a, b, c) == true) {
    Draw_line(a, b);
    Draw_line(a, c);
    Draw_line(b, c);   
  }                 
}

void Redraw_blackboard() {
  if(RS.Grid_working || RS.Grid_ready) Redraw_grid();
  Draw_boundary_and_points();
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
   XmNsensitive, false,
   NULL
  );
  XtVaSetValues(
   RS.gridsave,
   XmNsensitive, false,
   NULL
  );     
  XtVaSetValues(
   RS.quit,
   XmNsensitive, false,
   NULL
  );     
  XtVaSetValues(
   RS.blackboard,
   XmNsensitive, false,
   NULL
  );
}

void ActivateWindow() {
  XtVaSetValues(
   RS.initforget,
   XmNsensitive, true,
   NULL
  );     
  XtVaSetValues(
   RS.gridsave,
   XmNsensitive, true,
   NULL
  );     
  XtVaSetValues(
   RS.quit,
   XmNsensitive, true,
   NULL
  );     
  XtVaSetValues(
   RS.blackboard,
   XmNsensitive, true,
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
   XmNsensitive, true,
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
   XmNsensitive, true,
   NULL
  );
  XtVaSetValues(
   RS.initd_cancel,
   XmNsensitive, true,
   NULL
  );
  XtVaSetValues(
   RS.initd_set,
   XmNsensitive, true,
   NULL
  );
  XtVaSetValues(
   RS.initd_label,
   XmNsensitive, true,
   NULL
  );
}

void SavePointsErrorDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtUnmanageChild((Widget)client_data);

/* DO NOT REMOVE THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savepointsdialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, true,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savepointsdialog,
   XmNsensitive, true,
   NULL
  );
}

void InitDiskErrorDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  //w, client_data, call_data,
  XtUnmanageChild((Widget)client_data);

/* DO NOT REMOVE THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.initd_load_d, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, true,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.initd_load_d,
   XmNsensitive, true,
   NULL
  );
}

void SaveDiskErrorDialogOK(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  //w, client_data, call_data,
  XtUnmanageChild((Widget)client_data);

/* DO NOT REMOVE THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, true,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savedialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, true,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savedialog,
   XmNsensitive, true,
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
  Draw_boundary_and_points();
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
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savepointsdialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savepointsdialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, false,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savepointsdialog,
   XmNsensitive, false,
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

/* DO NOT REMOVE THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.initd_load_d, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.initd_load_d, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, false,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.initd_load_d,
   XmNsensitive, false,
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

/* DO NOT REMOVE THIS
  Widget child;
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_APPLY_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_CANCEL_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DIR_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_DEFAULT_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_FILTER_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_LIST_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_OK_BUTTON
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SELECTION_LABEL
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_SEPARATOR
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
  child = XmFileSelectionBoxGetChild(
   RS.savedialog, 
   XmDIALOG_TEXT
  );
  XtVaSetValues(
   child,
   XmNsensitive, false,
   NULL
  );
//  child = XmFileSelectionBoxGetChild(
//   RS.savedialog, 
//   XmDIALOG_WORK_AREA
//  );
//  XtVaSetValues(
//   child,
//   XmNsensitive, false,
//   NULL
//  );
*/

  XtVaSetValues(
   RS.savedialog,
   XmNsensitive, false,
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
     XmNsensitive, true,
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
            int ninterior = RS.RPtr->XgGiveInteriorPtsNum();
            strcpy(points_info, RS.info_ninterior_str);
            sprintf(number, "%d", ninterior);
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
            Draw_boundary_and_points();
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
    Draw_boundary_and_points();
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
   XmNsensitive, false,
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
   XmNsensitive, false,
   NULL
  );
  XtVaSetValues(
   RS.initd_set,
   XmNsensitive, false,
   NULL
  );
  XtVaSetValues(
   RS.initd_cancel,
   XmNsensitive, false,
   NULL
  );
  XtVaSetValues(
   RS.initd_label,
   XmNsensitive, false,
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
  RS.RPtr->XgSetPointsOverlay();
  int npoin = RS.RPtr->XgGiveNpoin();
  char points_info[255];
  strcpy(points_info, RS.info_npoin_str);
  char number[20];
  sprintf(number, "%d", npoin);
  strcat(points_info, number);
  XmString info_l2_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
  int ninterior = RS.RPtr->XgGiveInteriorPtsNum();
  strcpy(points_info, RS.info_ninterior_str);
  sprintf(number, "%d", ninterior);
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
  Draw_boundary_and_points();
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
   XmNsensitive, true,
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
   XmNsensitive, false,
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
     XmNsensitive, false,
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
     XmNsensitive, true,
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
     XmNsensitive, true,
     XmNlabelString, RS.grid_name,
     NULL
    );     
    XtVaSetValues(
     RS.blackboard,
     XmNsensitive, true,
     NULL
    );
    XClearWindow(
      XtDisplay(RS.blackboard),  
      XtWindow(RS.blackboard)    
    );    
    Draw_boundary_and_points();
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
     XmNsensitive, false,
     NULL
    );          
    XtVaSetValues(
     RS.blackboard,
     XmNsensitive, false,
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
      Draw_point(P);
      int npoin = RS.RPtr->XgGiveNpoin();
      char points_info[255];
      strcpy(points_info, RS.info_npoin_str);
      char number[20];
      sprintf(number, "%d", npoin);
      strcat(points_info, number);
      XmString info_l2_label = 
       XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
      int ninterior = RS.RPtr->XgGiveInteriorPtsNum();
      strcpy(points_info, RS.info_ninterior_str);
      sprintf(number, "%d", ninterior);
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
      Undraw_point(Ret);
      int npoin = RS.RPtr->XgGiveNpoin();
      char points_info[255];
      strcpy(points_info, RS.info_npoin_str);
      char number[20];
      sprintf(number, "%d", npoin);
      strcat(points_info, number);
      XmString info_l2_label = 
       XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);
      int ninterior = RS.RPtr->XgGiveInteriorPtsNum();
      strcpy(points_info, RS.info_ninterior_str);
      sprintf(number, "%d", ninterior);
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
  if(RS.Grid_ready || RS.Grid_failed) return false;
  if(RS.Grid_working) {
    Point a, b, c;
    bool finished;
    if(RS.RPtr->XgNextTriangle(a, b, c, finished)) {
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
      if(finished == true) {
        RS.Grid_working = 0;
        RS.Grid_ready = 1;
        REMOVE_BDY_PTS_ACTIVE = true;
        if(!RS.IsInfo) {
          XtVaSetValues(
           RS.info_l1,
           XmNlabelString, RS.info_l1_label3,
           NULL
          );
        }
        XtVaSetValues(
         RS.gridsave,
         XmNsensitive, true,
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
        Undraw_point(P1);
        Draw_point(P2);
      }
    }
    if(++Bound_count > RS.BoundaryRedrawInterval) {
      Redraw_boundary(); 
      Bound_count = 0;
    }      
  }
  return(false); 
}

bool Xgen::CreateMeshBatchMode() {
 
  Point a, b, c;
  bool finished, success;
  do {
    success = this->XgNextTriangle(a, b, c, finished);
    if (!success) return false;
  }
  while (finished == false);

  return true;
}

void MenuCancel(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtUnmanageChild(RS.menu);
  XtVaSetValues(
   RS.menu_button,
   XmNsensitive, true,
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
   XmNsensitive, true,
   NULL
  );
  RS.IsInfo = 1;
}

void Info(
 Widget w, XtPointer client_data, XtPointer call_data
) {
  XtVaSetValues(
   RS.info,
   XmNsensitive, false,
   NULL
  );
  RS.IsInfo = 0;
  Widget ok;  

  int i = 0;
  Arg args[20];
  
  XtSetArg(args[i], XtNtitle, "Info"); i++;
  XtSetArg(args[i], XmNminWidth, 300); i++;
  XtSetArg(args[i], XmNminHeight, 230); i++;
  XtSetArg(args[i], XmNwidth, 300); i++;
  XtSetArg(args[i], XmNheight, 230); i++;
  XtSetArg(args[i], XmNresizable, false); i++;

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

  int ninterior = RS.RPtr->XgGiveInteriorPtsNum();
  strcpy(points_info, RS.info_ninterior_str);
  sprintf(number, "%d", ninterior);
  strcat(points_info, number);
  XmString info_l3_label = 
   XmStringCreate(points_info, XmFONTLIST_DEFAULT_TAG);

  strcpy(points_info, "Boundary points: ");
  sprintf(number, "%d", npoin - ninterior);
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
  XtSetArg(args[i], XmNautoUnmanage, false); i++;
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
   XmNsensitive, RS.IsAbout ? true : false,
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
   XmNsensitive, RS.IsInfo ? true : false,
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
   XmNsensitive, false,
   NULL
  );
}

/*
 *    Graphic setup:     
 */

static void SetUpGraphics() {
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
  strcpy((char*)msg, "The meshing algorithm failed.");
  RS.info_l1_label4 = XmStringCreate(
    msg, XmFONTLIST_DEFAULT_TAG);
  strcpy(RS.info_npoin_str, "Number of points: ");
  strcpy(RS.info_ninterior_str, "Interior points: ");
  strcpy(RS.info_nelem_str, "Number of elements: ");
  strcpy(RS.info_timestep_str, "Timestep: ");

  //initialization of area extrems:
  Point min, max;
  User -> XgGiveLimits(min, max);
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
  SetUpGraphics();

  //running the XGen-Window:
  XtAppMainLoop(app_context);
}

















