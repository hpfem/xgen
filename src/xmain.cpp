/*
 *   xmain.cpp
 */

# include <stdio.h>
# include <string.h>
# include <strings.h>

/** 
 **   simple example applications
 **   29.6.1995, Prague
 **/

# include "xgen.h"
# include "disc.h"

/*
 *   class xsquare:
 */

class xsquare: public Xgen {
  public:
  xsquare(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xsquare::XgReadData(FILE *f) {
  int N; 
  if(!Get(f, &N)) XgError("Couldn't read N.");
  XgAddEdge(3, 0, 0, N);
  XgAddEdge(2, 1, 0, N);
  XgAddEdge(4, 1, 1, N);
  XgAddEdge(1, 0, 1, N);
}

/*
 *   class xhole:
 */

class xhole: public Xgen {
  public:
  xhole(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xhole::XgReadData(FILE *f) {
  double H;
  if(!Get(f, &H)) XgError("Couldn't read H.");

  double R;
  if(!Get(f, &R)) XgError("Couldn't read R.");
  int N2 = (int)(2*M_PI*R + 0.5); 

  int N;
  if(!Get(f, &N)) XgError("Couldn't read N.");
  XgAddEdge(3, 0, 0, 2*N);
  XgAddEdge(2, 1.5*N*H, 0, N);
  XgAddEdge(4, 1.5*N*H, N*H, 2*N);
  XgAddEdge(1, 0, N*H, N);
  XgNewComponent();
  for(int i=0; i<N2; i++) {
    XgAddEdge(5, 0.5*H*N + R*H*cos(i*2*M_PI/N2), 0.5*H*N - R*H*sin(i*2*M_PI/N2));
  }
}

/*
 *   class xcirc:
 */

class xcirc: public Xgen {
  public:
  xcirc(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xcirc::XgReadData(FILE *f) {
  double Rout; int Nout; 
  if(!Get(f, &Rout)) XgError("Couldn't read Rout.");
  if(!Get(f, &Nout)) XgError("Couldn't read Nout.");

  for(int i=0; i<Nout; i++) {
    XgAddEdge(1, Rout*cos(i*2*M_PI/Nout), Rout*sin(i*2*M_PI/Nout));
  }

  double H = 2*M_PI*Rout/Nout;

  double Rins;
  if(!Get(f, &Rins)) XgError("Couldn't read Rins.");
  int Nins = (int)(2*M_PI*Rins/H + 0.5);

  XgNewComponent();
  for(int i=0; i<Nins; i++) {
    XgAddEdge(2, Rins*cos(i*2*M_PI/Nins), -Rins*sin(i*2*M_PI/Nins));
  }
}

/*
 *   class xgamm:
 */

class xgamm: public Xgen {
  public:
  xgamm(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xgamm::XgReadData(FILE *f) {
  double H;
  if(!Get(f, &H)) XgError("Couldn't read H.");

  int N1, N2, N3, N4; double X5;
  if(!Get(f, &N1)) XgError("Couldn't read N1.");
  if(!Get(f, &N2)) XgError("Couldn't read N2.");
  if(!Get(f, &N3)) XgError("Couldn't read N3.");
  if(!Get(f, &N4)) XgError("Couldn't read N4.");
  if(!Get(f, &X5)) XgError("Couldn't read X5.");

  double k = (0.125*N3*N3/X5 - 0.5*X5)*H;
  double r = k + X5*H;

  XgAddEdge(3, 0, 0, N2);
  for(int i=0; i<N3; i++) {
    double x = (i - 0.5*N3)*H;
    XgAddEdge(3, (N2+i)*H, sqrt(r*r - x*x) - k);
  }
  XgAddEdge(3, (N2+N3)*H, 0, N4);
  XgAddEdge(2, (N2+N3+N4)*H, 0, N1);
  XgAddEdge(4, (N2+N3+N4)*H, N1*H, N2+N3+N4);
  XgAddEdge(1, 0, N1*H, N1);
}

/*
 *   class xnozzle0:
 */

class xnozzle0: public Xgen {
  public:
  xnozzle0(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xnozzle0::XgReadData(FILE *f) {
  double R1, R2, H1, H2, L, Dh;

  if(!Get(f, &R1)) XgError("Couldn't read R1.");
  if(!Get(f, &R2)) XgError("Couldn't read R2.");
  if(!Get(f, &H1)) XgError("Couldn't read H1.");
  if(!Get(f, &H2)) XgError("Couldn't read H2.");
  if(!Get(f, &L)) XgError("Couldn't read L.");
  if(!Get(f, &Dh)) XgError("Couldn't read Delta_h.");

  int N1, N2, N3, N4, N5, N6;
  N1 = (int)((R1 + R2 + L)/Dh + 0.5);
  N2 = (int)((H2)/Dh + 0.5);
  N3 = (int)(sqrt(L*L + (H1 - H2)*(H1 - H2))/Dh + 0.5);
  N4 = (int)((M_PI*R2/2.0)/Dh + 0.5);
  N5 = (int)((M_PI*R1/2.0)/Dh + 0.5);
  N6 = (int)((H1 + R1 + R2)/Dh + 0.5);
  
  double S1_x = 0;
  double S1_y = H1 + R2;
  double S2_x = R1 + R2;
  double S2_y = H1 + R2;

  XgAddEdge(3, 0, 0, N1);
  XgAddEdge(2, R1 + R2 + L, 0, N2);
  XgAddEdge(4, R1 + R2 + L, H2, N3);
  
  double angle2 = M_PI/(2.0*(double)N4);
  for(int i=0; i<N4; i++) {
    double px, py;
    px = S2_x - R2*sin(i*angle2);
    py = S2_y - R2*cos(i*angle2);
    XgAddEdge(4, px, py, 1);
  }
 
  double angle1 = M_PI/(2.0*(double)N5);
  for(int i=0; i<N5; i++) {
    double px, py;
    px = S1_x + R1*cos(i*angle1);
    py = S1_y + R1*sin(i*angle1);
    XgAddEdge(4, px, py, 1);
  }
 
  XgAddEdge(1, 0, R1 + R2 + H1, N6);
}

/*
 *   class xduese:
 */

class xduese: public Xgen {
  public:
  xduese(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xduese::XgReadData(FILE *f) {
  double R1, R2, H1, H2, L, Dh;
  double D1, D2, D3, D4, D5;

  if(!Get(f, &R1)) XgError("Couldn't read R1.");
  if(!Get(f, &R2)) XgError("Couldn't read R2.");
  if(!Get(f, &H1)) XgError("Couldn't read H1.");
  if(!Get(f, &H2)) XgError("Couldn't read H2.");
  if(!Get(f, &L)) XgError("Couldn't read L.");
  if(!Get(f, &D1)) XgError("Couldn't read D1.");
  if(!Get(f, &D2)) XgError("Couldn't read D2.");
  if(!Get(f, &D3)) XgError("Couldn't read D3.");
  if(!Get(f, &D4)) XgError("Couldn't read D4.");
  if(!Get(f, &D5)) XgError("Couldn't read D5.");
  if(!Get(f, &Dh)) XgError("Couldn't read Delta_h.");

  int N2, N3, N4, N5, N6;
  //N1 = (int)((R1 + R2 + L)/Dh + 0.5);
  N2 = (int)((H2)/Dh + 0.5);
  N3 = (int)(sqrt(L*L + (H1 - H2)*(H1 - H2))/Dh + 0.5);
  N4 = (int)((M_PI*R2/2.0)/Dh + 0.5);
  N5 = (int)((M_PI*R1/2.0)/Dh + 0.5);
  N6 = (int)((H1 + R1 + R2)/Dh + 0.5);
  
  int ND1, ND2, ND3, ND4, ND5;
  ND1 = (int)(D1/Dh + 0.5);
  ND2 = (int)(D2/Dh + 0.5);
  ND3 = (int)(D3/Dh + 0.5);
  ND4 = (int)(D4/Dh + 0.5);
  ND5 = (int)(D5/Dh + 0.5);

  double S1_x = 0;
  double S1_y = H1 + R2;
  double S2_x = R1 + R2;
  double S2_y = H1 + R2;

  //spodek rezervoaru s tryskou
  double angle1 = M_PI/(2.0*(double)N5);
  for(int i=0; i<N5; i++) {
    double px, py;
    px = S1_x + R1*sin(i*angle1);
    py = S1_y - 2*(H1 + R2) - R1*cos(i*angle1);
    XgAddEdge(3, px, py, 1);
  }

  double angle2 = M_PI/(2.0*(double)N4);
  for(int i=0; i<N4; i++) {
    double px, py;
    px = S2_x - R2*cos(i*angle2);
    py = S2_y - 2*(H1 + R2) + R2*sin(i*angle2);
    XgAddEdge(3, px, py, 1);
  }

  XgAddEdge(3, R1 + R2, -H1, N3);

  //prostor za tryskou
  XgAddEdge(3, R1 + R2 + L, -H2, ND1);
  XgAddEdge(3, R1 + R2 + L, -H2 - D1, ND2);
  XgAddEdge(3, R1 + R2 + L + D2, -H2 - D1, ND3);
  //toto je hranice dolniho rezervoaru:
  XgAddEdge(5, R1 + R2 + L + D2, -H2 - D1 - D3, ND4);
  XgAddEdge(3, R1 + R2 + L + D2 + D4, -H2 - D1 - D3, ND3);
  XgAddEdge(3, R1 + R2 + L + D2 + D4, -H2 - D1, ND5);
  //toto je vystup na pravem konci oblasti
  XgAddEdge(2, R1 + R2 + L + D2 + D4 + D5, -H2 - D1, 2*N2 + 2*ND1);
  XgAddEdge(4, R1 + R2 + L + D2 + D4 + D5, H2 + D1, ND2 + ND4 + ND5);
  XgAddEdge(4, R1 + R2 + L, +H2 + D1, ND1);

  //vrsek rezervoaru s tryskou
  XgAddEdge(4, R1 + R2 + L, H2, N3);
  
  for(int i=0; i<N4; i++) {
    double px, py;
    px = S2_x - R2*sin(i*angle2);
    py = S2_y - R2*cos(i*angle2);
    XgAddEdge(4, px, py, 1);
  }
 
  for(int i=0; i<N5; i++) {
    double px, py;
    px = S1_x + R1*cos(i*angle1);
    py = S1_y + R1*sin(i*angle1);
    XgAddEdge(4, px, py, 1);
  }
 
  XgAddEdge(1, 0, R1 + R2 + H1, 2*N6);
}

/*
 *   class xduese2:
 */

class xduese2: public Xgen {
  public:
  xduese2(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xduese2::XgReadData(FILE *f) {
  double R1, R2, H1, H2, L, Dh;
  double D1, D2, D3, D4, D5;

  if(!Get(f, &R1)) XgError("Couldn't read R1.");
  if(!Get(f, &R2)) XgError("Couldn't read R2.");
  if(!Get(f, &H1)) XgError("Couldn't read H1.");
  if(!Get(f, &H2)) XgError("Couldn't read H2.");
  if(!Get(f, &L)) XgError("Couldn't read L.");
  if(!Get(f, &D1)) XgError("Couldn't read D1.");
  if(!Get(f, &D2)) XgError("Couldn't read D2.");
  if(!Get(f, &D3)) XgError("Couldn't read D3.");
  if(!Get(f, &D4)) XgError("Couldn't read D4.");
  if(!Get(f, &D5)) XgError("Couldn't read D5.");
  if(!Get(f, &Dh)) XgError("Couldn't read Delta_h.");

  int N2, N3, N4, N5, N6;
  //N1 = (int)((R1 + R2 + L)/Dh + 0.5);
  N2 = (int)((H2)/Dh + 0.5);
  N3 = (int)(sqrt(L*L + (H1 - H2)*(H1 - H2))/Dh + 0.5);
  N4 = (int)((M_PI*R2/2.0)/Dh + 0.5);
  N5 = (int)((M_PI*R1/2.0)/Dh + 0.5);
  N6 = (int)((H1 + R1 + R2)/Dh + 0.5);
  
  int ND1, ND2, ND3, ND4, ND5;
  ND1 = (int)(D1/Dh + 0.5);
  ND2 = (int)(D2/Dh + 0.5);
  ND3 = (int)(D3/Dh + 0.5);
  ND4 = (int)((2*R1 + 2*R2 + D4)/Dh + 0.5);
  ND5 = (int)(D5/Dh + 0.5);

  double S1_x = 0;
  double S1_y = H1 + R2;
  double S2_x = R1 + R2;
  double S2_y = H1 + R2;

  //spodek rezervoaru s tryskou
  double angle1 = M_PI/(2.0*(double)N5);
  for(int i=0; i<N5; i++) {
    double px, py;
    px = S1_x + R1*sin(i*angle1);
    py = S1_y - 2*(H1 + R2) - R1*cos(i*angle1);
    XgAddEdge(3, px, py, 1);
  }

  double angle2 = M_PI/(2.0*(double)N4);
  for(int i=0; i<N4; i++) {
    double px, py;
    px = S2_x - R2*cos(i*angle2);
    py = S2_y - 2*(H1 + R2) + R2*sin(i*angle2);
    XgAddEdge(3, px, py, 1);
  }

  XgAddEdge(3, R1 + R2, -H1, N3);

  //prostor za tryskou
  XgAddEdge(3, R1 + R2 + L, -H2, ND1);
  XgAddEdge(3, R1 + R2 + L, -H2 - D1, ND2);
  XgAddEdge(3, R1 + R2 + L + D2, -H2 - D1, ND3);

  double S1_xa, S1_ya, S2_xa, S2_ya;

  S2_xa =  R1 + L + D2;
  S2_ya = -H2 - D1 - D3;

  for(int i=0; i<N4; i++) {
    double px, py;
    px = S2_xa + R2*cos(i*angle2);
    py = S2_ya - R2*sin(i*angle2);
    XgAddEdge(3, px, py, 1);
  }

  S1_xa = S2_xa;
  S1_ya = S2_ya - R1 - R2;

  for(int i=0; i<N5; i++) {
    double px, py;
    px = S1_xa - R1*sin(i*angle1);
    py = S1_ya + R1*cos(i*angle1);
    XgAddEdge(3, px, py, 1);
  }

  //toto je hranice dolniho rezervoaru:
  XgAddEdge(5, S1_xa - R1, S1_ya, ND4);

  for(int i=0; i<N5; i++) {
    double px, py;
    px = S1_xa + 2*R2 + D4 + R1*cos(i*angle1);
    py = S1_ya + R1*sin(i*angle1);
    XgAddEdge(3, px, py, 1);
  }

  for(int i=0; i<N4; i++) {
    double px, py;
    px = S2_xa + 2*R2 + D4 - R2*sin(i*angle2);
    py = S2_ya - R2*cos(i*angle2);
    XgAddEdge(3, px, py, 1);
  }

  XgAddEdge(3, R1 + R2 + L + D2 + D4, -H2 - D1 - D3, ND3);
  XgAddEdge(3, R1 + R2 + L + D2 + D4, -H2 - D1, ND5);
  //toto je vystup na pravem konci oblasti
  XgAddEdge(2, R1 + R2 + L + D2 + D4 + D5, -H2 - D1, 2*N2 + 2*ND1);
  XgAddEdge(4, R1 + R2 + L + D2 + D4 + D5, H2 + D1, ND2 + ND4 + ND5);
  XgAddEdge(4, R1 + R2 + L, +H2 + D1, ND1);

  //vrsek rezervoaru s tryskou
  XgAddEdge(4, R1 + R2 + L, H2, N3);
  
  for(int i=0; i<N4; i++) {
    double px, py;
    px = S2_x - R2*sin(i*angle2);
    py = S2_y - R2*cos(i*angle2);
    XgAddEdge(4, px, py, 1);
  }
 
  for(int i=0; i<N5; i++) {
    double px, py;
    px = S1_x + R1*cos(i*angle1);
    py = S1_y + R1*sin(i*angle1);
    XgAddEdge(4, px, py, 1);
  }
 
  XgAddEdge(1, 0, R1 + R2 + H1, 2*N6);
}

/*
 *   class xstep:
 */

class xstep: public Xgen {
  public:
  xstep(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xstep::XgReadData(FILE *f) {
  double H;
  if(!Get(f, &H)) XgError("Couldn't read H.");

  int N1, N2, N3, N4;
  if(!Get(f, &N1)) XgError("Couldn't read N1.");
  if(!Get(f, &N2)) XgError("Couldn't read N2.");
  if(!Get(f, &N3)) XgError("Couldn't read N3.");
  if(!Get(f, &N4)) XgError("Couldn't read N4.");

  XgAddEdge(3, 0, 0, N2);
  XgAddEdge(3, N2*H, 0, N1-N4);
  XgAddEdge(3, N2*H, (N1-N4)*H, N3);
  XgAddEdge(2, (N2+N3)*H, (N1-N4)*H, N4);
  XgAddEdge(4, (N2+N3)*H, N1*H, N2+N3);
  XgAddEdge(1, 0, N1*H, N1);
}

/*
 *   class xlist:
 */

class xlist: public Xgen {
  public:
  xlist(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xlist::XgReadData(FILE *f) {
  int index, lines, end = 0;
  Point A, B, C;
  char test[20];

  while(Get(f, &index)) {
    if(!Get(f, &A)) XgError("Couldn't read a point.");
    if(!Get(f, &lines)) XgError("Couldn't read a lines number.");
    XgNewComponent();
    XgAddEdge(index, A, lines);
    C = A;
    while(end = !Get(f, test), (end || test[0] == '#') ? 0:1) {
      index = atoi(test);
      if(!Get(f, &B)) XgError("Couldn't read a point.");
      if(!Get(f, &lines)) XgError("Couldn't read a lines number.");
      XgAddEdge(index, B, lines); 
      A=B;
    } 
  }
}

/*
 *   class xnozzle:
 */

double Radius(double x) {
  //NOTE: this function defines the radius r(x) of the nozzle

  //YOUR DEFINITION OF THE RADIUS
  double radius;

  radius = -sin(M_PI*10.*x)/250 + x/100. + 0.01;
  if(x >= -0.05 && x < 0.05) 
     radius = -cos(M_PI*10.*(x-0.05))/50 + 0.0265; //+ 0.0165;

  //END OF YOUR RADIUS DEFINITION

  return radius;
}

class xnozzle: public Xgen {
  public:
  xnozzle(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xnozzle::XgReadData(FILE *f) {
  int N;
  double x_left, x_right;
 
  if(!Get(f, &x_left)) XgError("Couldn't read x_left.");
  if(x_left >= 0.05) XgError("x_left must be < 0.05");

  if(!Get(f, &x_right)) XgError("Couldn't read x_right.");
  if(x_right <= 0.05) XgError("x_right must be > 0.05");

  if(!Get(f, &N)) XgError("Couldn't read N.");

  //number of inlet abscissas
  int N1;
  N1 = (int)(N*Radius(x_left)/(x_right - x_left) + 0.5);

  //number of outlet abscissas
  int N2;
  N2 = (int)(N*Radius(x_right)/(x_right - x_left) + 0.5);

  XgAddEdge(4, x_left, 0, N);
  XgAddEdge(3, x_right, 0, N2);
  for(int i=N; i>0; i--) {
    double x;
    x = x_left + (double)i*(x_right - x_left)/N;
    XgAddEdge(4, x, Radius(x), 1);
  }
  XgAddEdge(1, x_left, Radius(x_left), N1);
}

class xspir2d: public Xgen {
  public:
  xspir2d(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xspir2d::XgReadData(FILE *f) {
  double x_left = 0, x_right = 2;
 
  //division in the x-direction
  int Nx;
  if(!Get(f, &Nx)) XgError("Couldn't read Nx.");

  //number of inlet abscissas
  int Ny;
  Ny = Nx/4;

  //definition of the shape
  double hx = (x_right - x_left)/Nx;
  double x;
  for(int i=0; i<Nx; i++) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x));
    XgAddEdge(3, x, 0.25*cos(M_PI*x), 1);
  }
  x = x_right;
  //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x));
  XgAddEdge(2, x_right, 0.25*cos(M_PI*x_right), Ny);
  for(int i=Nx; i>0; i--) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x));
    XgAddEdge(4, x, 0.25*cos(M_PI*x) + 0.5, 1);
  }
  x = x_left;
  //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x) + 0.5);
  XgAddEdge(1, x, 0.25*cos(M_PI*x) + 0.5, Ny);
}

double sep_low(double x) {
  return 0.2*x*cos(M_PI*x) + x/10;
}

double sep_up(double x) {
  return 0.2*cos(M_PI*x-0.2) + 0.4 + x/10;
}

class xsep2d: public Xgen {
  public:
  xsep2d(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xsep2d::XgReadData(FILE *f) {
  double x_left = 0.2, x_right = 2.1;
 
  //division in the x-direction
  int Nx;
  if(!Get(f, &Nx)) XgError("Couldn't read Nx.");

  //number of inlet abscissas
  int Ny1;
  Ny1 = (int)(Nx*(sep_up(x_left) - sep_low(x_left))/(x_right - x_left) + 0.5);

  //number of outlet abscissas
  int Ny2;
  Ny2 = (int)(Nx*(sep_up(x_right) - sep_low(x_right))/(x_right - x_left) + 0.5);

  //definition of the shape
  double hx = (x_right - x_left)/Nx;
  double x;
  for(int i=0; i<Nx; i++) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, sep_low(x));
    XgAddEdge(3, x, sep_low(x), 1);
  }
  x = x_right;
  //printf("corner x_left low %g, %g\n", x_left, sep_low(x_left));
  //printf("corner x_left up%g, %g\n", x_left, sep_up(x_left));
  //printf("corner x_right low %g, %g\n", x_right, sep_low(x_right));
  //printf("corner x_right up  %g, %g\n", x_right, sep_up(x_right));
  XgAddEdge(2, x_right, sep_low(x), Ny2);
  for(int i=Nx; i>0; i--) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, sep_up(x));
    XgAddEdge(4, x, sep_up(x), 1);
  }
  x = x_left;
  //printf("starting point %g, %g\n", x, sep_up(x));
  XgAddEdge(1, x, sep_up(x), Ny1);
}

/*
 *   class xmunich:
 */

class xmunich: public Xgen {
  public:
  xmunich(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xmunich::XgReadData(FILE *f) {
  double R, L1, L2, H;
  if(!Get(f, &R)) XgError("Couldn't read R.");
  if(!Get(f, &L1)) XgError("Couldn't read L1.");
  if(!Get(f, &L2)) XgError("Couldn't read L2.");
  if(!Get(f, &H)) XgError("Couldn't read H.");

  double L3 = L1/20;
  double L4 = L2/10;

  int N1 = (int)(L1/H);
  int N2 = (int)(L2/H);
  int NR = (int)(2*M_PI*R/H);
  int N3 = (int)(L3/H);
  int N4 = (int)(L4/H);
  double L6 = 4*L4;
  int N6 = (int)(L6/H);
  int N5 = (int)(sqrt(L3*L3 + 0.25*(L2-L6)*(L2-L6))/H);

  XgAddEdge(3, -L3, 0, N1 + N3);
  XgAddEdge(2, L1, 0, N5);
  XgAddEdge(2, L1+2*L3, 0.5*(L2-L6), N3);
  XgAddEdge(2, L1+L3+2*L3, 0.5*(L2-L6), N6);
  XgAddEdge(2, L1+L3+2*L3, 0.5*(L2-L6)+L6, N3);
  XgAddEdge(2, L1+2*L3, 0.5*(L2-L6)+L6, N5);
  XgAddEdge(4, L1, L2, N1 + N3);
  XgAddEdge(4, -L3, L2, N4);
  XgAddEdge(4, -L3, L2 - L4, N3);

  XgAddEdge(1, 0, L2 - L4, N2 - 2* N4);
  XgAddEdge(1, 0, L4, N3);
  XgAddEdge(1, -L3, L4, N4);

  XgNewComponent();
  for(int i=0; i<NR; i++) {
    XgAddEdge(5, 0.3*L1 + R*cos(i*2*M_PI/NR), 0.5*L2 - R*sin(i*2*M_PI/NR));
  }
  XgNewComponent();
  for(int i=0; i<NR; i++) {
    XgAddEdge(6, 0.6*L1 + R*cos(i*2*M_PI/NR), 0.5*L2 - R*sin(i*2*M_PI/NR));
  }
  XgNewComponent();
  for(int i=0; i<NR; i++) {
    XgAddEdge(6, 0.9*L1 + R*cos(i*2*M_PI/NR), 0.5*L2 - R*sin(i*2*M_PI/NR));
  }

  double D1 = L2/3;
  double D2 = H;
  int ND1 = (int)(D1/H);
  int ND2 = (int)(D2/H);
  
  double alpha1 = M_PI/5.0;

  double S1 = 0.15*L1;
  double S2 = 0.25*L2;
  XgNewComponent();
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha1) + 0.5*D2*sin(alpha1), 
            S2 - 0.5*D1*sin(alpha1) - 0.5*D2*cos(alpha1), 
            ND2);
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha1) - 0.5*D2*sin(alpha1), 
            S2 - 0.5*D1*sin(alpha1) + 0.5*D2*cos(alpha1), 
            ND1);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha1) - 0.5*D2*sin(alpha1), 
            S2 + 0.5*D1*sin(alpha1) + 0.5*D2*cos(alpha1), 
            ND2);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha1) + 0.5*D2*sin(alpha1), 
            S2 + 0.5*D1*sin(alpha1) - 0.5*D2*cos(alpha1), 
            ND1);
  
  S1 = 0.45*L1;
  S2 = 0.25*L2;
  XgNewComponent();
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha1) + 0.5*D2*sin(alpha1), 
            S2 - 0.5*D1*sin(alpha1) - 0.5*D2*cos(alpha1), 
            ND2);
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha1) - 0.5*D2*sin(alpha1), 
            S2 - 0.5*D1*sin(alpha1) + 0.5*D2*cos(alpha1), 
            ND1);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha1) - 0.5*D2*sin(alpha1), 
            S2 + 0.5*D1*sin(alpha1) + 0.5*D2*cos(alpha1), 
            ND2);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha1) + 0.5*D2*sin(alpha1), 
            S2 + 0.5*D1*sin(alpha1) - 0.5*D2*cos(alpha1), 
            ND1);
  
  S1 = 0.75*L1;
  S2 = 0.25*L2;
  XgNewComponent();
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha1) + 0.5*D2*sin(alpha1), 
            S2 - 0.5*D1*sin(alpha1) - 0.5*D2*cos(alpha1), 
            ND2);
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha1) - 0.5*D2*sin(alpha1), 
            S2 - 0.5*D1*sin(alpha1) + 0.5*D2*cos(alpha1), 
            ND1);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha1) - 0.5*D2*sin(alpha1), 
            S2 + 0.5*D1*sin(alpha1) + 0.5*D2*cos(alpha1), 
            ND2);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha1) + 0.5*D2*sin(alpha1), 
            S2 + 0.5*D1*sin(alpha1) - 0.5*D2*cos(alpha1), 
            ND1);
  

  double alpha2 = M_PI - alpha1;

  S1 = 0.15*L1;
  S2 = 0.75*L2;
  XgNewComponent();
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha2) + 0.5*D2*sin(alpha2), 
            S2 - 0.5*D1*sin(alpha2) - 0.5*D2*cos(alpha2), 
            ND2);
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha2) - 0.5*D2*sin(alpha2), 
            S2 - 0.5*D1*sin(alpha2) + 0.5*D2*cos(alpha2), 
            ND1);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha2) - 0.5*D2*sin(alpha2), 
            S2 + 0.5*D1*sin(alpha2) + 0.5*D2*cos(alpha2), 
            ND2);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha2) + 0.5*D2*sin(alpha2), 
            S2 + 0.5*D1*sin(alpha2) - 0.5*D2*cos(alpha2), 
            ND1);

  S1 = 0.45*L1;
  S2 = 0.75*L2;
  XgNewComponent();
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha2) + 0.5*D2*sin(alpha2), 
            S2 - 0.5*D1*sin(alpha2) - 0.5*D2*cos(alpha2), 
            ND2);
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha2) - 0.5*D2*sin(alpha2), 
            S2 - 0.5*D1*sin(alpha2) + 0.5*D2*cos(alpha2), 
            ND1);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha2) - 0.5*D2*sin(alpha2), 
            S2 + 0.5*D1*sin(alpha2) + 0.5*D2*cos(alpha2), 
            ND2);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha2) + 0.5*D2*sin(alpha2), 
            S2 + 0.5*D1*sin(alpha2) - 0.5*D2*cos(alpha2), 
            ND1);

  S1 = 0.75*L1;
  S2 = 0.75*L2;
  XgNewComponent();
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha2) + 0.5*D2*sin(alpha2), 
            S2 - 0.5*D1*sin(alpha2) - 0.5*D2*cos(alpha2), 
            ND2);
  XgAddEdge(20, 
            S1 - 0.5*D1*cos(alpha2) - 0.5*D2*sin(alpha2), 
            S2 - 0.5*D1*sin(alpha2) + 0.5*D2*cos(alpha2), 
            ND1);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha2) - 0.5*D2*sin(alpha2), 
            S2 + 0.5*D1*sin(alpha2) + 0.5*D2*cos(alpha2), 
            ND2);
  XgAddEdge(20, 
            S1 + 0.5*D1*cos(alpha2) + 0.5*D2*sin(alpha2), 
            S2 + 0.5*D1*sin(alpha2) - 0.5*D2*cos(alpha2), 
            ND1);

}
/*
 *   class xsquare_circ:
 */

class xsquare_circ: public Xgen {
  public:
  xsquare_circ(char *cfg_filename) : Xgen() {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xsquare_circ::XgReadData(FILE *f) {
  int N;
  if(!Get(f, &N)) XgError("Couldn't read N.");
  XgAddEdge(1, -24, -24, N);
  XgAddEdge(2, 24, -24, N);
  XgAddEdge(3, 24, 24, N);
  XgAddEdge(4, -24, 24, N);

  double Rcirc; int Ncirc;
  if(!Get(f, &Rcirc)) XgError("Couldn't read Rout.");
  if(!Get(f, &Ncirc)) XgError("Couldn't read Nout.");

  XgNewComponent();
  for(int i=0; i<Ncirc; i++) {
    XgAddEdge(5, Rcirc*cos(i*2*M_PI/Ncirc), -Rcirc*sin(i*2*M_PI/Ncirc));
  }

}


/* PROGRAM BODY*/

main(int argc, char *argv[]) {
  if(argc < 2) {
    printf("Missing application name.\n");
    exit(0);
  }
  if(!strcmp(argv[1], "xsquare")) {
    xsquare X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xhole")) {
    xhole X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xcirc")) {
    xcirc X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xgamm")) {
    xgamm X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xstep")) {
    xstep X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xlist")) {
    xlist X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xduese")) {
    xduese X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xduese2")) {
    xduese2 X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xnozzle0")) {
    xnozzle0 X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xnozzle")) {
    xnozzle X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xspir2d")) {
    xspir2d X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xsep2d")) {
    xsep2d X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xmunich")) {
    xmunich X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xsquare_circ")) {
    xsquare_circ X(argv[2]);
    XgMainLoop(&X, argc, argv);
  }
  printf("Unknown application name.\n");
}












