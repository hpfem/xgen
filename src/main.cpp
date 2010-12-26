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
  xsquare(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xsquare::XgReadData(FILE *f) {
  int N; 
  if(!Get(f, &N)) XgError("Couldn't read N.");
  double alpha = 0;
  XgAddBoundarySegment(3, 0, 0, N, alpha);
  XgAddBoundarySegment(2, 1, 0, N, alpha);
  XgAddBoundarySegment(4, 1, 1, N, alpha);
  XgAddBoundarySegment(1, 0, 1, N, alpha);
}

/*
 *   class xhole:
 */

class xhole: public Xgen {
  public:
  xhole(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xhole::XgReadData(FILE *f) {
  double H;
  if(!Get(f, &H)) XgError("Couldn't read H.");

  double R;
  if(!Get(f, &R)) XgError("Couldn't read R.");
  int N2 = (int)(2*M_PI*R + 0.5); 

  int N;
  double alpha = 0;
  if(!Get(f, &N)) XgError("Couldn't read N.");
  XgAddBoundarySegment(3, 0, 0, 2*N, alpha);
  XgAddBoundarySegment(2, 1.5*N*H, 0, N, alpha);
  XgAddBoundarySegment(4, 1.5*N*H, N*H, 2*N, alpha);
  XgAddBoundarySegment(1, 0, N*H, N, alpha);
  XgCreateNewBoundaryComponent();
  for(int i=0; i<N2; i++) {
    XgAddBoundarySegment(5, 0.5*H*N + R*H*cos(i*2*M_PI/N2), 0.5*H*N - R*H*sin(i*2*M_PI/N2), 1, alpha);
  }
}

/*
 *   class xcirc:
 */

class xcirc: public Xgen {
  public:
  xcirc(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xcirc::XgReadData(FILE *f) {
  double Rout; int Nout; 
  if(!Get(f, &Rout)) XgError("Couldn't read Rout.");
  if(!Get(f, &Nout)) XgError("Couldn't read Nout.");

  double alpha = 0;
  for(int i=0; i<Nout; i++) {
    XgAddBoundarySegment(1, Rout*cos(i*2*M_PI/Nout), Rout*sin(i*2*M_PI/Nout), 1, alpha);
  }

  double H = 2*M_PI*Rout/Nout;

  double Rins;
  if(!Get(f, &Rins)) XgError("Couldn't read Rins.");
  int Nins = (int)(2*M_PI*Rins/H + 0.5);

  XgCreateNewBoundaryComponent();
  for(int i=0; i<Nins; i++) {
    XgAddBoundarySegment(2, Rins*cos(i*2*M_PI/Nins), -Rins*sin(i*2*M_PI/Nins), 1, alpha);
  }
}

/*
 *   class xgamm:
 */

class xgamm: public Xgen {
  public:
  xgamm(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
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

  double alpha = 0;
  XgAddBoundarySegment(3, 0, 0, N2, alpha);
  for(int i=0; i<N3; i++) {
    double x = (i - 0.5*N3)*H;
    XgAddBoundarySegment(3, (N2+i)*H, sqrt(r*r - x*x) - k, 1, alpha);
  }
  XgAddBoundarySegment(3, (N2+N3)*H, 0, N4, alpha);
  XgAddBoundarySegment(2, (N2+N3+N4)*H, 0, N1, alpha);
  XgAddBoundarySegment(4, (N2+N3+N4)*H, N1*H, N2+N3+N4, alpha);
  XgAddBoundarySegment(1, 0, N1*H, N1, alpha);
}

/*
 *   class xstep:
 */

class xstep: public Xgen {
  public:
  xstep(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
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

  double alpha = 0;
  XgAddBoundarySegment(3, 0, 0, N2, alpha);
  XgAddBoundarySegment(3, N2*H, 0, N1-N4, alpha);
  XgAddBoundarySegment(3, N2*H, (N1-N4)*H, N3, alpha);
  XgAddBoundarySegment(2, (N2+N3)*H, (N1-N4)*H, N4, alpha);
  XgAddBoundarySegment(4, (N2+N3)*H, N1*H, N2+N3, alpha);
  XgAddBoundarySegment(1, 0, N1*H, N1, alpha);
}

/*
 *   class xlist:
 */

class xlist: public Xgen {
  public:
  xlist(char *cfg_filename, bool nogui, int nsteps, bool overlay) 
        : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

// Boundary edges have the following format:
// <marker start_x, start_y, subdivision, alpha>
// where alpha = 0 for straight edges and alpha != 0 for curved.
// New component is introduced with '='.
void xlist::XgReadData(FILE *f) {
  int marker, subdiv, end = 0;
  Point a, b, c;
  double alpha;
  char test[20];

  while(Get(f, &marker)) {
    if(!Get(f, &a.x)) XgError("Couldn't read a point.");
    if(!Get(f, &a.y)) XgError("Couldn't read a point.");
    if(!Get(f, &subdiv)) XgError("Couldn't read a subdivision number.");
    if(!Get(f, &alpha)) XgError("Couldn't read a boundary segment's angle.");
    XgCreateNewBoundaryComponent();
    XgAddBoundarySegment(marker, a.x, a.y, subdiv, alpha);
    c = a;
    while(end = !Get(f, test), (end || test[0] == '=') ? 0:1) {
      marker = atoi(test);
      if(!Get(f, &b.x)) XgError("Couldn't read a point.");
      if(!Get(f, &b.y)) XgError("Couldn't read a point.");
      if(!Get(f, &subdiv)) XgError("Couldn't read a subdivision number.");
      if(!Get(f, &alpha)) XgError("Couldn't read a boundary segment's angle.");
      XgAddBoundarySegment(marker, b.x, b.y, subdiv, alpha); 
      a=b;
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
  xnozzle(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
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

  double alpha = 0;
  XgAddBoundarySegment(4, x_left, 0, N, alpha);
  XgAddBoundarySegment(3, x_right, 0, N2, alpha);
  for(int i=N; i>0; i--) {
    double x;
    x = x_left + (double)i*(x_right - x_left)/N;
    XgAddBoundarySegment(4, x, Radius(x), 1, alpha);
  }
  XgAddBoundarySegment(1, x_left, Radius(x_left), N1, alpha);
}

class xspir2d: public Xgen {
  public:
  xspir2d(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
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
  double alpha = 0;
  double hx = (x_right - x_left)/Nx;
  double x;
  for(int i=0; i<Nx; i++) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x));
    XgAddBoundarySegment(3, x, 0.25*cos(M_PI*x), 1, alpha);
  }
  x = x_right;
  //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x));
  XgAddBoundarySegment(2, x_right, 0.25*cos(M_PI*x_right), Ny, alpha);
  for(int i=Nx; i>0; i--) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x));
    XgAddBoundarySegment(4, x, 0.25*cos(M_PI*x) + 0.5, 1, alpha);
  }
  x = x_left;
  //printf("starting point %g, %g\n", x, 0.25*cos(M_PI*x) + 0.5);
  XgAddBoundarySegment(1, x, 0.25*cos(M_PI*x) + 0.5, Ny, alpha);
}

double sep_low(double x) {
  return 0.2*x*cos(M_PI*x) + x/10;
}

double sep_up(double x) {
  return 0.2*cos(M_PI*x-0.2) + 0.4 + x/10;
}

class xsep2d: public Xgen {
  public:
  xsep2d(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
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
  double alpha = 0;
  for(int i=0; i<Nx; i++) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, sep_low(x));
    XgAddBoundarySegment(3, x, sep_low(x), 1, alpha);
  }
  x = x_right;
  //printf("corner x_left low %g, %g\n", x_left, sep_low(x_left));
  //printf("corner x_left up%g, %g\n", x_left, sep_up(x_left));
  //printf("corner x_right low %g, %g\n", x_right, sep_low(x_right));
  //printf("corner x_right up  %g, %g\n", x_right, sep_up(x_right));
  XgAddBoundarySegment(2, x_right, sep_low(x), Ny2, alpha);
  for(int i=Nx; i>0; i--) {
    x = x_left + i*hx;
    //printf("starting point %g, %g\n", x, sep_up(x));
    XgAddBoundarySegment(4, x, sep_up(x), 1, alpha);
  }
  x = x_left;
  //printf("starting point %g, %g\n", x, sep_up(x));
  XgAddBoundarySegment(1, x, sep_up(x), Ny1, alpha);
}

/*
 *   class xsquare_circ:
 */

class xsquare_circ: public Xgen {
  public:
  xsquare_circ(char *cfg_filename, bool nogui, int nsteps, bool overlay) : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xsquare_circ::XgReadData(FILE *f) {
  double A;  // horizontal length
  double B;  // vertical length
  double S1; // x-coordinate of circle midpoint
  double S2; // y-coordinate of circle midpoint
  double R;  // radius of circle
  int N;     // division of circle
  if(!Get(f, &A)) XgError("Couldn't read A (rectangle width).");
  if(!Get(f, &B)) XgError("Couldn't read B (rectangle height).");
  if(!Get(f, &S1)) XgError("Couldn't read S1 (circle center x-coordinate).");
  if(!Get(f, &S2)) XgError("Couldn't read S2 (circle center y-coordinate).");
  if(!Get(f, &R)) XgError("Couldn't read R (circle radous).");
  if(!Get(f, &N)) XgError("Couldn't read N (circle subdivision).");
  double h = 2*M_PI*R / N;
  int NA = (int) (A / h + 0.5);
  int NB = (int) (B / h + 0.5);
  double alpha = 0;
  XgAddBoundarySegment(1, 0, 0, NA, alpha);
  XgAddBoundarySegment(2, A, 0, NB, alpha);
  XgAddBoundarySegment(3, A, B, NA, alpha);
  XgAddBoundarySegment(4, 0, B, NB, alpha);

  XgCreateNewBoundaryComponent();
  for(int i = 0; i < N; i++) {
    XgAddBoundarySegment(5, S1 + R * cos(i*2*M_PI/N), 
		    S2 - R * sin(i*2*M_PI/N), 1, alpha);
  }

}


/* PROGRAM BODY*/

main(int argc, char *argv[]) {

  // Reading command-line parameters.
  if(argc < 3) {
    printf("Overview of command-line parameters:\n");
    printf("1. application (project) name such as \"xgamm\" - mandatory.\n");
    printf("2. text configuration file such as \"cfg/xgamm.cfg\" - mandatory.\n");
    printf("3. \"-nogui N\", where N is the number of relaxation steps - optional.\n");
    printf("4. \"-overlay\" regular overlay pattern of points will be used, as opposed\n");
    printf("   to the default random distribution - optional.\n");
    exit(0);
  }
  bool nogui = false;
  bool overlay = false;
  int nsteps = -1;
  if (argc > 3) {
    if (!strcmp(argv[3], "-nogui")) {
      nogui = true;
      if (argc <= 4) {
        printf("The parameter \"-nogui\" must be followed by the number of relaxation steps.\n");
        exit(0);
      }     
      nsteps = atoi(argv[4]);
      if (nsteps < 0) {
        printf("Number of relaxation steps must be nonnegative.\n");
        exit(0);
      }
      if (argc > 5) {
        if (!strcmp(argv[5], "-overlay")) {
          overlay = true;
        }
        else {
          printf("Invalid command-line parameter. Did you mean \"-overlay\"?\n");
          exit(0);
        }
      }
    }
    else {
      if (!strcmp(argv[3], "-overlay")) {
        overlay = true;
      }
      else {
        printf("Invalid command-line parameter. Did you mean \"-nogui\"?\n");
        exit(0);
      }
    }
  }

  if(!strcmp(argv[1], "xsquare")) {
    xsquare X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xhole")) {
    xhole X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xcirc")) {
    xcirc X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xgamm")) {
    xgamm X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xstep")) {
    xstep X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xlist")) {
    xlist X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xnozzle")) {
    xnozzle X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xspir2d")) {
    xspir2d X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xsep2d")) {
    xsep2d X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  if(!strcmp(argv[1], "xsquare_circ")) {
    xsquare_circ X(argv[2], nogui, nsteps, overlay);
    XgMainLoop(&X, argc, argv);
  }
  printf("Unknown application name.\n");
}












