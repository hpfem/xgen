/*
 *   xmain.cpp
 */

# include <stdio.h>
# include <string.h>
# include <strings.h>

/** 
 **   simple sample applications
 **/

# include "xgen.h"
# include "disc.h"

/*
 *   class xsquare:
 *   unit square whose edges are subdivided into N intervals
 *   boundary does not contain circular arcs
 */

class xsquare: public Xgen {
  public:
  xsquare(char *cfg_filename, bool nogui, int nsteps, bool overlay) 
          : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
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
 *   rectangular domain containing a circular hole
 *   boundary contains circular arcs
 */

class xhole: public Xgen {
  public:
  xhole(char *cfg_filename, bool nogui, int nsteps, bool overlay) 
        : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xhole::XgReadData(FILE *f) {
  double length, height, center_x, center_y, radius;
  int subdiv_height; 
  if(!Get(f, &length)) XgError("Couldn't read domain length.");
  if(!Get(f, &height)) XgError("Couldn't read domain height.");
  if(!Get(f, &center_x)) XgError("Couldn't read circle center coordinate.");
  if(!Get(f, &center_y)) XgError("Couldn't read circle center coordinate.");
  if(!Get(f, &radius)) XgError("Couldn't read circle radius.");
  if(!Get(f, &subdiv_height)) XgError("Couldn't read subdivision.");

  int subdiv_length = (int) (length / height * subdiv_height + 0.5);
  double half_perimeter = M_PI*radius;
  int subdiv_half_circle = (int) (half_perimeter / height * subdiv_height + 0.5);

  // outer rectangle
  double alpha = 0;
  XgAddBoundarySegment(3, 0, 0, subdiv_length, alpha);
  XgAddBoundarySegment(2, length, 0, subdiv_height, alpha);
  XgAddBoundarySegment(4, length, height, subdiv_length, alpha);
  XgAddBoundarySegment(1, 0, height, subdiv_height, alpha);

  // closing the loop
  XgCreateNewBoundaryComponent();

  // circle as two segments
  alpha = -180;
  XgAddBoundarySegment(5, center_x, center_y - radius, subdiv_half_circle, alpha);
  XgAddBoundarySegment(5, center_x, center_y + radius, subdiv_half_circle, alpha);
}

/*
 *   class xcirc:
 *   domain spans the area between two concentric circles
 *   boundary contains circular arcs
 */

class xcirc: public Xgen {
  public:
  xcirc(char *cfg_filename, bool nogui, int nsteps, bool overlay) 
        : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xcirc::XgReadData(FILE *f) {
  double r_ext, r_int;
  int n_ext;
  if(!Get(f, &r_ext)) XgError("Couldn't read radius of exterior circle.");
  if(!Get(f, &n_ext)) XgError("Couldn't read division of exterior circle.");
  if(!Get(f, &r_int)) XgError("Couldn't read radius of interior circle.");
  int n_int = (int)(r_int / r_ext * n_ext + 0.5);

  // outer circle as two segments
  double alpha = 180.;
  XgAddBoundarySegment(1, 0, 0 - r_ext, n_ext / 2, alpha);
  XgAddBoundarySegment(1, 0, 0 + r_ext, n_ext / 2, alpha);

  // closing the loop
  XgCreateNewBoundaryComponent();

  // interior circle as two segments
  alpha = -180.;
  XgAddBoundarySegment(2, 0, 0 - r_int, n_int / 2, alpha);
  XgAddBoundarySegment(2, 0, 0 + r_int, n_int / 2, alpha);
}

/*
 *   class xgamm: 
 *   represents a rectangular channel with a circular hump on the bottom
 *   boundary contains circular arcs
 */

class xgamm: public Xgen {
  public:
  xgamm(char *cfg_filename, bool nogui, int nsteps, bool overlay) 
        : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xgamm::XgReadData(FILE *f) {
  double domain_height, length_front, length_hump, length_back, hump_central_angle;
  int inlet_division;
  if(!Get(f, &domain_height)) XgError("Couldn't read domain_height.");
  if(!Get(f, &inlet_division)) XgError("Couldn't read inlet_division.");
  if(!Get(f, &length_front)) XgError("Couldn't read length_front.");
  if(!Get(f, &length_hump)) XgError("Couldn't read length_hump.");
  if(!Get(f, &hump_central_angle)) XgError("Couldn't read hump_central_angle.");
  if(!Get(f, &length_back)) XgError("Couldn't read length_back.");

  double h = domain_height / inlet_division;
  int n_front = (int) (length_front / h + 0.5);
  int n_hump = (int) (length_hump / h + 0.5);
  int n_back = (int) (length_back / h + 0.5);

  // straight segment in front of hump
  double alpha = 0;
  XgAddBoundarySegment(3, 0, 0, n_front, alpha);
  // hump
  alpha = -hump_central_angle;
  XgAddBoundarySegment(3, length_front, 0, n_hump, alpha);
  // straight segment behind hump
  alpha = 0;
  XgAddBoundarySegment(3, length_front + length_hump, 0, n_back, alpha);
  // vertical edge on the right
  alpha = 0;
  XgAddBoundarySegment(2, length_front + length_hump + length_back, 0, inlet_division, alpha);
  // ceiling
  alpha = 0;
  XgAddBoundarySegment(4, length_front + length_hump + length_back, domain_height, 
                       n_front + n_back + n_back, alpha);
  // vertical edge on the left (inlet)
  alpha = 0;
  XgAddBoundarySegment(1, 0, domain_height, inlet_division, alpha);
 
}

/*
 *   class xstep:
 *   forward facing step
 *   only contains straight lines
 */

class xstep: public Xgen {
  public:
  xstep(char *cfg_filename, bool nogui, int nsteps, bool overlay) 
        : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

void xstep::XgReadData(FILE *f) {
  double H;
  int N1, N2, N3, N4;
  if(!Get(f, &H)) XgError("Couldn't read H.");
  if(!Get(f, &N1)) XgError("Couldn't read N1.");
  if(!Get(f, &N2)) XgError("Couldn't read N2.");
  if(!Get(f, &N3)) XgError("Couldn't read N3.");
  if(!Get(f, &N4)) XgError("Couldn't read N4.");

  // boundary is just one loop
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
 *   general polygonal domain with possibly curved edges
 *   see cfg/xlist_curved_6.cfg for an example of input file
 */

class xlist: public Xgen {
  public:
  xlist(char *cfg_filename, bool nogui, int nsteps, bool overlay) 
        : Xgen(nogui, nsteps, overlay) {XgInit(cfg_filename);}
  virtual void XgReadData(FILE *f);
};

// Boundary edges have the following format:
// <marker start_x, start_y, subdivision, angle>
// where angle = 0 for straight edges and angle != 0 for curved.
// New component is introduced with '='.
void xlist::XgReadData(FILE *f) {
  int marker, subdiv, end = 0;
  Point a, b, c;
  double angle;
  char test[20];

  while(Get(f, &marker)) {
    if(!Get(f, &a.x)) XgError("Couldn't read a starting point.");
    if(!Get(f, &a.y)) XgError("Couldn't read a starting point.");
    if(!Get(f, &subdiv)) XgError("Couldn't read a subdivision number.");
    if(!Get(f, &angle)) XgError("Couldn't read an angle.");

    // Closing a boundary loop.
    XgCreateNewBoundaryComponent();

    XgAddBoundarySegment(marker, a.x, a.y, subdiv, angle);
    c = a;
    while(end = !Get(f, test), (end || test[0] == '=') ? false : true) {
      marker = atoi(test);
      if(!Get(f, &b.x)) XgError("Couldn't read a starting point.");
      if(!Get(f, &b.y)) XgError("Couldn't read a starting point.");
      if(!Get(f, &subdiv)) XgError("Couldn't read a subdivision number.");
      if(!Get(f, &angle)) XgError("Couldn't read an angle.");
      XgAddBoundarySegment(marker, b.x, b.y, subdiv, angle); 
      a = b;
    } 
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
  printf("Unknown application name.\n");
}












