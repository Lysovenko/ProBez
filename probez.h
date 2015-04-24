#include <vmath.h>
#include <pgl.h>
#include <stdio.h>

/*
 *  нехай вся модель розміщена так, що проектування відбувається на площину X Y
 *  всі елементи моделі мають від'ємну Z
 */
typedef struct
{
  vector vp;
  double Z;
} VIEWPOINT;
typedef struct
{
  vector a, b;
  char mode;
} LINE3D;

// -typedef struct {vector o; double r;char mode;}SPHERE;
// -typedef struct {int atoms,types,N_symp;MF size;}CFGHEADC;
typedef struct
{
  double b, e;
} PARBE;
typedef struct
{
  double x, y;
} LINVEC;

// /////
// об'єкти
// елменти 3d графіки
typedef struct
{
  vector n, o;
  double r, l;
  int figure0, figure1, type;
} CYLINDER;
typedef struct
{
  vector n, o;
  double r;
  int figure, fig_id;
} HOLE;
typedef struct
{
  vector o;
  double r;
  int nh, mir_id;
  HOLE *holes;
} SPHERE;

typedef struct
{
  vector a, b;
} LINE;
typedef struct
{
  vector a, b, c;
} BEZIER;
typedef struct
{
  LINE *lines;
  BEZIER *beziers;
  int nlines, nbeziers, id, isRotate;
} MIRAGE;
typedef struct
{
  CYLINDER *cyls;
  SPHERE *sphers;
  MIRAGE *mirages;
  int ncyls, nsphers, nmirages;
} KUPA3D;
typedef struct
{
  KUPA3D *kupas;
  MIRAGE *mirages;
  int nkupas, nmirages;
} KUPAS3D;
typedef struct
{
  KUPAS3D ks3d;
  int position;
} KUPOS;
// проекційне представлення
typedef struct
{
  double max_x, max_y, min_x, min_y;
} SQRAREA;
typedef struct
{
  LINVEC a, b;
  char vis;
  int nbe;
  PARBE *be;
} LINEP;
typedef struct
{
  LINVEC a, b, c;
  char vis;
  int nbe;
  PARBE *be;
} BEZIERP;
typedef struct
{
  LINEP *lines;
  BEZIERP *beziers;
  int nlines, nbeziers, id;
} P_MIRAGE;
typedef struct
{
  LINVEC a, b;
  double lv1, lv2;
  BEZIERP *b1, *b2;
  LINEP l1, l2;
  double lvis;
  SQRAREA sqa;
} P_CYLINDER;
typedef struct
{
  BEZIERP *bs;
  int fig_id, figure;
} P_HOLE;
typedef struct
{
  BEZIERP *bs;
  double lvis;
  P_HOLE *holes;
  int nholes;
  P_MIRAGE mir;
  SQRAREA sqa;
} P_SPHERE;
typedef struct
{
  P_CYLINDER *cyls;
  P_SPHERE *sphers;
  LINEP *lines;
  int ncyls, nsphers, nlines, ncor;
} P_KUPA;

// методи
vector *Poligon (vector O, vector norm, double rpol, int ndot);
double alph (LINVEC lv);
LINVEC Projection (VIEWPOINT, vector v);
PARBE *add_be (PARBE * arr, double b, double e, int *nbe);
// math
int m_bez_lin_intersection (LINVEC A, LINVEC B, LINVEC C, LINVEC E, LINVEC F,
			    double *t1, double *u1, double *t2, double *u2);
double lv_prod (LINVEC a, LINVEC b);
LINVEC lv_ineq (LINVEC a, LINVEC b);
LINVEC lv_sum (LINVEC a, LINVEC b);
LINVEC lv_scal (double a, LINVEC b);
int m_bez_plan_intersection (vector A, vector B, vector C, vector O, vector N,
			     double *t1, double *t2);
LINVEC lv_tri_per (LINVEC A, LINVEC B, LINVEC C, int *err);
int lv_pinpol (LINVEC pt, LINVEC * per, int sper, LINVEC * os, int so,
	       int pers);
int lv_pinbezpol (LINVEC pt, BEZIERP * ABC, int ncor);
int lv_bez_lin_inters (LINVEC A, LINVEC B, LINVEC C, LINVEC E, LINVEC F,
		       double *t1, double *t2);
int intersect_lin_lin (LINVEC a1, LINVEC b1, LINVEC a2, LINVEC b2,
		       double *ts);
int intersect_lin_bsect (LINVEC la, LINVEC lb, LINVEC ba, LINVEC bb,
			 LINVEC bc, double *ts);
int intersect_bsect_lin (LINVEC ba, LINVEC bb, LINVEC bc, LINVEC la,
			 LINVEC lb, double *ts);
LINVEC m_bezier_point (LINVEC a, LINVEC b, LINVEC c, double t);
PrimBuf m_bezier_cut (PrimBuf prb, LINVEC a, LINVEC b, LINVEC c, double t1,
		      double t2);
vector ArbPer (vector norm);
int intersect_bsect_bsect (LINVEC ba, LINVEC bb, LINVEC bc,
			   LINVEC ba2, LINVEC bb2, LINVEC bc2,
			   double *ts, int nst);
// PROJECTION
SQRAREA start_area (LINVEC st);
SQRAREA enlarge_area (SQRAREA in, LINVEC en);
P_CYLINDER proection_cylinder (CYLINDER cyl, int ncor, VIEWPOINT VP);
P_SPHERE proection_sphere (SPHERE sph, int ncor, VIEWPOINT VP,
			   MIRAGE * init_mirages);
// /////
pr_point lv2prp (LINVEC v);
///
KUPA3D interpret_3dvl (FILE * fp);
KUPA3D interpret_3d_xml (char *docname);
KUPAS3D interpret_kupas3d_xml (char *docname);
PrimBuf plot_projection_all (PrimBuf pb, P_KUPA all);
P_KUPA project_all (KUPA3D all, tensor tens, int ncor, VIEWPOINT VP);
void mask_all (P_KUPA all);
void kupa3d_del (KUPA3D k3d);
void kupas3d_del (KUPAS3D ks3d);
/*====================*/
void InitRequests ();
void *GetRequest (int n);
void SetRequest (int n, void *what);
void UnsetRequest (int n);

/*=======================*/
void projection_all_del (P_KUPA all);

/// FIGURES
#define FIG_CYLINDER 1
#define FIG_SPHERE 2
