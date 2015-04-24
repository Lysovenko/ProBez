/*
 *      (C) Serhii Lysovenko <lisovenko.s at the Gmail>
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 3 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */

#include <vmath.h>
#include <pgl.h>
#include <stdio.h>

typedef struct
{
  vector vp;
  double Z;
} Viewpoint;
typedef struct
{
  Viewpoint *view_pt;
  tensor *tens;
  int n_corners;
} ProjPar;
typedef struct
{
  vector a, b;
  char mode;
} LINE3D;

typedef struct
{
  double b, e;
} PARBE;
typedef struct
{
  double x, y;
} LINVEC;

/* elements of 3d graphics */
typedef struct
{
  vector n, o;
  double r, l;
  int figure0, figure1, type;
} Cylinder;
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
} Sphere;

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
} Mirage;
typedef struct
{
  Cylinder *cyls;
  Sphere *sphers;
  Mirage *mirages;
  int ncyls, nsphers, nmirages;
} Elements3D;
typedef struct
{
  Elements3D *kupas;
  Mirage *mirages;
  int nkupas, nmirages;
} KUPAS3D;
typedef struct
{
  KUPAS3D ks3d;
  int position;
} KUPOS;
/* projections */
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

/* methods */
vector *Poligon (vector O, vector norm, double rpol, int ndot);
double alph (LINVEC lv);
LINVEC Projection (Viewpoint, vector v);
PARBE *add_be (PARBE * arr, double b, double e, int *nbe);
/* math */
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
/* PROJECTION */
SQRAREA start_area (LINVEC st);
SQRAREA enlarge_area (SQRAREA in, LINVEC en);
P_CYLINDER proection_cylinder (Cylinder cyl, int ncor, Viewpoint VP);
P_SPHERE proection_sphere (Sphere sph, int ncor, Viewpoint VP,
			   Mirage * init_mirages, tensor * tens);
pr_point lv2prp (LINVEC v);

Elements3D interpret_3dvl (FILE * fp);
Elements3D interpret_3d_xml (char *docname);
KUPAS3D interpret_kupas3d_xml (char *docname);
PrimBuf plot_projection_all (PrimBuf pb, P_KUPA all);
P_KUPA project_all (const Elements3D * all, tensor tens, int ncor,
		    const Viewpoint * VP);
void mask_all (P_KUPA all);
void kupa3d_del (Elements3D k3d);
void kupas3d_del (KUPAS3D ks3d);
PrimBuf image_generator (const Elements3D * k3d, const Viewpoint * VP,
			 const tensor * tens);
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
