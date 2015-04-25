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
} Line3D;

typedef struct
{
  double b, e;
} ParBE;
typedef struct
{
  double x, y;
} LinVec;

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
} Hole;
typedef struct
{
  vector o;
  double r;
  int nh, mir_id;
  Hole *holes;
} Sphere;

typedef struct
{
  vector a, b;
} Line;
typedef struct
{
  vector a, b, c;
} SBezier;
typedef struct
{
  Line *lines;
  SBezier *beziers;
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
} Sets3D;
typedef struct
{
  Sets3D ks3d;
  int position;
} KUPOS;
/* projections */
typedef struct
{
  double max_x, max_y, min_x, min_y;
} SQRAREA;
typedef struct
{
  LinVec a, b;
  char vis;
  int nbe;
  ParBE *be;
} LINEP;
typedef struct
{
  LinVec a, b, c;
  char vis;
  int nbe;
  ParBE *be;
} BEZIERP;
typedef struct
{
  LINEP *lines;
  BEZIERP *beziers;
  int nlines, nbeziers, id;
} P_MIRAGE;
typedef struct
{
  LinVec a, b;
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
double alph (LinVec lv);
LinVec Projection (Viewpoint, vector v);
ParBE *add_be (ParBE * arr, double b, double e, int *nbe);
/* math */
int m_bez_lin_intersection (LinVec A, LinVec B, LinVec C, LinVec E, LinVec F,
			    double *t1, double *u1, double *t2, double *u2);
double lv_prod (LinVec a, LinVec b);
LinVec lv_ineq (LinVec a, LinVec b);
LinVec lv_sum (LinVec a, LinVec b);
LinVec lv_scal (double a, LinVec b);
int m_bez_plan_intersection (vector A, vector B, vector C, vector O, vector N,
			     double *t1, double *t2);
LinVec lv_tri_per (LinVec A, LinVec B, LinVec C, int *err);
int lv_pinpol (LinVec pt, LinVec * per, int sper, LinVec * os, int so,
	       int pers);
int lv_pinbezpol (LinVec pt, BEZIERP * ABC, int ncor);
int lv_bez_lin_inters (LinVec A, LinVec B, LinVec C, LinVec E, LinVec F,
		       double *t1, double *t2);
int intersect_lin_lin (LinVec a1, LinVec b1, LinVec a2, LinVec b2,
		       double *ts);
int intersect_lin_bsect (LinVec la, LinVec lb, LinVec ba, LinVec bb,
			 LinVec bc, double *ts);
int intersect_bsect_lin (LinVec ba, LinVec bb, LinVec bc, LinVec la,
			 LinVec lb, double *ts);
LinVec m_bezier_point (LinVec a, LinVec b, LinVec c, double t);
PrimBuf m_bezier_cut (PrimBuf prb, LinVec a, LinVec b, LinVec c, double t1,
		      double t2);
vector ArbPer (vector norm);
int intersect_bsect_bsect (LinVec ba, LinVec bb, LinVec bc,
			   LinVec ba2, LinVec bb2, LinVec bc2,
			   double *ts, int nst);
/* PROJECTION */
SQRAREA start_area (LinVec st);
SQRAREA enlarge_area (SQRAREA in, LinVec en);
P_CYLINDER proection_cylinder (Cylinder cyl, int ncor, Viewpoint VP);
P_SPHERE proection_sphere (Sphere sph, int ncor, Viewpoint VP,
			   Mirage * init_mirages, tensor * tens);
pr_point lv2prp (LinVec v);

Elements3D interpret_3dvl (FILE * fp);
Elements3D interpret_3d_xml (char *docname);
Sets3D interpret_sets3d_xml (char *docname);
PrimBuf plot_projection_all (PrimBuf pb, P_KUPA all);
P_KUPA project_all (const Elements3D * all, tensor tens, int ncor,
		    const Viewpoint * VP);
void mask_all (P_KUPA all);
void kupa3d_del (Elements3D k3d);
void sets3d_del (Sets3D ks3d);
PrimBuf image_generator (const Elements3D * k3d, const Viewpoint * VP,
			 const tensor * tens);

/*=======================*/
void projection_all_del (P_KUPA all);

/// FIGURES
#define FIG_CYLINDER 1
#define FIG_SPHERE 2
