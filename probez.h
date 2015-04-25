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

#ifndef __PROBEZ_H__
#define __PROBEZ_H__

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
  Elements3D *set;
  Mirage *mirages;
  int ncorners, nmirages;
  tensor *rot;
  Viewpoint *vp;
} Model;
typedef struct
{
  Elements3D *sets;
  Mirage *mirages;
  int nsets, nmirages;
} Sets3D;
/* projections */
typedef struct
{
  double max_x, max_y, min_x, min_y;
} SqrArea;
typedef struct
{
  LinVec a, b;
  char vis;
  int nbe;
  ParBE *be;
} LineP;
typedef struct
{
  LinVec a, b, c;
  char vis;
  int nbe;
  ParBE *be;
} SBezierP;
typedef struct
{
  LineP *lines;
  SBezierP *beziers;
  int nlines, nbeziers, id;
} MirageP;
typedef struct
{
  LinVec a, b;
  double lv1, lv2;
  SBezierP *b1, *b2;
  LineP l1, l2;
  double lvis;
  SqrArea sqa;
} CylinderP;
typedef struct
{
  SBezierP *bs;
  int fig_id, figure;
} HoleP;
typedef struct
{
  SBezierP *bs;
  double lvis;
  HoleP *holes;
  int nholes;
  MirageP mir;
  SqrArea sqa;
} SphereP;
typedef struct
{
  CylinderP *cyls;
  SphereP *sphers;
  LineP *lines;
  int ncyls, nsphers, nlines, ncor;
} SetP;

/* methods */
vector *poligon (vector O, vector norm, double rpol, int ndot);
LinVec project (Viewpoint, vector v);
ParBE *add_be (ParBE * arr, double b, double e, int *nbe);
/* math */
double lv_prod (LinVec a, LinVec b);
LinVec lv_ineq (LinVec a, LinVec b);
LinVec lv_sum (LinVec a, LinVec b);
LinVec lv_scal (double a, LinVec b);
int m_bez_plan_intersection (vector A, vector B, vector C, vector O, vector N,
			     double *t1, double *t2);
LinVec lv_tri_per (LinVec A, LinVec B, LinVec C, int *err);
int lv_pinpol (LinVec pt, LinVec * per, int sper, LinVec * os, int so,
	       int pers);
int lv_pinbezpol (LinVec pt, SBezierP * ABC, int ncor);
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
vector arb_perpendicular (vector norm);
int intersect_bsect_bsect (LinVec ba, LinVec bb, LinVec bc,
			   LinVec ba2, LinVec bb2, LinVec bc2,
			   double *ts, int nst);
/* PROJECTION */
SqrArea start_area (LinVec st);
SqrArea enlarge_area (SqrArea in, LinVec en);
CylinderP proection_cylinder (Cylinder cyl, int ncor, Viewpoint VP);
SphereP proection_sphere (Sphere sph, int ncor, Viewpoint VP,
			  Mirage * init_mirages, tensor * tens);
pr_point lv2prp (LinVec v);

Elements3D interpret_3dvl (FILE * fp);
Elements3D interpret_3d_xml (char *docname);
Sets3D interpret_sets3d_xml (char *docname);
PrimBuf plot_projection_all (PrimBuf pb, SetP all);
SetP project_all (const Elements3D * all, tensor tens, int ncor,
		  const Viewpoint * VP);
void mask_all (SetP all);
void elements3d_del (Elements3D k3d);
void sets3d_del (Sets3D sets);
PrimBuf image_generator (Model * mod);

/*=======================*/
void projection_all_del (SetP all);

enum Figures
{
  FIG_CYLINDER = 1,
  FIG_SPHERE
};

#endif
