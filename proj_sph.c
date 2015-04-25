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
#include "probez.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_poly.h>
#define Allocator(x) assert(x)

int
is_hole_vis (vector v2v, vector n, double R, double r)
{
  double sing, sinb, cosg, cosb, l, cosum;

  sing = r / R;
  l = sqrt (VecAbs2 (v2v));
  cosb = R / l;
  cosg = sqrt (1. - sing * sing);
  sinb = sqrt (1. - cosg * cosg);
  cosum = cosg * cosb - sing * sinb;
  return cosum < DotProd (v2v, n) / l;
}

int
is_hole_cut (vector v2v, vector n, double R, double r)
{
  double sing, sinb, cosg, cosb, l, cosum;

  sing = r / R;
  l = sqrt (VecAbs2 (v2v));
  cosb = R / l;
  cosg = sqrt (1. - sing * sing);
  sinb = sqrt (1. - cosg * cosg);
  cosum = cosg * cosb + sing * sinb;
  return cosum > DotProd (v2v, n) / l;
}

MirageP
proection_mirage (Mirage init_mir, Sphere sph, int sphID, Viewpoint VP,
		  tensor * tens)
{
  MirageP res;
  int i;

  res.id = sphID;
  res.nlines = init_mir.nlines;
  res.nbeziers = init_mir.nbeziers;
  if (res.nlines)
    res.lines = calloc (res.nlines, sizeof (LineP));
  if (res.nbeziers)
    res.beziers = calloc (res.nbeziers, sizeof (SBezierP));
  for (i = 0; i < res.nlines; i++)
    {
      res.lines[i].a = project (VP, VecSum (sph.o, init_mir.lines[i].a));
      res.lines[i].b = project (VP, VecSum (sph.o, init_mir.lines[i].b));
      res.lines[i].nbe = 0;
      res.lines[i].be = 0;
      res.lines[i].vis = 1;
    }
  for (i = 0; i < res.nbeziers; i++)
    {
      if (init_mir.isRotate)
	{
	  vector a, b, c, n;

	  double an, bn, cn;

	  a = NewVector (init_mir.beziers[i].a, *tens);
	  b = NewVector (init_mir.beziers[i].b, *tens);
	  c = NewVector (init_mir.beziers[i].c, *tens);
	  res.beziers[i].a = project (VP, VecSum (sph.o, a));
	  res.beziers[i].b = project (VP, VecSum (sph.o, b));
	  res.beziers[i].c = project (VP, VecSum (sph.o, c));
	  n = VecIneq (VP.vp, sph.o);
	  an = DotProd (n, a);
	  bn = DotProd (n, b);
	  cn = DotProd (n, c);

	  if (an >= 0. && bn >= 0. && cn >= 0.)
	    {
	      res.beziers[i].vis = 1;
	    }
	  else if (!(an <= 0. && bn <= 0. && cn <= 0.))
	    {
	      double A = an - 2. * bn + cn, B = 2. * (bn - an), C =
		an, X1, X2, st, en;

	      int nroots;

	      res.beziers[i].vis = 1;

	      nroots = gsl_poly_solve_quadratic (A, B, C, &X1, &X2);
	      if (nroots == 2 && X1 < 0.)
		X1 = X2;
	      if (an < 0)
		{
		  st = 0.;
		  en = X1;
		}
	      else
		{
		  st = X1;
		  en = 1.;
		}
	      res.beziers[i].be =
		add_be (res.beziers[i].be, st, en, &res.beziers[i].nbe);
	    }
	}
      else
	{
	  res.beziers[i].a =
	    project (VP, VecSum (sph.o, init_mir.beziers[i].a));
	  res.beziers[i].b =
	    project (VP, VecSum (sph.o, init_mir.beziers[i].b));
	  res.beziers[i].c =
	    project (VP, VecSum (sph.o, init_mir.beziers[i].c));
	  res.beziers[i].vis = 1;
	}
    }
  return res;
}

SphereP
proection_sphere (Sphere sph, int ncor, Viewpoint VP,
		  Mirage * init_mirages, tensor * tens)
{
  vector *polig, nv, O;
  SphereP res;
  int i, nholes;
  double R;

  nv = VecIneq (VP.vp, sph.o);
  memset (&res, 0, sizeof (SphereP));
  res.lvis = sqrt (VecAbs2 (nv));
  nv = ProdScal (1. / res.lvis, nv);
  R = sqrt (res.lvis * res.lvis - sph.r * sph.r) * sph.r / res.lvis;
  O = VecSum (sph.o, ProdScal (R * R / res.lvis, nv));
  polig = Poligon (O, nv, R / cos (M_PI / (double) ncor), ncor);
  Allocator (res.bs = calloc (ncor, sizeof (SBezierP)));
  for (i = 0; i < ncor; i++)
    {
      vector mp1, mp2;

      mp1 =
	ProdScal (0.5, VecSum (polig[i], polig[(i > 0) ? i - 1 : ncor - 1]));
      mp2 =
	ProdScal (0.5, VecSum (polig[i], polig[(i < ncor - 1) ? i + 1 : 0]));
      res.bs[i].a = project (VP, mp1);
      res.bs[i].b = project (VP, polig[i]);
      res.bs[i].c = project (VP, mp2);
      res.bs[i].vis = 1;
      if (i == 0)
	res.sqa = start_area (res.bs[i].a);
      else
	res.sqa = enlarge_area (res.sqa, res.bs[i].a);
      res.sqa = enlarge_area (res.sqa, res.bs[i].b);
      res.sqa = enlarge_area (res.sqa, res.bs[i].c);
    }
  nholes = 0;
  res.holes = 0;
  for (i = 0; i < sph.nh; i++)
    if (is_hole_vis
	(VecIneq (VP.vp, sph.o), sph.holes[i].n, sph.r, sph.holes[i].r))
      {
	vector *polig2;
	HoleP holelin;
	int j, cut =
	  is_hole_cut (VecIneq (VP.vp, sph.o), sph.holes[i].n, sph.r,
		       sph.holes[i].r);
	holelin.fig_id = sph.holes[i].fig_id;
	holelin.figure = sph.holes[i].figure;
	polig2 =
	  Poligon (sph.holes[i].o, sph.holes[i].n,
		   sph.holes[i].r / cos (M_PI / (double) ncor), ncor);
	polig2 = realloc (polig2, ncor * 2 * sizeof (*polig2));
	for (j = 0; j < ncor; j++)
	  polig2[j + ncor] = ProdScal (0.5,
				       VecSum (polig2[j],
					       polig2[(j <
						       (ncor - 1)) ? j +
						      1 : 0]));
	Allocator (res.holes =
		   realloc (res.holes, sizeof (*res.holes) * (nholes + 1)));
	Allocator (holelin.bs = calloc (sizeof (*holelin.bs), ncor));
	for (j = 0; j < ncor; j++)
	  {
	    vector mp1, mp2;
	    mp1 = polig2[ncor + (j > 0 ? j - 1 : ncor - 1)];
	    mp2 = polig2[ncor + j];
	    holelin.bs[j].a = project (VP, mp1);
	    holelin.bs[j].b = project (VP, polig2[j]);
	    holelin.bs[j].c = project (VP, mp2);
	    holelin.bs[j].vis = 1;
	  }

	if (cut)
	  {
	    vector tpi, nnp, v1, v2, pv1, pv2;
	    int err;
	    double l2p;

	    nnp = CrossProd (nv, sph.holes[i].n);
	    nnp = ProdScal (1. / sqrt (VecAbs2 (nnp)), nnp);
	    tpi =
	      TriPlane (nv, DotProd (nv, O), sph.holes[i].n,
			DotProd (sph.holes[i].n, sph.holes[i].o), nnp,
			DotProd (nnp, O), &err);
	    l2p = sqrt (R * R - VecAbs2 (VecIneq (tpi, O)));
	    v1 = VecIneq (VecSum (tpi, ProdScal (l2p, nnp)), O);
	    v2 = VecIneq (VecSum (tpi, ProdScal (-l2p, nnp)), O);
	    pv1 =
	      VecIneq (v2, ProdScal (DotProd (v1, v2) / VecAbs2 (v1), v1));
	    pv2 =
	      VecIneq (v1, ProdScal (DotProd (v1, v2) / VecAbs2 (v2), v2));

	    for (j = 0; j < ncor; j++)
	      {
		int nroots, nroots2, case1, case2;
		vector mp1, mp2, p;
		double t1, t2, t_1, t_2;

		mp1 =
		  ProdScal (0.5,
			    VecSum (polig[j],
				    polig[(j > 0) ? j - 1 : ncor - 1]));
		mp2 =
		  ProdScal (0.5,
			    VecSum (polig[j],
				    polig[(j < ncor - 1) ? j + 1 : 0]));

		case1 = 1;
		nroots =
		  m_bez_plan_intersection (mp1, polig[j], mp2, O, pv1, &t1,
					   &t2);
		if (nroots == 2 && t1 < 0.)
		  t1 = t2;

		if (nroots && t1 > 0. && t1 < 1.)
		  {
		    p =
		      VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
			      VecSum (ProdScal
				      (2. * (1. - t1) * t1, polig[j]),
				      ProdScal (t1 * t1, mp2)));
		    if (DotProd (pv2, VecIneq (p, O)) > 0.)
		      {
			case1 = 0;
			t_1 = t1;
		      }
		  }
		case2 = 1;
		nroots2 =
		  m_bez_plan_intersection (mp1, polig[j], mp2, O, pv2, &t1,
					   &t2);
		if (nroots2 == 2 && t1 < 0.)
		  t1 = t2;
		if (nroots2 && t1 > 0. && t1 < 1.)
		  {
		    p =
		      VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
			      VecSum (ProdScal
				      (2. * (1. - t1) * t1, polig[j]),
				      ProdScal (t1 * t1, mp2)));
		    if (DotProd (pv1, VecIneq (p, O)) > 0.)
		      {
			case2 = 0;
			t_2 = t1;
		      }
		  }
		if ((!case1) && (!case2))
		  {
		    if (t_1 < t_2)
		      res.bs[j].be =
			add_be (res.bs[j].be, t_1, t_2, &res.bs[j].nbe);
		    else
		      res.bs[j].be =
			add_be (res.bs[j].be, t_2, t_1, &res.bs[j].nbe);
		  }
		else
		  {
		    if (!case1)
		      {
			if (DotProd (pv1, VecIneq (mp1, O)) > 0.)
			  {
			    res.bs[j].be =
			      add_be (res.bs[j].be, 0., t_1, &res.bs[j].nbe);
			  }
			else
			  {
			    res.bs[j].be =
			      add_be (res.bs[j].be, t_1, 1., &res.bs[j].nbe);
			  }
		      }

		    if (!case2)
		      {
			if (DotProd (pv2, VecIneq (mp1, O)) > 0.)
			  {
			    res.bs[j].be =
			      add_be (res.bs[j].be, 0., t_2, &res.bs[j].nbe);
			  }
			else
			  {
			    res.bs[j].be =
			      add_be (res.bs[j].be, t_2, 1., &res.bs[j].nbe);
			  }
		      }
		  }

		if (case1 && case2
		    && (DotProd (pv2, VecIneq (mp1, O)) > 0.
			&& DotProd (pv1, VecIneq (mp1, O)) > 0.)
		    && (DotProd (pv2, VecIneq (mp2, O)) > 0.
			&& DotProd (pv1, VecIneq (mp2, O)) > 0.))
		  {
		    res.bs[j].vis = 0;
		  }

	      }
	    v1 = VecIneq (VecSum (tpi, ProdScal (l2p, nnp)), sph.holes[i].o);
	    v2 = VecIneq (VecSum (tpi, ProdScal (-l2p, nnp)), sph.holes[i].o);
	    pv1 =
	      VecIneq (v2, ProdScal (DotProd (v1, v2) / VecAbs2 (v1), v1));
	    pv2 =
	      VecIneq (v1, ProdScal (DotProd (v1, v2) / VecAbs2 (v2), v2));
	    if (DotProd (sph.holes[i].n, VecIneq (VP.vp, sph.holes[i].o)) <
		0.)
	      for (j = 0; j < ncor; j++)
		{
		  vector mp1, mp2, p;
		  int nroots, nroots2, case1, case2;
		  double t1, t2;

		  mp1 = polig2[ncor + (j > 0 ? j - 1 : ncor - 1)];
		  mp2 = polig2[ncor + j];

		  case1 = 1;
		  nroots =
		    m_bez_plan_intersection (mp1, polig2[j], mp2,
					     sph.holes[i].o, pv1, &t1, &t2);
		  if (nroots == 2 && t1 < 0.)
		    t1 = t2;
		  if (nroots && t1 > 0. && t1 < 1.)
		    {
		      p =
			VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
				VecSum (ProdScal
					(2. * (1. - t1) * t1, polig2[j]),
					ProdScal (t1 * t1, mp2)));
		      if (DotProd (pv2, VecIneq (p, sph.holes[i].o)) > 0.)
			{
			  case1 = 0;
			  if (DotProd (pv1, VecIneq (mp1, sph.holes[i].o)) <=
			      0.)
			    {
			      holelin.bs[j].be =
				add_be (holelin.bs[j].be, 0., t1,
					&holelin.bs[j].nbe);
			    }
			  else
			    {
			      holelin.bs[j].be =
				add_be (holelin.bs[j].be, t1, 1.,
					&holelin.bs[j].nbe);
			    }
			}
		    }
		  case2 = 1;
		  nroots2 =
		    m_bez_plan_intersection (mp1, polig2[j], mp2,
					     sph.holes[i].o, pv2, &t1, &t2);
		  if (nroots == 2 && t1 < 0.)
		    t1 = t2;
		  if (nroots2 && t1 > 0. && t1 < 1.)
		    {
		      p =
			VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
				VecSum (ProdScal
					(2. * (1. - t1) * t1, polig2[j]),
					ProdScal (t1 * t1, mp2)));
		      if (DotProd (pv1, VecIneq (p, sph.holes[i].o)) > 0.)
			{
			  case2 = 0;
			  if (DotProd (pv2, VecIneq (mp1, sph.holes[i].o)) <=
			      0.)
			    {
			      holelin.bs[j].be =
				add_be (holelin.bs[j].be, 0., t1,
					&holelin.bs[j].nbe);
			    }
			  else
			    {
			      holelin.bs[j].be =
				add_be (holelin.bs[j].be, t1, 1.,
					&holelin.bs[j].nbe);
			    }
			}
		    }
		  if (case1 && case2
		      && !(DotProd (pv2, VecIneq (mp1, sph.holes[i].o)) > 0.
			   && DotProd (pv1,
				       VecIneq (mp1, sph.holes[i].o)) > 0.)
		      && !(DotProd (pv2, VecIneq (mp2, sph.holes[i].o)) > 0.
			   && DotProd (pv1,
				       VecIneq (mp2, sph.holes[i].o)) > 0.))
		    {
		      holelin.bs[j].vis = 0;
		    }
		}

	  }
	// /////////////////////////
	if (sph.holes[i].figure == FIG_CYLINDER)
	  {
	    vector n3, v1, v2, pv1, pv2, pvp, pvp1;
	    double l, alpha, beta;

	    pvp =
	      VecIneq (VP.vp,
		       ProdScal (DotProd
				 (sph.holes[i].n,
				  VecIneq (VP.vp, sph.holes[i].o)),
				 sph.holes[i].n));
	    pvp1 = VecIneq (pvp, sph.holes[i].o);
	    l = sqrt (VecAbs2 (pvp1));
	    pvp1 = ProdScal (1. / l, pvp1);
	    if (l > sph.holes[i].r)
	      {
		n3 = CrossProd (sph.holes[i].n, pvp1);
		beta = sph.holes[i].r / l;
		alpha = sqrt (1. - beta * beta);
		v1 = VecSum (ProdScal (beta, pvp1), ProdScal (alpha, n3));
		v2 =
		  VecSum (ProdScal (beta, pvp1),
			  ProdScal (alpha, ProdScal (-1., n3)));
		pv1 =
		  VecIneq (v2,
			   ProdScal (DotProd (v1, v2) / VecAbs2 (v1), v1));
		pv2 =
		  VecIneq (v1,
			   ProdScal (DotProd (v1, v2) / VecAbs2 (v2), v2));
		for (j = 0; j < ncor; j++)
		  {
		    vector mp1, mp2, p;
		    int nroots, nroots2, case1, case2;
		    double t1, t2;

		    mp1 = polig2[ncor + (j > 0 ? j - 1 : ncor - 1)];
		    mp2 = polig2[ncor + j];

		    case1 = 1;
		    nroots =
		      m_bez_plan_intersection (mp1, polig2[j], mp2,
					       sph.holes[i].o, pv1, &t1, &t2);
		    if (nroots == 2 && t1 < 0.)
		      t1 = t2;
		    if (nroots && t1 > 0. && t1 < 1.)
		      {
			p =
			  VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
				  VecSum (ProdScal
					  (2. * (1. - t1) * t1, polig2[j]),
					  ProdScal (t1 * t1, mp2)));
			if (DotProd (pv2, VecIneq (p, sph.holes[i].o)) > 0.)
			  {
			    case1 = 0;
			    if (DotProd (pv1, VecIneq (mp1, sph.holes[i].o))
				<= 0.)
			      {
				holelin.bs[j].be =
				  add_be (holelin.bs[j].be, 0., t1,
					  &holelin.bs[j].nbe);
			      }
			    else
			      {
				holelin.bs[j].be =
				  add_be (holelin.bs[j].be, t1, 1.,
					  &holelin.bs[j].nbe);
			      }
			  }
		      }
		    case2 = 1;
		    nroots2 =
		      m_bez_plan_intersection (mp1, polig2[j], mp2,
					       sph.holes[i].o, pv2, &t1, &t2);
		    if (nroots == 2 && t1 < 0.)
		      t1 = t2;
		    if (nroots2 && t1 > 0. && t1 < 1.)
		      {
			p =
			  VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
				  VecSum (ProdScal
					  (2. * (1. - t1) * t1, polig2[j]),
					  ProdScal (t1 * t1, mp2)));
			if (DotProd (pv1, VecIneq (p, sph.holes[i].o)) > 0.)
			  {
			    case2 = 0;
			    if (DotProd (pv2, VecIneq (mp1, sph.holes[i].o))
				<= 0.)
			      {
				holelin.bs[j].be =
				  add_be (holelin.bs[j].be, 0., t1,
					  &holelin.bs[j].nbe);
			      }
			    else
			      {
				holelin.bs[j].be =
				  add_be (holelin.bs[j].be, t1, 1.,
					  &holelin.bs[j].nbe);
			      }
			  }
		      }
		    if (case1 && case2
			&& !(DotProd (pv2, VecIneq (mp1, sph.holes[i].o)) > 0.
			     && DotProd (pv1,
					 VecIneq (mp1, sph.holes[i].o)) > 0.)
			&& !(DotProd (pv2, VecIneq (mp2, sph.holes[i].o)) > 0.
			     && DotProd (pv1,
					 VecIneq (mp2, sph.holes[i].o)) > 0.))
		      {
			holelin.bs[j].vis = 0;
		      }
		  }

	      }
	    else
	      memset (holelin.bs, 0, sizeof (SBezierP) * ncor);
	  }
	else
	  for (j = 0; j < ncor; j++)
	    {
	      vector mp1, mp2;
	      int nroots;
	      double t1, t2;
	      mp1 = polig2[ncor + (j > 0 ? j - 1 : ncor - 1)];
	      mp2 = polig2[ncor + j];


	      if (DotProd (VecIneq (polig2[j], O), nv) < 0.
		  && DotProd (VecIneq (mp1, O), nv) < 0.
		  && DotProd (VecIneq (mp2, O), nv) < 0.)
		holelin.bs[j].vis = 0;
	      // place for accelerator
	      nroots =
		m_bez_plan_intersection (mp1, polig2[j], mp2, O, nv, &t1,
					 &t2);
	      if (nroots == 2 && t1 < 0.)
		{
		  t1 = t2;
		  nroots = 1;
		}
	      if (nroots == 2 && t1 == t2)
		nroots = 0;
	      if (nroots == 2 && t1 > 0. && t1 < 1. && t2 > 0. && t2 < 1.)
		{
		  if (DotProd (VecIneq (polig2[j], O), nv) < 0.)
		    {
		      holelin.bs[j].be =
			add_be (holelin.bs[j].be, t1, t2, &holelin.bs[j].nbe);
		    }
		  else
		    {
		      holelin.bs[j].be =
			add_be (holelin.bs[j].be, 0., t1, &holelin.bs[j].nbe);
		      holelin.bs[j].be =
			add_be (holelin.bs[j].be, t2, 1., &holelin.bs[j].nbe);
		    }
		}
	      if (nroots == 1 && t1 > 0. && t1 < 1.)
		{
		  if (DotProd (VecIneq (mp1, O), nv) < 0.)
		    holelin.bs[j].be =
		      add_be (holelin.bs[j].be, 0., t1, &holelin.bs[j].nbe);
		  else
		    holelin.bs[j].be =
		      add_be (holelin.bs[j].be, t1, 1., &holelin.bs[j].nbe);
		}
	    }
	res.holes[nholes] = holelin;
	nholes++;
	free (polig2);
      }
  free (polig);
  res.nholes = nholes;
  if (sph.mir_id)
    {
      res.mir =
	proection_mirage (init_mirages[sph.mir_id - 1], sph, 0, VP, tens);
    }
  return res;
}
