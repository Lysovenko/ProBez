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
#define Allocator(x) assert(x)

CylinderP
proection_cylinder (Cylinder cyl, int ncor, Viewpoint VP)
{
  CylinderP res;

  vector pvp, pvp1, *polig, longn;

  double l, beta, alpha, orient, orient2;

  int i;
  if (cyl.type == 1)
    if (DotProd (VecIneq (VP.vp, cyl.o), cyl.n) < 0.)
      {
	memset (&res, 0, sizeof (res));
	return res;
      }
  if (cyl.figure0 != 1 && cyl.r > 0)
    Allocator (res.b1 = calloc (ncor, sizeof (SBezierP)));
  else
    res.b1 = NULL;
  if (cyl.figure1 != 1 && cyl.r > 0)
    Allocator (res.b2 = calloc (ncor, sizeof (SBezierP)));
  else
    res.b2 = NULL;
  for (i = 0; i < ncor; i++)
    {
      if (res.b2)
	res.b2[i].vis = 1;
      if (res.b1)
	res.b1[i].vis = 1;
    }
  if (res.b2 || res.b1)
    polig = poligon (cyl.o, cyl.n, cyl.r / cos (M_PI / (double) ncor), ncor);
  else
    polig = 0;
  longn = ProdScal (cyl.l, cyl.n);
  pvp =
    VecIneq (VP.vp,
	     ProdScal (DotProd (cyl.n, VecIneq (VP.vp, cyl.o)), cyl.n));
  orient = DotProd (cyl.n, VecIneq (VP.vp, cyl.o));
  orient2 = DotProd (cyl.n, VecIneq (VP.vp, VecIneq (cyl.o, longn)));
  pvp1 = VecIneq (pvp, cyl.o);
  l = sqrt (VecAbs2 (pvp1));
  pvp1 = ProdScal (1. / l, pvp1);
  res.lvis =
    sqrt (VecAbs2
	  (VecIneq (VP.vp, VecSum (ProdScal (-cyl.l / 2., cyl.n), cyl.o))));

  if (l > cyl.r)
    {
      vector n3, v1, v2, pv1, pv2;

      n3 = CrossProd (cyl.n, pvp1);
      beta = cyl.r / l;
      alpha = sqrt (1. - beta * beta);
      v1 = VecSum (ProdScal (beta, pvp1), ProdScal (alpha, n3));
      v2 =
	VecSum (ProdScal (beta, pvp1), ProdScal (alpha, ProdScal (-1., n3)));
      pv1 = VecIneq (v2, ProdScal (DotProd (v1, v2) / VecAbs2 (v1), v1));
      pv2 = VecIneq (v1, ProdScal (DotProd (v1, v2) / VecAbs2 (v2), v2));
      if (res.b1 || res.b2)
	for (i = 0; i < ncor; i++)
	  {
	    vector mp1, mp2, p;

	    int nroots, nroots2, case1, case2;

	    double t1, t2;

	    mp1 =
	      ProdScal (0.5,
			VecSum (polig[i], polig[(i > 0) ? i - 1 : ncor - 1]));
	    mp2 =
	      ProdScal (0.5,
			VecSum (polig[i],
				polig[(i < (ncor - 1)) ? i + 1 : 0]));
	    case1 = 1;
	    nroots =
	      m_bez_plan_intersection (mp1, polig[i], mp2, cyl.o, pv1, &t1,
				       &t2);
	    if (nroots == 2 && t1 < 0.)
	      t1 = t2;
	    if (nroots && t1 > 0. && t1 < 1.)
	      {
		p =
		  VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
			  VecSum (ProdScal (2. * (1. - t1) * t1, polig[i]),
				  ProdScal (t1 * t1, mp2)));
		if (DotProd (pv2, VecIneq (p, cyl.o)) > 0.)
		  {
		    case1 = 0;
		    if (DotProd (pv1, VecIneq (mp1, cyl.o)) < 0.)
		      {
			if (res.b1 && orient <= 0.)
			  res.b1[i].be =
			    add_be (res.b1[i].be, 0., t1, &res.b1[i].nbe);
			if (res.b2 && orient2 >= 0.)
			  res.b2[i].be =
			    add_be (res.b2[i].be, 0., t1, &res.b2[i].nbe);
		      }
		    else
		      {
			if (res.b1 && orient <= 0.)
			  res.b1[i].be =
			    add_be (res.b1[i].be, t1, 1., &res.b1[i].nbe);
			if (res.b2 && orient2 >= 0.)
			  res.b2[i].be =
			    add_be (res.b2[i].be, t1, 1., &res.b2[i].nbe);
		      }
		  }
	      }
	    case2 = 1;
	    nroots2 =
	      m_bez_plan_intersection (mp1, polig[i], mp2, cyl.o, pv2, &t1,
				       &t2);
	    if (nroots == 2 && t1 < 0.)
	      t1 = t2;
	    if (nroots2 && t1 > 0. && t1 < 1.)
	      {
		p =
		  VecSum (ProdScal ((1. - t1) * (1. - t1), mp1),
			  VecSum (ProdScal (2. * (1. - t1) * t1, polig[i]),
				  ProdScal (t1 * t1, mp2)));
		if (DotProd (pv1, VecIneq (p, cyl.o)) > 0.)
		  {
		    case2 = 0;
		    if (DotProd (pv2, VecIneq (mp1, cyl.o)) < 0.)
		      {
			if (res.b1 && orient <= 0.)
			  res.b1[i].be =
			    add_be (res.b1[i].be, 0., t1, &res.b1[i].nbe);
			if (res.b2 && orient2 >= 0.)
			  res.b2[i].be =
			    add_be (res.b2[i].be, 0., t1, &res.b2[i].nbe);
		      }
		    else
		      {
			if (res.b1 && orient <= 0.)
			  res.b1[i].be =
			    add_be (res.b1[i].be, t1, 1., &res.b1[i].nbe);
			if (res.b2 && orient2 >= 0.)
			  res.b2[i].be =
			    add_be (res.b2[i].be, t1, 1., &res.b2[i].nbe);
		      }
		  }
	      }
	    if (case1 && case2
		&& !(DotProd (pv2, VecIneq (mp1, cyl.o)) > 0.
		     && DotProd (pv1, VecIneq (mp1, cyl.o)) > 0.)
		&& !(DotProd (pv2, VecIneq (mp2, cyl.o)) > 0.
		     && DotProd (pv1, VecIneq (mp2, cyl.o)) > 0.))
	      {
		if (res.b1 && orient <= 0.)
		  res.b1[i].vis = 0;
		if (res.b2 && orient2 >= 0.)
		  res.b2[i].vis = 0;
	      }
	  }

      /* Lines */
      v1 = ProdScal (cyl.r, v1);
      v2 = ProdScal (cyl.r, v2);
      res.l1.a = project (VP, VecSum (cyl.o, v1));
      res.sqa = start_area (res.l1.a);
      res.l1.b = project (VP, VecSum (cyl.o, VecIneq (v1, longn)));
      res.sqa = enlarge_area (res.sqa, res.l1.b);
      res.l1.vis = 1;
      res.l1.be = 0;
      res.l1.nbe = 0;
      res.l2.a = project (VP, VecSum (cyl.o, v2));
      res.sqa = enlarge_area (res.sqa, res.l2.a);
      res.l2.b = project (VP, VecSum (cyl.o, VecIneq (v2, longn)));
      res.sqa = enlarge_area (res.sqa, res.l2.b);
      if (cyl.r > 0)
	res.l2.vis = 1;
      res.l2.be = 0;
      res.l2.nbe = 0;
    }
  else
    {
      res.l1.vis = 0;
      res.l1.nbe = 0;
      res.l2.vis = 0;
      res.l2.nbe = 0;
      if (res.b2 && orient > 0.)
	for (i = 0; i < ncor; i++)
	  res.b2[i].vis = 0;
      if (res.b1 && orient < 0.)
	for (i = 0; i < ncor; i++)
	  res.b1[i].vis = 0;
    }
  if (res.b1 || res.b2)
    for (i = 0; i < ncor; i++)
      {
	vector mp1, mp2;

	mp1 =
	  ProdScal (0.5,
		    VecSum (polig[i],
			    (i > 0) ? polig[i - 1] : polig[ncor - 1]));
	mp2 =
	  ProdScal (0.5,
		    VecSum (polig[i],
			    (i < ncor - 1) ? polig[i + 1] : polig[0]));
	if (res.b1)
	  {
	    res.b1[i].a = project (VP, mp1);
	    res.b1[i].b = project (VP, polig[i]);
	    res.b1[i].c = project (VP, mp2);
	    res.sqa = enlarge_area (res.sqa, res.b1[i].a);
	    res.sqa = enlarge_area (res.sqa, res.b1[i].b);
	    res.sqa = enlarge_area (res.sqa, res.b1[i].c);
	  }
	if (res.b2)
	  {
	    res.b2[i].a = project (VP, VecIneq (mp1, longn));
	    res.b2[i].b = project (VP, VecIneq (polig[i], longn));
	    res.b2[i].c = project (VP, VecIneq (mp2, longn));
	    res.sqa = enlarge_area (res.sqa, res.b2[i].a);
	    res.sqa = enlarge_area (res.sqa, res.b2[i].b);
	    res.sqa = enlarge_area (res.sqa, res.b2[i].c);
	  }

      }

  if (polig)
    free (polig);
  return res;
}
