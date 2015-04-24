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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vmath.h>
#include "probez.h"

static BEZIERP mask_sph_sph_helper (BEZIERP bp, BEZIERP * ABC, int ncor,
				    int mode, int m2);



/* Help function for next only ;) */
static LINEP
mask_cyl_cyl_helper (LINEP lp, LinVec * pts, LinVec * prps)
{
  int i, j, nt = 0;
  double m, start, ts[4], t[2];
  LINEP res = lp;

  res = lp;
  for (i = 0; i < 4; i++)
    {
      LinVec a, b;
      int nroots;

      a = pts[i];
      b = pts[i < 3 ? i + 1 : 0];
      nroots = intersect_lin_lin (lp.a, lp.b, a, b, t);
      memcpy (ts + nt, t, nroots * sizeof (double));
      nt += nroots;
    }
  for (i = 0; i < nt; i++)
    for (j = i + 1; j < nt; j++)
      if (ts[i] > ts[j])
	{
	  m = ts[i];
	  ts[i] = ts[j];
	  ts[j] = m;
	}
  start = 0.;
  for (i = 0; i < nt; i++)
    {
      double med;

      med = (start + ts[i]) / 2.;
      if (lv_pinpol
	  (lv_sum (lp.a, lv_scal (med, lv_ineq (lp.b, lp.a))), prps, 1, pts,
	   1, 4))
	res.be = add_be (res.be, start, ts[i], &res.nbe);
      start = ts[i];
    }
  if (start != 1.)
    {
      double med;

      med = (start + 1.) / 2.;
      if (lv_pinpol
	  (lv_sum (lp.a, lv_scal (med, lv_ineq (lp.b, lp.a))), prps, 1, pts,
	   1, 4))
	{
	  if (start == 0.)
	    res.vis = 0;
	  else
	    res.be = add_be (res.be, start, 1., &res.nbe);
	}
    }

  return res;
}

static BEZIERP mask_sph_cyl_helper (BEZIERP bp, LinVec * pts, LinVec * prps);
void
mask_cyl_cyl (P_CYLINDER * who, P_CYLINDER whom, int ncor)
{
  if (who->lvis > whom.lvis)
    {
      int i, j, err;
      LinVec pts[4], perp[4], tri[3];

      pts[0] = whom.l1.a;
      pts[1] = whom.l2.a;
      pts[2] = whom.l2.b;
      pts[3] = whom.l1.b;
      for (i = 0; i < 4; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      int pos;

	      pos = i + j;
	      if (pos >= 4)
		pos -= 4;
	      tri[j] = pts[pos];
	    }
	  perp[i] = lv_tri_per (tri[0], tri[1], tri[2], &err);
	  if (err)
	    return;
	}
      if (who->l1.vis)
	who->l1 = mask_cyl_cyl_helper (who->l1, pts, perp);
      if (who->l2.vis)
	who->l2 = mask_cyl_cyl_helper (who->l2, pts, perp);
      for (i = 0; i < ncor; i++)
	{
	  if (who->b1 && who->b1[i].vis)
	    who->b1[i] = mask_sph_cyl_helper (who->b1[i], pts, perp);
	  if (who->b2 && who->b2[i].vis)
	    who->b2[i] = mask_sph_cyl_helper (who->b2[i], pts, perp);
	}
    }
}

static LINEP
mask_cyl_sph_helper (LINEP lp, BEZIERP * ABC, int ncor, int mode)
{
  LINEP res = lp;
  int i, rr = 0;
  double ts[2], tsr[2];

  for (i = 0; i < ncor; i++)
    {
      int nroots;

      nroots =
	intersect_lin_bsect (lp.a, lp.b, ABC[i].a, ABC[i].b, ABC[i].c, tsr);
      memcpy (ts + rr, tsr, nroots * sizeof (double));
      rr += nroots;
    }
  if (rr == 0)
    if (lv_pinbezpol (lp.a, ABC, ncor) && mode)
      res.vis = 0;
  if (rr == 1)
    {
      if (lv_pinbezpol (lp.a, ABC, ncor) && mode)
	res.be = add_be (res.be, 0., ts[0], &res.nbe);
      else
	res.be = add_be (res.be, ts[0], 1., &res.nbe);
    }
  if (rr == 2)
    {
      if (ts[0] < ts[1])
	{
	  if (mode)
	    res.be = add_be (res.be, ts[0], ts[1], &res.nbe);
	  else
	    {
	      res.be = add_be (res.be, 0., ts[0], &res.nbe);
	      res.be = add_be (res.be, ts[1], 1., &res.nbe);
	    }
	}
      else
	{
	  if (mode)
	    res.be = add_be (res.be, ts[1], ts[0], &res.nbe);
	  else
	    {
	      res.be = add_be (res.be, 0., ts[1], &res.nbe);
	      res.be = add_be (res.be, ts[0], 1., &res.nbe);
	    }
	}
    }
  return res;
}

void
mask_cyl_sph (P_CYLINDER * who, P_SPHERE whom, int ncor)
{
  if (who->lvis > whom.lvis)
    {
      int i;

      if (who->l1.vis)
	who->l1 = mask_cyl_sph_helper (who->l1, whom.bs, ncor, 1);
      if (who->l2.vis)
	who->l2 = mask_cyl_sph_helper (who->l2, whom.bs, ncor, 1);

      for (i = 0; i < ncor; i++)	// start accelerator
	if (who->b1
	    && (lv_pinbezpol (who->b1[i].a, whom.bs, ncor)
		|| lv_pinbezpol (who->b1[i].b, whom.bs, ncor)))
	  break;
      if (i < ncor)
	{			// end accelerator

	  for (i = 0; i < ncor; i++)
	    {
	      if (who->b1[i].vis)
		if (lv_pinbezpol (who->b1[i].a, whom.bs, ncor)
		    && lv_pinbezpol (who->b1[i].b, whom.bs, ncor)
		    && lv_pinbezpol (who->b1[i].c, whom.bs, ncor))
		  who->b1[i].vis = 0;
	    }
	  for (i = 0; i < ncor; i++)
	    who->b1[i] =
	      mask_sph_sph_helper (who->b1[i], whom.bs, ncor, 1, 1);
	}
      for (i = 0; i < ncor; i++)	// start accelerator
	if (who->b2
	    && (lv_pinbezpol (who->b2[i].a, whom.bs, ncor)
		|| lv_pinbezpol (who->b2[i].b, whom.bs, ncor)))
	  break;
      if (i < ncor)
	{			// end accelerator
	  for (i = 0; i < ncor; i++)
	    if (who->b2[i].vis)
	      if (lv_pinbezpol (who->b2[i].a, whom.bs, ncor)
		  && lv_pinbezpol (who->b2[i].b, whom.bs, ncor)
		  && lv_pinbezpol (who->b2[i].c, whom.bs, ncor))
		who->b2[i].vis = 0;
	  for (i = 0; i < ncor; i++)
	    if (who->b2[i].vis)
	      who->b2[i] =
		mask_sph_sph_helper (who->b2[i], whom.bs, ncor, 1, 1);
	}

    }
}

static BEZIERP
mask_sph_sph_helper (BEZIERP bp, BEZIERP * ABC, int ncor, int mode, int m2)
{
  BEZIERP res = bp;

  int i, j, nt = 0;

  double ts[20], tsr[8], m, start;

  for (i = 0; i < ncor; i++)
    if (m2 || ABC[i].vis)
      {
	int nroots;

	nroots =
	  intersect_bsect_bsect (bp.a, bp.b, bp.c, ABC[i].a, ABC[i].b,
				 ABC[i].c, tsr, 10);

	memcpy (ts + nt, tsr, nroots * sizeof (double));
	nt += nroots;
      }

  for (i = 0; i < nt; i++)
    for (j = i + 1; j < nt; j++)
      if (ts[i] > ts[j])
	{
	  m = ts[i];
	  ts[i] = ts[j];
	  ts[j] = m;
	}
  start = 0.;
  for (i = 0; i < nt; i++)
    {
      double med;

      med = (start + ts[i]) / 2.;
      if (lv_pinbezpol (m_bezier_point (bp.a, bp.b, bp.c, med), ABC, ncor)
	  && mode)

	res.be = add_be (res.be, start, ts[i], &res.nbe);
      start = ts[i];
    }
  if (start != 1.)
    {
      double med;

      med = (start + 1.) / 2.;
      if (lv_pinbezpol (m_bezier_point (bp.a, bp.b, bp.c, med), ABC, ncor)
	  && mode)
	{
	  if (start == 0.)
	    res.vis = 0;
	  else
	    res.be = add_be (res.be, start, 1., &res.nbe);
	}
    }

  return res;
}

void
mask_mirage_sph (P_MIRAGE who, P_SPHERE whom, int ncor)
{
  int i;

  for (i = 0; i < who.nlines; i++)
    if (who.lines[i].vis)
      {
	who.lines[i] = mask_cyl_sph_helper (who.lines[i], whom.bs, ncor, 1);
      }

  for (i = 0; i < who.nbeziers; i++)
    {
      who.beziers[i] =
	mask_sph_sph_helper (who.beziers[i], whom.bs, ncor, 1, 1);
    }
}

void
mask_sph_sph (P_SPHERE * who, P_SPHERE whom, int ncor)
{
  if (who->lvis > whom.lvis)
    {
      int i, j;
      {				// end accelerator
	for (i = 0; i < ncor; i++)
	  if (who->bs[i].vis)
	    if (lv_pinbezpol (who->bs[i].a, whom.bs, ncor)
		&& lv_pinbezpol (who->bs[i].b, whom.bs, ncor)
		&& lv_pinbezpol (who->bs[i].c, whom.bs, ncor))
	      who->bs[i].vis = 0;
	for (i = 0; i < ncor; i++)
	  if (who->bs[i].vis)
	    who->bs[i] =
	      mask_sph_sph_helper (who->bs[i], whom.bs, ncor, 1, 1);
      }

      for (i = 0; i < who->nholes; i++)
	for (j = 0; j < ncor; j++)
	  if (who->holes[i].bs[j].vis)
	    if (lv_pinbezpol (who->holes[i].bs[j].a, whom.bs, ncor)
		&& lv_pinbezpol (who->holes[i].bs[j].b, whom.bs, ncor)
		&& lv_pinbezpol (who->holes[i].bs[j].c, whom.bs, ncor))
	      who->holes[i].bs[j].vis = 0;
      for (i = 0; i < who->nholes; i++)
	for (j = 0; j < ncor; j++)
	  if (who->holes[i].bs[j].vis)
	    who->holes[i].bs[j] =
	      mask_sph_sph_helper (who->holes[i].bs[j], whom.bs, ncor, 1, 1);
      mask_mirage_sph (who->mir, whom, ncor);
    }
}

/* Help function for next only ;) */
static BEZIERP
mask_sph_cyl_helper (BEZIERP bp, LinVec * pts, LinVec * prps)
{
  int i, j, nt = 0, nroots;
  double ts[4], tsr[2], m, start /* , end */ ;
  BEZIERP res = bp;

  for (i = 0; i < 4; i++)
    {
      LinVec a, b;

      a = pts[i];
      b = pts[i < 3 ? i + 1 : 0];
      nroots = intersect_bsect_lin (bp.a, bp.b, bp.c, a, b, tsr);
      memcpy (ts + nt, tsr, nroots * sizeof (double));
      nt += nroots;
    }
  for (i = 0; i < nt; i++)
    for (j = i + 1; j < nt; j++)
      if (ts[i] > ts[j])
	{
	  m = ts[i];
	  ts[i] = ts[j];
	  ts[j] = m;
	}
  start = 0.;
  for (i = 0; i < nt; i++)
    {
      double med;

      med = (start + ts[i]) / 2.;
      if (lv_pinpol
	  (m_bezier_point (bp.a, bp.b, bp.c, med), prps, 1, pts, 1, 4))
	res.be = add_be (res.be, start, ts[i], &res.nbe);
      start = ts[i];
    }
  if (start != 1.)
    {
      double med;

      med = (start + 1.) / 2.;
      if (lv_pinpol
	  (m_bezier_point (bp.a, bp.b, bp.c, med), prps, 1, pts, 1, 4))
	{
	  if (start == 0.)
	    res.vis = 0;
	  else
	    res.be = add_be (res.be, start, 1., &res.nbe);
	}
    }

  return res;
}

P_MIRAGE
mask_mirage_cyl (P_MIRAGE who, P_CYLINDER whom, int ncor, LinVec * pts,
		 LinVec * perp)
{
  int i;

  P_MIRAGE res;

  res = who;
  for (i = 0; i < who.nlines; i++)
    if (res.lines[i].vis)
      {
	res.lines[i] = mask_cyl_cyl_helper (who.lines[i], pts, perp);
	if (whom.b1 != NULL)
	  res.lines[i] = mask_cyl_sph_helper (who.lines[i], whom.b1, ncor, 1);
	if (whom.b2 != NULL)
	  res.lines[i] = mask_cyl_sph_helper (who.lines[i], whom.b2, ncor, 1);
      }

  for (i = 0; i < who.nbeziers; i++)
    {
      res.beziers[i] = mask_sph_cyl_helper (who.beziers[i], pts, perp);
      if (whom.b1)
	res.beziers[i] =
	  mask_sph_sph_helper (who.beziers[i], whom.b1, ncor, 1, 1);
      if (whom.b2)
	res.beziers[i] =
	  mask_sph_sph_helper (who.beziers[i], whom.b2, ncor, 1, 1);
    }
  return res;
}

void
mask_mirage_holes (P_MIRAGE who, BEZIERP * holeline, int nholes, int ncor)
{
  int i, j;

  for (i = 0; i < nholes; i++)
    {
      for (j = 0; j < ncor; j++)
	if (holeline[i * ncor + j].vis)
	  break;
      if (j < ncor)
	{
	  int k;

	  for (k = 0; k < who.nlines; k++)
	    if (who.lines[k].vis)
	      who.lines[k] =
		mask_cyl_sph_helper (who.lines[k], holeline + (i * ncor),
				     ncor, 1);
	  for (k = 0; k < who.nbeziers; k++)
	    who.beziers[k] =
	      mask_sph_sph_helper (who.beziers[k], holeline + (i * ncor),
				   ncor, 1, 0);
	}
    }
}

void
mask_mirage_self (P_MIRAGE who, BEZIERP * bs, int ncor)
{

  int k;

  for (k = 0; k < who.nlines; k++)
    if (who.lines[k].vis)
      who.lines[k] = mask_cyl_sph_helper (who.lines[k], bs, ncor, 0);
  for (k = 0; k < who.nbeziers; k++)
    who.beziers[k] = mask_sph_sph_helper (who.beziers[k], bs, ncor, 0, 1);
}

void
mask_sph_cyl (P_SPHERE * who, P_CYLINDER whom, int ncor, int cyl_id)
{
  if (who->lvis > whom.lvis)
    {
      int i, j, err;
      LinVec pts[4], perp[4], tri[3];

      pts[0] = whom.l1.a;
      pts[1] = whom.l2.a;
      pts[2] = whom.l2.b;
      pts[3] = whom.l1.b;
      for (i = 0; i < 4; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      int pos;

	      pos = i + j;
	      if (pos >= 4)
		pos -= 4;
	      tri[j] = pts[pos];
	    }
	  perp[i] = lv_tri_per (tri[0], tri[1], tri[2], &err);
	  if (err)
	    return;
	}
      for (i = 0; i < ncor; i++)
	if (who->bs[i].vis)
	  who->bs[i] = mask_sph_cyl_helper (who->bs[i], pts, perp);
      for (i = 0; i < who->nholes; i++)
	if (who->holes[i].figure != FIG_CYLINDER
	    || who->holes[i].fig_id != cyl_id)
	  for (j = 0; j < ncor; j++)
	    if (who->holes[i].bs[j].vis)
	      who->holes[i].bs[j] =
		mask_sph_cyl_helper (who->holes[i].bs[j], pts, perp);
      who->mir = mask_mirage_cyl (who->mir, whom, ncor, pts, perp);
    }
}

void
mask_all (P_KUPA all)
{
  int i, j;

  for (i = 0; i < all.ncyls; i++)
    for (j = 0; j < all.ncyls; j++)
      if (i != j)
	mask_cyl_cyl (all.cyls + i, all.cyls[j], all.ncor);
  for (i = 0; i < all.ncyls; i++)
    for (j = 0; j < all.nsphers; j++)
      mask_cyl_sph (all.cyls + i, all.sphers[j], all.ncor);
  for (i = 0; i < all.nsphers; i++)
    for (j = 0; j < all.nsphers; j++)
      if (i != j)
	mask_sph_sph (all.sphers + i, all.sphers[j], all.ncor);
  for (i = 0; i < all.nsphers; i++)
    for (j = 0; j < all.ncyls; j++)
      mask_sph_cyl (all.sphers + i, all.cyls[j], all.ncor, j);
  for (i = 0; i < all.nsphers; i++)
    {
      mask_mirage_self (all.sphers[i].mir, all.sphers[i].bs, all.ncor);
      for (j = 0; j < all.sphers[i].nholes; j++)
	mask_mirage_holes (all.sphers[i].mir, all.sphers[i].holes[j].bs,
			   1, all.ncor);
    }
}
