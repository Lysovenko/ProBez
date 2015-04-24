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
#include "probez.h"
pr_point
lv2prp (LinVec v)
{
  pr_point res;

  res.x = v.x + 1.;
  res.y = v.y + 1.;
  return res;
}

int
intersection (ParBE a, ParBE b, ParBE * m)
{
  if (a.b < b.b)
    if (b.b > a.e)
      return 0;
  if (a.b > b.b && a.b > b.e)
    return 0;
  m->e = a.e > b.e ? a.e : b.e;
  m->b = a.b < b.b ? a.b : b.b;
  return 1;
}

int
par_short (ParBE * arrbe, int nbe)
{
  int i, j, k, n, havedo;

  n = nbe;
  do
    {
      havedo = 0;
      for (i = 0; i < n; i++)
	{
	  for (j = i + 1; j < n; j++)
	    if (!(arrbe[i].b > arrbe[j].e || arrbe[j].b > arrbe[i].e))
	      {
		if (arrbe[i].b > arrbe[j].b)
		  arrbe[i].b = arrbe[j].b;
		if (arrbe[i].e < arrbe[j].e)
		  arrbe[i].e = arrbe[j].e;
		for (k = j + 1; k < n; k++)
		  arrbe[k - 1] = arrbe[k];
		n--;
		havedo = 1;
	      }
	}
    }
  while (havedo);
  return n;

}

ParBE *
be_lin_inv (ParBE * arbe, int *np)
{
  int size, i, j, k, Np;
  ParBE *res, med;
  double b0;

  res = 0;
  size = 0;
  Np = *np;
  for (i = 0; i < Np; i++)
    {
      double m;

      m = arbe[i].b;
      k = i;
      for (j = i + 1; j < Np; j++)
	if (m > arbe[j].b)
	  {
	    k = j;
	    m = arbe[j].b;
	  }
      if (k > i)
	{
	  med = arbe[i];
	  arbe[i] = arbe[k];
	  arbe[k] = med;
	}
    }

  b0 = 0.;
  for (i = 0; i < Np; i++)
    {
      if (arbe[i].b > b0)
	{
	  res = realloc (res, (size + 1) * sizeof (ParBE));
	  res[size].b = b0;
	  res[size].e = arbe[i].b;
	  size++;
	}
      b0 = arbe[i].e;

    }

  if (b0 < 1.)
    {
      res = realloc (res, (size + 1) * sizeof (ParBE));
      res[size].b = b0;
      res[size].e = 1.;
      size++;
    }
  *np = size;
  free (arbe);
  return res;
}

PrimBuf
plot_bezierp_be (PrimBuf prb, BEZIERP bp)
{
  int i;

  PrimBuf res;

  res = prb;
  for (i = 0; i < bp.nbe; i++)
    res = m_bezier_cut (res, bp.a, bp.b, bp.c, bp.be[i].b, bp.be[i].e);
  return res;
}

PrimBuf
plot_bezierp (PrimBuf prb, BEZIERP bp)
{
  return pri_sqr_bezier (prb, lv2prp (bp.a), lv2prp (bp.b), lv2prp (bp.c),
			 0x000000);

}

PrimBuf
plot_linep_be (PrimBuf prb, LINEP lp)
{
  int i;
  PrimBuf res;

  res = prb;
  for (i = 0; i < lp.nbe; i++)
    res =
      pri_line (res,
		lv2prp (lv_sum
			(lv_scal (1. - lp.be[i].b, lp.a),
			 lv_scal (lp.be[i].b, lp.b))),
		lv2prp (lv_sum
			(lv_scal (1. - lp.be[i].e, lp.a),
			 lv_scal (lp.be[i].e, lp.b))), 0x000000);

  return res;
}

PrimBuf
plot_linep (PrimBuf prb, LINEP lp)
{
  return pri_line (prb, lv2prp (lp.a), lv2prp (lp.b), 0x000000);
}

PrimBuf
plot_proection_cylinder (PrimBuf pb, P_CYLINDER syl, int ncor)
{
  PrimBuf res;
  int i;

  res = pb;
  if (syl.b1 || syl.b2)
    for (i = 0; i < ncor; i++)
      {
	if (syl.b1 && syl.b1[i].vis)
	  {
	    if (syl.b1[i].nbe)
	      {
		syl.b1[i].nbe = par_short (syl.b1[i].be, syl.b1[i].nbe);
		syl.b1[i].be = be_lin_inv (syl.b1[i].be, &syl.b1[i].nbe);
		res = plot_bezierp_be (res, syl.b1[i]);
	      }
	    else
	      res = plot_bezierp (res, syl.b1[i]);
	  }
	if (syl.b2 && syl.b2[i].vis)
	  {
	    if (syl.b2[i].nbe)
	      {
		syl.b2[i].nbe = par_short (syl.b2[i].be, syl.b2[i].nbe);
		syl.b2[i].be = be_lin_inv (syl.b2[i].be, &syl.b2[i].nbe);
		res = plot_bezierp_be (res, syl.b2[i]);
	      }
	    else
	      res = plot_bezierp (res, syl.b2[i]);
	  }
      }
  if (syl.l1.vis)
    {
      if (syl.l1.nbe)
	{
	  syl.l1.nbe = par_short (syl.l1.be, syl.l1.nbe);
	  syl.l1.be = be_lin_inv (syl.l1.be, &syl.l1.nbe);
	  res = plot_linep_be (res, syl.l1);
	}
      else
	res = plot_linep (res, syl.l1);
    }
  if (syl.l2.vis)
    {
      if (syl.l2.nbe)
	{
	  syl.l2.nbe = par_short (syl.l2.be, syl.l2.nbe);
	  syl.l2.be = be_lin_inv (syl.l2.be, &syl.l2.nbe);
	  res = plot_linep_be (res, syl.l2);
	}
      else
	res = plot_linep (res, syl.l2);
    }
  return res;
}

PrimBuf
plot_proection_sphere (PrimBuf pb, P_SPHERE sph, int ncor)
{
  PrimBuf res;
  int i, j;

  res = pb;
  for (i = 0; i < ncor; i++)
    if (sph.bs[i].vis)
      {
	if (sph.bs[i].nbe)
	  {
	    sph.bs[i].nbe = par_short (sph.bs[i].be, sph.bs[i].nbe);
	    sph.bs[i].be = be_lin_inv (sph.bs[i].be, &sph.bs[i].nbe);
	    res = plot_bezierp_be (res, sph.bs[i]);
	  }
	else
	  res = plot_bezierp (res, sph.bs[i]);
      }
  for (i = 0; i < sph.nholes; i++)
    for (j = 0; j < ncor; j++)
      if (sph.holes[i].bs[j].vis)
	{
	  if (sph.holes[i].bs[j].nbe)
	    {
	      sph.holes[i].bs[j].nbe =
		par_short (sph.holes[i].bs[j].be, sph.holes[i].bs[j].nbe);
	      sph.holes[i].bs[j].be =
		be_lin_inv (sph.holes[i].bs[j].be, &sph.holes[i].bs[j].nbe);
	      res = plot_bezierp_be (res, sph.holes[i].bs[j]);
	    }
	  else
	    res = plot_bezierp (res, sph.holes[i].bs[j]);
	}
  // mirages can be nowhere except spheres
  for (i = 0; i < sph.mir.nlines; i++)
    if (sph.mir.lines[i].vis)
      {
	if (sph.mir.lines[i].nbe)
	  {
	    sph.mir.lines[i].nbe =
	      par_short (sph.mir.lines[i].be, sph.mir.lines[i].nbe);
	    sph.mir.lines[i].be =
	      be_lin_inv (sph.mir.lines[i].be, &sph.mir.lines[i].nbe);
	    res = plot_linep_be (res, sph.mir.lines[i]);
	  }
	else
	  res = plot_linep (res, sph.mir.lines[i]);
      }
  for (i = 0; i < sph.mir.nbeziers; i++)
    if (sph.mir.beziers[i].vis)
      {
	if (sph.mir.beziers[i].nbe)
	  {
	    sph.mir.beziers[i].nbe =
	      par_short (sph.mir.beziers[i].be, sph.mir.beziers[i].nbe);
	    sph.mir.beziers[i].be =
	      be_lin_inv (sph.mir.beziers[i].be, &sph.mir.beziers[i].nbe);
	    res = plot_bezierp_be (res, sph.mir.beziers[i]);
	  }
	else
	  res = plot_bezierp (res, sph.mir.beziers[i]);
      }
  return res;
}

PrimBuf
plot_projection_all (PrimBuf pb, P_KUPA all)
{
  PrimBuf res;
  int i;

  res = pb;
  for (i = 0; i < all.nsphers; i++)
    res = plot_proection_sphere (res, all.sphers[i], all.ncor);
  for (i = 0; i < all.ncyls; i++)
    res = plot_proection_cylinder (res, all.cyls[i], all.ncor);
  return res;
}
