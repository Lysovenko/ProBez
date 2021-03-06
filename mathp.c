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
#include <gsl/gsl_poly.h>

vector
arb_perpendicular (vector norm)
{
  vector res;

  if (norm.x != 0.)
    {
      res.y = 1.;
      res.z = 1.;
      res.x = -(norm.y + norm.z) / norm.x;
      return ProdScal (1. / sqrt (VecAbs2 (res)), res);
    }
  if (norm.y != 0.)
    {
      res.x = 1.;
      res.z = 1.;
      res.x = -(norm.x + norm.z) / norm.y;
      return ProdScal (1. / sqrt (VecAbs2 (res)), res);
    }
  if (norm.z != 0.)
    {
      res.y = 1.;
      res.x = 1.;
      res.x = -(norm.y + norm.x) / norm.z;
      return ProdScal (1. / sqrt (VecAbs2 (res)), res);
    }
  return res;
}

double
lv_prod (LinVec a, LinVec b)
{
  return a.x * b.x + a.y * b.y;
}

LinVec
lv_ineq (LinVec a, LinVec b)
{
  LinVec res;

  res.x = a.x - b.x;
  res.y = a.y - b.y;
  return res;
}

LinVec
lv_sum (LinVec a, LinVec b)
{
  LinVec res;

  res.x = a.x + b.x;
  res.y = a.y + b.y;
  return res;
}

LinVec
lv_scal (double a, LinVec b)
{
  LinVec res;

  res.x = a * b.x;
  res.y = a * b.y;
  return res;
}

double
lv_abs2 (LinVec v)
{
  return v.x * v.x + v.y * v.y;
}

LinVec
lv_arb_per (LinVec lv)
{
  LinVec res;

  res.x = -lv.y;
  res.y = lv.x;
  return res;
}

LinVec
lv_tri_per (LinVec A, LinVec B, LinVec C, int *err)
{
  LinVec ac, ab, ad;
  double ab2;

  *err = 0;
  ac = lv_ineq (C, A);
  ab = lv_ineq (B, A);
  ab2 = lv_abs2 (ab);
  if (ab2 == 0.)
    *err = 1;
  ad = lv_scal (lv_prod (ab, ac) / ab2, ab);
  return lv_ineq (ac, ad);
}

int
lv_bez_lin_inters (LinVec A, LinVec B, LinVec C, LinVec la, LinVec lb,
		   double *t1, double *t2)
{
  double va, vb, vc, vo, a, b, c;
  LinVec N;

  N = lv_arb_per (lv_ineq (la, lb));
  va = lv_prod (A, N);
  vb = lv_prod (B, N);
  vc = lv_prod (C, N);
  vo = lv_prod (la, N);
  a = va - 2. * vb + vc;
  b = 2. * (vb - va);
  c = va - vo;
  return gsl_poly_solve_quadratic (a, b, c, t1, t2);
}

double
lv_u4plane (LinVec A, LinVec B, LinVec P)
{
  LinVec AP, AB;
  double sig, ap2, ab2;

  AP = lv_ineq (P, A);
  AB = lv_ineq (B, A);
  ap2 = lv_abs2 (AP);
  ab2 = lv_abs2 (AB);
  sig = lv_prod (AP, AB) / sqrt (ap2 * ab2);

  return sqrt (ap2 / ab2) * sig;
}

int
intersect_lin_bsect (LinVec la, LinVec lb, LinVec ba, LinVec bb, LinVec bc,
		     double *ts)
{
  int nroots, res;
  double t1, t2, u;

  res = 0;
  nroots = lv_bez_lin_inters (ba, bb, bc, la, lb, &t1, &t2);
  if (nroots == 2 && t1 == t2)
    nroots = 1;
  if (nroots == 2 && t1 < 0.)
    {
      nroots = 1;
      t1 = t2;
    }
  if (nroots == 2 && t2 >= 0. && t2 <= 1.)
    {
      u = lv_u4plane (la, lb, m_bezier_point (ba, bb, bc, t2));
      if (u >= 0. && u <= 1.)
	{
	  ts[res] = u;
	  res++;
	}
    }
  if (nroots && t1 >= 0. && t1 <= 1.)
    {
      u = lv_u4plane (la, lb, m_bezier_point (ba, bb, bc, t1));
      if (u >= 0. && u <= 1.)
	{
	  ts[res] = u;
	  res++;
	}
    }
  return res;
}

int
intersect_bsect_lin (LinVec ba, LinVec bb, LinVec bc, LinVec la, LinVec lb,
		     double *ts)
{
  int nroots, res;
  double t1, t2, u;

  res = 0;
  nroots = lv_bez_lin_inters (ba, bb, bc, la, lb, &t1, &t2);
  if (nroots == 2 && t1 == t2)
    nroots = 0;
  if (nroots == 2 && t1 < 0.)
    {
      nroots = 1;
      t1 = t2;
    }
  if (nroots == 2 && t2 >= 0. && t2 <= 1.)
    {
      u = lv_u4plane (la, lb, m_bezier_point (ba, bb, bc, t2));
      if (u >= 0. && u <= 1.)
	{
	  ts[res] = t2;
	  res++;
	}
    }
  if (nroots && t1 >= 0. && t1 <= 1.)
    {
      u = lv_u4plane (la, lb, m_bezier_point (ba, bb, bc, t1));
      if (u >= 0. && u <= 1.)
	{
	  ts[res] = t1;
	  res++;
	}
    }
  return res;
}

int
intersect_bsect_bsect (LinVec ba, LinVec bb, LinVec bc, LinVec ba2,
		       LinVec bb2, LinVec bc2, double *ts, int nst)
{
  double pos, stp, out[2];
  int size = 0, i, j;

  stp = 1. / (double) nst;
  pos = 0.;
  for (i = 0; i < nst; i++)
    {
      int nroots;
      LinVec la, lb;

      la = m_bezier_point (ba, bb, bc, pos);
      lb = m_bezier_point (ba, bb, bc, pos + stp);
      nroots = intersect_lin_bsect (la, lb, ba2, bb2, bc2, out);

      for (j = 0; j < nroots; j++)
	out[j] = pos + out[j] * stp;
      memcpy (ts + size, out, sizeof (double) * nroots);
      size += nroots;
      pos += stp;
    }
  return size;
}

int
intersect_lin_lin (LinVec a1, LinVec b1, LinVec a2, LinVec b2, double *ts)
{
  double on, an, bn, z;
  LinVec N;

  N = lv_ineq (b2, a2);
  N = lv_arb_per (N);
  on = lv_prod (N, a2);
  an = lv_prod (N, a1);
  bn = lv_prod (N, b1);
  z = bn - an;
  if (z == 0.)
    return 0;
  ts[0] = (on - an) / z;
  if (ts[0] < 0. || ts[0] > 1.)
    return 0;
  return 1;
}

int
lv_pinpol (LinVec pt, LinVec * per, int sper, LinVec * os, int so, int pers)
{
  int i;

  for (i = 0; i < pers; i++)
    if (lv_prod (per[i * sper], lv_ineq (pt, os[i * so])) <= 0.)
      return 0;
  return 1;
}

int
lv_pinbezpol (LinVec pt, SBezierP * ABC, int ncor)
{
  int i;
  double ts[2];

  for (i = 0; i < ncor; i++)
    {

      if (!intersect_lin_bsect
	  (pt, ABC[i].b, ABC[i].a, ABC[i].b, ABC[i].c, ts))
	return 0;
    }
  return 1;
}

int
m_bez_plan_intersection (vector A, vector B, vector C, vector O, vector N,
			 double *t1, double *t2)
{
  double va, vb, vc, vo, a, b, c;

  va = DotProd (A, N);
  vb = DotProd (B, N);
  vc = DotProd (C, N);
  vo = DotProd (O, N);
  a = va - 2. * vb + vc;
  b = 2. * (vb - va);
  c = va - vo;
  return gsl_poly_solve_quadratic (a, b, c, t1, t2);
}

LinVec
m_bezier_point (LinVec a, LinVec b, LinVec c, double t)
{
  return lv_sum (lv_scal ((1. - t) * (1. - t), a),
		 lv_sum (lv_scal (2. * t * (1. - t), b), lv_scal (t * t, c)));
}

LinVec
m_bezier_grad (LinVec a, LinVec b, LinVec c, double t)
{
  return lv_sum (lv_scal ((2. * t - 2.), a),
		 lv_sum (lv_scal (2. * (1. - 2. * t), b),
			 lv_scal (2. * t, c)));
}

PrimBuf
m_bezier_cut (PrimBuf prb, LinVec a, LinVec b, LinVec c, double t1, double t2)
{
  LinVec oa, ob, oc;

  oa = m_bezier_point (a, b, c, t1);
  oc = m_bezier_point (a, b, c, t2);
  ob = m_bezier_point (a, b, c, (t2 + t1) / 2.);
  ob =
    lv_ineq (lv_scal (2., ob), lv_sum (lv_scal (.5, oa), lv_scal (0.5, oc)));
  return pri_sqr_bezier (prb, lv2prp (oa), lv2prp (ob), lv2prp (oc));
}
