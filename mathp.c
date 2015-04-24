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
ArbPer (vector norm)		// arbitrary perpendicular
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
lv_prod (LINVEC a, LINVEC b)
{
  return a.x * b.x + a.y * b.y;
}

LINVEC
lv_ineq (LINVEC a, LINVEC b)
{
  LINVEC res;

  res.x = a.x - b.x;
  res.y = a.y - b.y;
  return res;
}

LINVEC
lv_sum (LINVEC a, LINVEC b)
{
  LINVEC res;

  res.x = a.x + b.x;
  res.y = a.y + b.y;
  return res;
}

LINVEC
lv_scal (double a, LINVEC b)
{
  LINVEC res;

  res.x = a * b.x;
  res.y = a * b.y;
  return res;
}

double
lv_abs2 (LINVEC v)
{
  return v.x * v.x + v.y * v.y;
}

LINVEC
lv_arb_per (LINVEC lv)
{
  LINVEC res;

  res.x = -lv.y;
  res.y = lv.x;
  return res;
}

LINVEC
lv_tri_per (LINVEC A, LINVEC B, LINVEC C, int *err)
{
  LINVEC ac, ab, ad;
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
lv_bez_lin_inters (LINVEC A, LINVEC B, LINVEC C, LINVEC la, LINVEC lb,
		   double *t1, double *t2)
{
  double va, vb, vc, vo, a, b, c;
  LINVEC N;

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
lv_u4plane (LINVEC A, LINVEC B, LINVEC P)
{
  LINVEC AP, AB;
  double sig, ap2, ab2;

  AP = lv_ineq (P, A);
  AB = lv_ineq (B, A);
  ap2 = lv_abs2 (AP);
  ab2 = lv_abs2 (AB);
  sig = lv_prod (AP, AB) / sqrt (ap2 * ab2);

  return sqrt (ap2 / ab2) * sig;
}

int
intersect_lin_bsect (LINVEC la, LINVEC lb, LINVEC ba, LINVEC bb, LINVEC bc,
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
intersect_bsect_lin (LINVEC ba, LINVEC bb, LINVEC bc, LINVEC la, LINVEC lb,
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
intersect_bsect_bsect (LINVEC ba, LINVEC bb, LINVEC bc, LINVEC ba2,
		       LINVEC bb2, LINVEC bc2, double *ts, int nst)
{
  double pos, stp, out[2];
  int size = 0, i, j;

  stp = 1. / (double) nst;
  pos = 0.;
  for (i = 0; i < nst; i++)
    {
      int nroots;
      LINVEC la, lb;

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
intersect_lin_lin (LINVEC a1, LINVEC b1, LINVEC a2, LINVEC b2, double *ts)
{
  double on, an, bn, z;
  LINVEC N;

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
lv_pinpol (LINVEC pt, LINVEC * per, int sper, LINVEC * os, int so, int pers)
{
  int i;

  for (i = 0; i < pers; i++)
    if (lv_prod (per[i * sper], lv_ineq (pt, os[i * so])) <= 0.)
      return 0;
  return 1;
}

int
lv_pinbezpol (LINVEC pt, BEZIERP * ABC, int ncor)
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

LINVEC
m_bezier_point (LINVEC a, LINVEC b, LINVEC c, double t)
{
  return lv_sum (lv_scal ((1. - t) * (1. - t), a),
		 lv_sum (lv_scal (2. * t * (1. - t), b), lv_scal (t * t, c)));
}

LINVEC
m_bezier_grad (LINVEC a, LINVEC b, LINVEC c, double t)
{
  return lv_sum (lv_scal ((2. * t - 2.), a),
		 lv_sum (lv_scal (2. * (1. - 2. * t), b),
			 lv_scal (2. * t, c)));
}

PrimBuf
m_bezier_cut (PrimBuf prb, LINVEC a, LINVEC b, LINVEC c, double t1, double t2)
{
  LINVEC oa, ob, oc;

  oa = m_bezier_point (a, b, c, t1);
  oc = m_bezier_point (a, b, c, t2);
  ob = m_bezier_point (a, b, c, (t2 + t1) / 2.);
  ob =
    lv_ineq (lv_scal (2., ob), lv_sum (lv_scal (.5, oa), lv_scal (0.5, oc)));
  return pri_sqr_bezier (prb, lv2prp (oa), lv2prp (ob), lv2prp (oc),
			 0x000000);
}

int
m_bez_lin_intersection (LINVEC A, LINVEC B, LINVEC C, LINVEC E, LINVEC F,
			double *t1, double *u1, double *t2, double *u2)
{
  double Z, a, b, c, X[2], U[2];
  int nroots;

  Z = (F.y - E.y) / (F.x - E.x);
  a = A.y - 2. * B.y + C.y - A.x * Z + 2. * B.x * Z - C.x * Z;
  b = 2. * (B.y - A.y + A.x * Z - B.x * Z);
  c = A.y - A.x * Z + E.x * Z - E.y;
  nroots = gsl_poly_solve_quadratic (a, b, c, X, X + 1);
  if (nroots == 0)
    return 0;
  U[0] =
    ((1. - X[0]) * (1. - X[0]) * A.x + 2. * X[0] * (1. - X[0]) * B.x +
     X[0] * X[0] * C.x - E.x) / (F.x - E.x);
  *t1 = X[0];
  *u1 = U[0];
  if (nroots == 2)
    {
      U[1] =
	((1. - X[1]) * (1. - X[1]) * A.x + 2. * X[1] * (1. - X[1]) * B.x +
	 X[1] * X[1] * C.x - E.x) / (F.x - E.x);
      *t2 = X[1];
      *u2 = U[1];
    }
  return nroots;
}
