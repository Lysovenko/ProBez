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

SqrArea
start_area (LinVec st)
{
  SqrArea res;

  res.max_x = res.min_x = st.x;
  res.max_y = res.min_y = st.y;
  return res;
}

SqrArea
enlarge_area (SqrArea in, LinVec en)
{
  SqrArea res;

  res.max_x = (in.max_x > en.x) ? in.max_x : en.x;
  res.max_y = (in.max_y > en.y) ? in.max_y : en.y;
  res.min_x = (in.min_x < en.x) ? in.min_x : en.x;
  res.min_y = (in.min_y < en.y) ? in.min_y : en.y;
  return res;
}

ParBE *
add_be (ParBE * arr, double b, double e, int *nbe)
{
  ParBE *res;
  int narr;

  narr = *nbe;
  Allocator (res = realloc (arr, (narr + 1) * sizeof (ParBE)));
  res[narr].b = b;
  res[narr].e = e;
  narr++;
  *nbe = narr;
  return res;
}

/*
 * Makes poligon arround point O in plane with normal vector 'norm',
 * radius 'rpol' and 'ndot' corners
 */
vector *
poligon (vector O, vector norm, double rpol, int ndot)
{
  vector *res;
  int i;
  double cosa, sina;

  Allocator (res = malloc (sizeof (vector) * ndot));
  res[0] = arb_perpendicular (norm);
  cosa = cos (2. * M_PI / (double) ndot);
  sina = sin (2. * M_PI / (double) ndot);
  for (i = 1; i < ndot; i++)
    {
      vector Y;

      Y = CrossProd (norm, res[i - 1]);
      res[i] = VecSum (ProdScal (cosa, res[i - 1]), ProdScal (sina, Y));
    }
  for (i = 0; i < ndot; i++)
    res[i] = VecSum (ProdScal (rpol, res[i]), O);
  return res;
}

/**projection to plane*/
LinVec
project (Viewpoint vp, vector v)
{
  LinVec res;
  double K;

  K = vp.Z / (vp.vp.z - v.z);
  res.x = (v.x - vp.vp.x) * K;
  res.y = (v.y - vp.vp.y) * K;
  return res;
}

SetP
project_all (const Elements3D * all, tensor tens, int ncor,
	     const Viewpoint * VP)
{
  int i;
  SetP res;

  res.ncor = ncor;
  res.ncyls = all->ncyls;
  Allocator (res.cyls = calloc (res.ncyls, sizeof (CylinderP)));
  for (i = 0; i < all->ncyls; i++)
    {
      Cylinder cyl;

      cyl = all->cyls[i];
      cyl.o = NewVector (cyl.o, tens);
      cyl.n = NewVector (cyl.n, tens);
      res.cyls[i] = proection_cylinder (cyl, ncor, *VP);
    }
  res.nsphers = all->nsphers;
  Allocator (res.sphers = calloc (res.nsphers, sizeof (SphereP)));
  for (i = 0; i < all->nsphers; i++)
    {
      Sphere sph;
      Hole *holes;
      int j;

      sph = all->sphers[i];
      Allocator (holes = malloc (sizeof (Hole) * sph.nh));
      memcpy (holes, sph.holes, sizeof (Hole) * sph.nh);
      sph.holes = holes;
      sph.o = NewVector (sph.o, tens);
      for (j = 0; j < sph.nh; j++)
	{
	  sph.holes[j].o = NewVector (sph.holes[j].o, tens);
	  sph.holes[j].n = NewVector (sph.holes[j].n, tens);
	}
      res.sphers[i] = proection_sphere (sph, ncor, *VP, all->mirages, &tens);
      free (holes);
    }
  return res;
}

void
bezierp_destruct (SBezierP b)
{
  if (b.nbe)
    free (b.be);
}

void
linep_destruct (LineP l)
{
  if (l.nbe)
    free (l.be);
}

void
proection_cylinder_del (CylinderP cyl, int ncor)
{
  int i;

  for (i = 0; i < ncor; i++)
    {
      if (cyl.b1)
	{
	  fprintf (stderr, "-");
	  bezierp_destruct (cyl.b1[i]);
	  free (cyl.b1);
	}
      if (cyl.b2)
	{
	  fprintf (stderr, "-");
	  bezierp_destruct (cyl.b2[i]);
	  free (cyl.b2);
	}
    }

}

void
proection_mirage_del (MirageP mir)
{
  int i;

  for (i = 0; i < mir.nlines; i++)
    linep_destruct (mir.lines[i]);
  if (mir.nlines)
    free (mir.lines);
  for (i = 0; i < mir.nbeziers; i++)
    bezierp_destruct (mir.beziers[i]);
  if (mir.nbeziers)
    free (mir.beziers);
}

void
projection_sphere_del (SphereP sph, int ncor)
{
  int i, j;

  for (i = 0; i < ncor; i++)
    bezierp_destruct (sph.bs[i]);
  free (sph.bs);

  for (i = 0; i < sph.nholes; i++)
    {
      for (j = 0; j < ncor; j++)
	bezierp_destruct (sph.holes[i].bs[j]);
      free (sph.holes[i].bs);
    }
  free (sph.holes);
  proection_mirage_del (sph.mir);
}

void
projection_all_del (SetP all)
{
  int i;

  if (all.ncyls)
    {
      free (all.cyls);
    }
  if (all.nsphers)
    {
      for (i = 0; i < all.nsphers; i++)
	projection_sphere_del (all.sphers[i], all.ncor);
      free (all.sphers);
    }
}

void
elements3d_del (Elements3D k3d)
{
  if (k3d.ncyls)
    free (k3d.cyls);
  if (k3d.nsphers)
    free (k3d.sphers);
  if (k3d.nmirages)
    free (k3d.mirages);
}

void
sets3d_del (Sets3D sets)
{
  int i;

  for (i = 0; i < sets.nsets; i++)
    {
      if (sets.mirages == sets.sets[i].mirages)
	sets.sets[i].nmirages = 0;
      elements3d_del (sets.sets[i]);
    }
  if (sets.nmirages)
    free (sets.mirages);
  if (sets.nsets)
    free (sets.sets);
}

PrimBuf
image_generator (Model * mod)
{
  PrimBuf prb = NULL;
  SetP pk;
  if (mod->set)
    {
      /* TODO: struct approximation replace "magic" numbers */
      pk = project_all (mod->set, *(mod->rot), mod->ncorners, mod->vp);
      mask_all (pk);
      prb = pri_init (2., 2.);
      prb = plot_projection_all (prb, pk);
      projection_all_del (pk);
    }
  return prb;
}
