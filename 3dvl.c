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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include "probez.h"

Cylinder
interspher_cyl (Sphere * sph1, Sphere * sph2, double r, int cyl_id, int type)
{
  double r1, r2, l;
  vector nap;
  Cylinder res;
  res.type = type;
  r1 = sqrt (sph1->r * sph1->r - r * r);
  r2 = sqrt (sph2->r * sph2->r - r * r);
  nap = VecIneq (sph1->o, sph2->o);
  l = sqrt (VecAbs2 (nap));
  res.r = r;
  res.n = ProdScal (1. / l, nap);
  res.o = VecIneq (sph1->o, ProdScal (r1, res.n));
  res.l = l - r1 - r2;
  sph2->holes = realloc (sph2->holes, (sph2->nh + 1) * sizeof (Hole));
  sph2->holes[sph2->nh].n = res.n;
  sph2->holes[sph2->nh].o = VecSum (sph2->o, ProdScal (r2, res.n));
  sph2->holes[sph2->nh].r = r;
  sph2->holes[sph2->nh].figure = FIG_CYLINDER;
  sph2->holes[sph2->nh].fig_id = cyl_id;
  sph2->nh++;
  sph1->holes = realloc (sph1->holes, (sph1->nh + 1) * sizeof (Hole));
  sph1->holes[sph1->nh].n = ProdScal (-1., res.n);
  sph1->holes[sph1->nh].o = res.o;
  sph1->holes[sph1->nh].r = r;
  sph1->holes[sph1->nh].figure = FIG_CYLINDER;
  sph1->holes[sph1->nh].fig_id = cyl_id;
  sph1->nh++;
  res.figure1 = 1;
  res.figure0 = 1;
  return res;
}

typedef struct
{
  int n1, n2, type;
  double r;
} INTERSPHERCYL;

INTERSPHERCYL
parse_intersphcyl (xmlDocPtr doc, xmlNodePtr cur)
{
  INTERSPHERCYL res;
  xmlChar *key;
  int neibours = 0;
  res.type = 0;
  key = xmlGetProp (cur, (xmlChar *) "type");
  if (key)
    {
      if (!xmlStrcmp (key, (const xmlChar *) "Projection"))
	res.type = 1;
      xmlFree (key);
    }
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
    {
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "radius"))
	  || (!xmlStrcmp (cur->name, (const xmlChar *) "r")))
	{
	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  res.r = atof ((char *) key);
	  xmlFree (key);
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "neibour")))
	{
	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  if (neibours == 0)
	    {
	      res.n1 = atoi ((char *) key);
	      neibours++;
	    }
	  else
	    res.n2 = atoi ((char *) key);
	  xmlFree (key);
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "neibours"))
	  || (!xmlStrcmp (cur->name, (const xmlChar *) "nbs")))
	{
	  char *endptr, *sptr;
	  long int val;

	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  sptr = (char *) key;
	  val = strtol (sptr, &endptr, 10);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.n1 = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, "error reading neibours numbers\n");
	  val = strtol (sptr, &endptr, 10);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.n2 = val;
	    }
	  else
	    fprintf (stderr, "error reading neibours numbers\n");
	  xmlFree (key);
	}
      cur = cur->next;
    }
  return res;
}

Sphere
parse_sphere (xmlDocPtr doc, xmlNodePtr cur)
{

  xmlChar *key;

  cur = cur->xmlChildrenNode;
  Sphere res;

  memset (&res, 0, sizeof (Sphere));
  while (cur != NULL)
    {
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "center"))
	  || (!xmlStrcmp (cur->name, (const xmlChar *) "c")))
	{
	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  if (!sscanf
	      ((char *) key, "%lf %lf %lf", &res.o.x, &res.o.y,
	       &res.o.z) == 3)
	    fprintf (stderr, "error reading atom coordinates\n");
	  xmlFree (key);
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "radius"))
	  || (!xmlStrcmp (cur->name, (const xmlChar *) "r")))
	{
	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  res.r = atof ((char *) key);
	  xmlFree (key);
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "mir")))
	{
	  char *endptr, *sptr;
	  int val;

	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  sptr = (char *) key;
	  val = strtol (sptr, &endptr, 10);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.mir_id = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, "%s:%d\terror reading mirage number \n",
		     __FILE__, __LINE__);
	  xmlFree (key);
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "hole")))
	{
	  vector n;
	  double r;

	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  if (sscanf ((char *) key, "%lf %lf %lf %lf", &n.x, &n.y, &n.z, &r)
	      == 4)
	    {
	      double r1 = sqrt (res.r * res.r - r * r);

	      res.holes =
		realloc (res.holes, (res.nh + 1) * sizeof (*res.holes));
	      n = ProdScal (1. / sqrt (VecAbs2 (n)), n);
	      res.holes[res.nh].o = VecSum (res.o, ProdScal (r1, n));
	      res.holes[res.nh].n = n;
	      res.holes[res.nh].r = r;
	      res.nh++;
	    }
	  xmlFree (key);
	}

      cur = cur->next;
    }

  return res;
}

Mirage
parse_mirage (xmlDocPtr doc, xmlNodePtr cur)
{

  xmlChar *key;

  Mirage res;

  memset (&res, 0, sizeof (res));
  key = xmlGetProp (cur, (xmlChar *) "rotate");
  if (key)
    {
      if (!xmlStrcmp (key, (const xmlChar *) "True"))
	res.isRotate = 1;
      xmlFree (key);
    }
  cur = cur->xmlChildrenNode;

  while (cur != NULL)
    {
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "line")))
	{
	  char *endptr, *sptr, *errmsg = "error reading line coordinates\n";

	  double val;

	  res.lines =
	    realloc (res.lines, sizeof (*res.lines) * (res.nlines + 1));
	  memset (res.lines + res.nlines, 0, sizeof (*res.lines));
	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);

	  sptr = (char *) key;
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.lines[res.nlines].a.x = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.lines[res.nlines].a.y = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.lines[res.nlines].a.z = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.lines[res.nlines].b.x = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.lines[res.nlines].b.y = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.lines[res.nlines].b.z = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  res.nlines++;
	  xmlFree (key);
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "radius"))
	  || (!xmlStrcmp (cur->name, (const xmlChar *) "r")))
	{
	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);
	  res.beziers[res.nbeziers].a.x = atof ((char *) key);
	  xmlFree (key);
	}
      if (!xmlStrcmp (cur->name, (const xmlChar *) "bezier"))
	{
	  char *endptr, *sptr, *errmsg = "error reading bezier coordinates\n";

	  double val;

	  res.beziers =
	    realloc (res.beziers, sizeof (*res.beziers) * (res.nbeziers + 1));
	  memset (res.beziers + res.nbeziers, 0, sizeof (*res.beziers));
	  key = xmlNodeListGetString (doc, cur->xmlChildrenNode, 1);

	  sptr = (char *) key;
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].a.x = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].a.y = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].a.z = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].b.x = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].b.y = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].b.z = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].c.x = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].c.y = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  val = strtod (sptr, &endptr);
	  if (sptr != endptr && *sptr != 0)
	    {
	      res.beziers[res.nbeziers].c.z = val;
	      sptr = endptr;
	    }
	  else
	    fprintf (stderr, errmsg);
	  res.nbeziers++;
	  xmlFree (key);
	}
      cur = cur->next;
    }

  return res;
}

Elements3D
parse_kupa3d (xmlDocPtr doc, xmlNodePtr cur)
{
  Elements3D res;

  cur = cur->xmlChildrenNode;
  memset (&res, 0, sizeof (res));
  while (cur != NULL)
    {
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "sphere")))
	{
	  Sphere sph = parse_sphere (doc, cur);

	  res.sphers =
	    realloc (res.sphers, sizeof (Sphere) * (res.nsphers + 1));
	  res.sphers[res.nsphers] = sph;
	  res.nsphers++;
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "intersphcyl")))
	{
	  INTERSPHERCYL isph = parse_intersphcyl (doc, cur);

	  if (isph.n1 >= res.nsphers || isph.n2 >= res.nsphers)
	    {
	      fprintf (stderr, __FILE__ ":%d: Spheres Overflow %d %d %d \n",
		       __LINE__, isph.n1, isph.n2, res.nsphers);
	      continue;
	    }
	  res.cyls = realloc (res.cyls, sizeof (Cylinder) * (res.ncyls + 1));
	  res.cyls[res.ncyls] =
	    interspher_cyl (res.sphers + isph.n1, res.sphers + isph.n2,
			    isph.r, res.ncyls, isph.type);
	  res.ncyls++;
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "mirage")))
	{
	  Mirage mir = parse_mirage (doc, cur);

	  res.mirages =
	    realloc (res.mirages, sizeof (MirageP) * (res.nmirages + 1));
	  res.mirages[res.nmirages] = mir;
	  res.nmirages++;
	}

      cur = cur->next;
    }
  return res;
}

Sets3D
interpret_sets3d_xml (char *docname)
{
  xmlDocPtr doc;
  xmlNodePtr cur;
  Sets3D res;

  memset (&res, 0, sizeof (res));
  doc = xmlParseFile (docname);
  if (doc == NULL)
    {
      fprintf (stderr, "Document not parsed successfully. \n");
      return res;
    }
  cur = xmlDocGetRootElement (doc);
  if (cur == NULL)
    {
      fprintf (stderr, "empty document\n");
      xmlFreeDoc (doc);
      return res;
    }
  if (xmlStrcmp (cur->name, (const xmlChar *) "multibox"))
    {
      fprintf (stderr, "document of the wrong type, root node != multibox\n");
      xmlFreeDoc (doc);
      return res;
    }
  cur = cur->xmlChildrenNode;
  setlocale (LC_NUMERIC, "C");
  while (cur != NULL)
    {
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "box")))
	{
	  Elements3D k3d = parse_kupa3d (doc, cur);

	  res.sets =
	    realloc (res.sets, sizeof (*res.sets) * (res.nsets + 1));
	  res.sets[res.nsets] = k3d;
	  res.nsets++;
	}
      if ((!xmlStrcmp (cur->name, (const xmlChar *) "mirage")))
	{
	  Mirage mir = parse_mirage (doc, cur);

	  res.mirages =
	    realloc (res.mirages, sizeof (MirageP) * (res.nmirages + 1));
	  res.mirages[res.nmirages] = mir;
	  res.nmirages++;
	}
      cur = cur->next;
    }

  setlocale (LC_NUMERIC, "");
  xmlFreeDoc (doc);
  if (res.nmirages)
    {
      int i;

      for (i = 0; i < res.nsets; i++)
	{
	  if (res.sets[i].nmirages == 0)
	    {
	      res.sets[i].nmirages = res.nmirages;
	      res.sets[i].mirages = res.mirages;
	    }
	}
    }
  return res;
}

Elements3D
interpret_3d_xml (char *docname)
{
  xmlDocPtr doc;
  xmlNodePtr cur;
  Elements3D res;
  int boxmode = 0, multiboxmode = 0;

  setlocale (LC_NUMERIC, "C");
  doc = xmlParseFile (docname);
  memset (&res, 0, sizeof (Elements3D));
  if (doc == NULL)
    {
      fprintf (stderr, "Document not parsed successfully. \n");
      return res;
    }

  cur = xmlDocGetRootElement (doc);
  if (cur == NULL)
    {
      fprintf (stderr, "empty document\n");
      xmlFreeDoc (doc);
      return res;
    }
  if (!xmlStrcmp (cur->name, (const xmlChar *) "box"))
    boxmode = 1;
  if (!xmlStrcmp (cur->name, (const xmlChar *) "multibox"))
    multiboxmode = 1;
  if (!(boxmode || multiboxmode))
    {
      fprintf (stderr, "document of the wrong type, root node != box");
      xmlFreeDoc (doc);
      return res;
    }
  res = parse_kupa3d (doc, cur);
  xmlFreeDoc (doc);
  setlocale (LC_NUMERIC, "");
  return res;
}
