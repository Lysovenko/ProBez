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

#include <gtk/gtk.h>
#include "probez.h"
#include <pgl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <locale.h>
#include "interf.h"

void
dialog_wiewpoint (GtkButton * button, gpointer user_data)
{
  GtkWidget *dialog, *hbox, *stock, *table, *entry_x,
    *entry_H, *entry_y, *entry_z, *label, *entry_c;
  gint response;

  dialog =
    gtk_dialog_new_with_buttons ("View point", GTK_WINDOW (user_data),
				 GTK_DIALOG_MODAL |
				 GTK_DIALOG_DESTROY_WITH_PARENT, GTK_STOCK_OK,
				 GTK_RESPONSE_OK, GTK_STOCK_CANCEL,
				 GTK_RESPONSE_CANCEL, NULL);

  hbox = gtk_hbox_new (FALSE, 8);
  Viewpoint *vp;
  char buf[20];

  vp = get_request (PT_OF_V);
  gtk_container_set_border_width (GTK_CONTAINER (hbox), 8);
  gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), hbox, FALSE, FALSE,
		      0);

  stock =
    gtk_image_new_from_stock (GTK_STOCK_DIALOG_QUESTION,
			      GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start (GTK_BOX (hbox), stock, FALSE, FALSE, 0);

  table = gtk_table_new (2, 4, FALSE);
  gtk_table_set_row_spacings (GTK_TABLE (table), 4);
  gtk_table_set_col_spacings (GTK_TABLE (table), 4);
  gtk_box_pack_start (GTK_BOX (hbox), table, TRUE, TRUE, 0);
  label = gtk_label_new_with_mnemonic ("point x");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 0, 1);
  entry_x = gtk_entry_new ();
  sprintf (buf, "%g", vp->vp.x);
  gtk_entry_set_text (GTK_ENTRY (entry_x), buf);
  gtk_table_attach_defaults (GTK_TABLE (table), entry_x, 1, 2, 0, 1);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_x);

  label = gtk_label_new_with_mnemonic ("point y");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 1, 2);

  entry_y = gtk_entry_new ();
  sprintf (buf, "%g", vp->vp.y);
  gtk_entry_set_text (GTK_ENTRY (entry_y), buf);
  gtk_table_attach_defaults (GTK_TABLE (table), entry_y, 1, 2, 1, 2);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_y);

  label = gtk_label_new_with_mnemonic ("point z");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 2, 3);
  entry_z = gtk_entry_new ();
  sprintf (buf, "%g", vp->vp.z);
  gtk_entry_set_text (GTK_ENTRY (entry_z), buf);
  gtk_table_attach_defaults (GTK_TABLE (table), entry_z, 1, 2, 2, 3);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_z);

  label = gtk_label_new_with_mnemonic ("point H");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 3, 4);
  entry_H = gtk_entry_new ();
  sprintf (buf, "%g", vp->Z);
  gtk_entry_set_text (GTK_ENTRY (entry_H), buf);
  gtk_table_attach_defaults (GTK_TABLE (table), entry_H, 1, 2, 3, 4);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_H);

  label = gtk_label_new_with_mnemonic ("Corners");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 4, 5);
  entry_c = gtk_entry_new ();
  sprintf (buf, "%g", vp->Z);
  gtk_entry_set_text (GTK_ENTRY (entry_c), buf);
  gtk_table_attach_defaults (GTK_TABLE (table), entry_c, 1, 2, 4, 5);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_c);

  gtk_widget_show_all (hbox);
  response = gtk_dialog_run (GTK_DIALOG (dialog));

  if (response == GTK_RESPONSE_OK)
    {

      double x, y, z, H;

      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_x)), "%lf", &x);
      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_y)), "%lf", &y);
      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_z)), "%lf", &z);
      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_H)), "%lf", &H);

      vp->vp.x = x;
      vp->vp.y = y;
      vp->vp.z = z;
      vp->Z = H;
    }

  gtk_widget_destroy (dialog);
}

void
dialog_turn (GtkButton * button, gpointer user_data)
{
  GtkWidget *dialog;
  GtkWidget *hbox;
  GtkWidget *stock;
  GtkWidget *table;
  GtkWidget *entry_x, *entry_H;
  GtkWidget *entry_y, *entry_z;
  GtkWidget *label;
  gint response;

  dialog =
    gtk_dialog_new_with_buttons ("View turn", GTK_WINDOW (user_data),
				 GTK_DIALOG_MODAL |
				 GTK_DIALOG_DESTROY_WITH_PARENT, GTK_STOCK_OK,
				 GTK_RESPONSE_OK, GTK_STOCK_CANCEL,
				 GTK_RESPONSE_CANCEL, NULL);

  hbox = gtk_hbox_new (FALSE, 8);
  Viewpoint *vp;
  char buf[20];

  vp = get_request (PT_OF_V);
  gtk_container_set_border_width (GTK_CONTAINER (hbox), 8);
  gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), hbox, FALSE, FALSE,
		      0);

  stock =
    gtk_image_new_from_stock (GTK_STOCK_DIALOG_QUESTION,
			      GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start (GTK_BOX (hbox), stock, FALSE, FALSE, 0);

  table = gtk_table_new (2, 4, FALSE);
  gtk_table_set_row_spacings (GTK_TABLE (table), 4);
  gtk_table_set_col_spacings (GTK_TABLE (table), 4);
  gtk_box_pack_start (GTK_BOX (hbox), table, TRUE, TRUE, 0);
  label = gtk_label_new_with_mnemonic ("rotation x");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 0, 1);
  entry_x = gtk_entry_new ();
  sprintf (buf, "%g", vp->vp.x);
  gtk_entry_set_text (GTK_ENTRY (entry_x), buf);
  gtk_table_attach_defaults (GTK_TABLE (table), entry_x, 1, 2, 0, 1);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_x);

  label = gtk_label_new_with_mnemonic ("rotation y");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 1, 2);

  entry_y = gtk_entry_new ();
  sprintf (buf, "%g", vp->vp.y);
  gtk_entry_set_text (GTK_ENTRY (entry_y), buf);

  gtk_table_attach_defaults (GTK_TABLE (table), entry_y, 1, 2, 1, 2);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_y);
  label = gtk_label_new_with_mnemonic ("point z");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 2, 3);
  entry_z = gtk_entry_new ();
  sprintf (buf, "%g", vp->vp.z);
  gtk_entry_set_text (GTK_ENTRY (entry_z), buf);

  gtk_table_attach_defaults (GTK_TABLE (table), entry_z, 1, 2, 2, 3);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_z);
  label = gtk_label_new_with_mnemonic ("point H");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 3, 4);
  entry_H = gtk_entry_new ();
  sprintf (buf, "%g", vp->Z);
  gtk_entry_set_text (GTK_ENTRY (entry_H), buf);

  gtk_table_attach_defaults (GTK_TABLE (table), entry_H, 1, 2, 3, 4);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), entry_H);

  gtk_widget_show_all (hbox);
  response = gtk_dialog_run (GTK_DIALOG (dialog));

  if (response == GTK_RESPONSE_OK)
    {

      double x, y, z, H;

      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_x)), "%lf", &x);
      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_y)), "%lf", &y);
      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_z)), "%lf", &z);
      sscanf (gtk_entry_get_text (GTK_ENTRY (entry_H)), "%lf", &H);

      vp->vp.x = x;
      vp->vp.y = y;
      vp->vp.z = z;
      vp->Z = H;
    }

  gtk_widget_destroy (dialog);
}

void
ps_to_stdout (GtkButton * button, gpointer user_data)
{
  GtkWidget *dialog, *hbox, *stock, *table, *local_entry2, *label;

  gint response;

  dialog =
    gtk_dialog_new_with_buttons ("Interactive Dialog", GTK_WINDOW (user_data),
				 GTK_DIALOG_MODAL |
				 GTK_DIALOG_DESTROY_WITH_PARENT, GTK_STOCK_OK,
				 GTK_RESPONSE_OK, GTK_STOCK_CANCEL,
				 GTK_RESPONSE_CANCEL, NULL);

  hbox = gtk_hbox_new (FALSE, 8);
  gtk_container_set_border_width (GTK_CONTAINER (hbox), 8);
  gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), hbox, FALSE, FALSE,
		      0);

  stock =
    gtk_image_new_from_stock (GTK_STOCK_DIALOG_QUESTION,
			      GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start (GTK_BOX (hbox), stock, FALSE, FALSE, 0);

  table = gtk_table_new (1, 1, FALSE);
  gtk_table_set_row_spacings (GTK_TABLE (table), 4);
  gtk_table_set_col_spacings (GTK_TABLE (table), 4);
  gtk_box_pack_start (GTK_BOX (hbox), table, TRUE, TRUE, 0);
  label = gtk_label_new_with_mnemonic ("Multiply");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 0, 1);
  local_entry2 = gtk_entry_new ();
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry2, 1, 2, 0, 1);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry2);
  gtk_widget_show_all (hbox);
  response = gtk_dialog_run (GTK_DIALOG (dialog));
  if (response == GTK_RESPONSE_OK)
    {
      PrimBuf prb;
      pr_scale psc;
      double koef;

      sscanf (gtk_entry_get_text (GTK_ENTRY (local_entry2)), "%lf", &koef);
      psc.K = koef;
      prb = get_request (PRIMITIVES);
      if (prb != NULL)
	{
	  GtkWidget *filew;

	  filew =
	    gtk_file_chooser_dialog_new ("Open xml file", NULL,
					 GTK_FILE_CHOOSER_ACTION_SAVE,
					 GTK_STOCK_CANCEL,
					 GTK_RESPONSE_CANCEL, GTK_STOCK_SAVE,
					 GTK_RESPONSE_ACCEPT, NULL);
	  if (gtk_dialog_run (GTK_DIALOG (filew)) == GTK_RESPONSE_ACCEPT)
	    {
	      char *filename;
	      FILE *fp;
	      filename =
		gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (filew));
	      fp = fopen (filename, "w");
	      if (fp)
		{
		  setlocale (LC_NUMERIC, "C");
		  prp_step_by_step_ps (fp, psc, prb);
		  fclose (fp);
		  setlocale (LC_NUMERIC, "");
		}
	      else
		fprintf (stderr, "Can\'t open file: %s\n", filename);
	    }
	  gtk_widget_destroy (filew);
	}
    }

  gtk_widget_destroy (dialog);
}

gboolean
load3d_xml_file (GtkWidget * widget, gpointer * data)
{
  GtkWidget *filew;

  filew =
    gtk_file_chooser_dialog_new ("Open xml file", NULL,
				 GTK_FILE_CHOOSER_ACTION_OPEN,
				 GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
				 GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT, NULL);
  if (gtk_dialog_run (GTK_DIALOG (filew)) == GTK_RESPONSE_ACCEPT)
    {
      char *filename;
      Elements3D *k3d = get_request (KUPA_3D);
      SetsContainer *setcon = get_request (KUPAS_3D);

      memset (k3d, 0, sizeof (Elements3D));
      sets3d_del (setcon->sets);

      filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (filew));
      setcon->sets = interpret_sets3d_xml (filename);
#ifndef NDEBUG
      if (setcon->sets.nsets)
	{
	  *k3d = setcon->sets.sets[0];
	  fprintf (stderr,
		   "%s:%d: params: ncyls=%d, nsphers=%d,nmirages=%d\n",
		   __FILE__, __LINE__, k3d->ncyls, k3d->nsphers,
		   k3d->nmirages);
	  setcon->position = 0;
	}
#endif
    }
  gtk_widget_destroy (filew);
  return FALSE;

}

void
ps_to_file2 (GtkButton * button, gpointer user_data)
{
  GtkWidget *dialog;
  GtkWidget *hbox;
  GtkWidget *stock;
  GtkWidget *table;
  GtkWidget *local_entry[4];
  GtkWidget *label;
  gint response;

  dialog =
    gtk_dialog_new_with_buttons ("Interactive Dialog", GTK_WINDOW (user_data),
				 GTK_DIALOG_MODAL |
				 GTK_DIALOG_DESTROY_WITH_PARENT, GTK_STOCK_OK,
				 GTK_RESPONSE_OK, GTK_STOCK_CANCEL,
				 GTK_RESPONSE_CANCEL, NULL);

  hbox = gtk_hbox_new (FALSE, 8);
  gtk_container_set_border_width (GTK_CONTAINER (hbox), 8);
  gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), hbox, FALSE, FALSE,
		      0);

  stock =
    gtk_image_new_from_stock (GTK_STOCK_DIALOG_QUESTION,
			      GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start (GTK_BOX (hbox), stock, FALSE, FALSE, 0);

  table = gtk_table_new (2, 4, FALSE);
  gtk_table_set_row_spacings (GTK_TABLE (table), 4);
  gtk_table_set_col_spacings (GTK_TABLE (table), 4);
  gtk_box_pack_start (GTK_BOX (hbox), table, TRUE, TRUE, 0);
  label = gtk_label_new_with_mnemonic ("first point xyz");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 0, 1);
  local_entry[0] = gtk_entry_new ();
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[0], 1, 2, 0, 1);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[0]);

  label = gtk_label_new_with_mnemonic ("second point xyz");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 1, 2);
  local_entry[1] = gtk_entry_new ();
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[1], 1, 2, 1, 2);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[1]);

  label = gtk_label_new_with_mnemonic ("Multiply, shift");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 2, 3);

  local_entry[2] = gtk_entry_new ();
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[2], 1, 2, 2, 3);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[2]);

  label = gtk_label_new_with_mnemonic ("file name");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 3, 4);
  local_entry[3] = gtk_entry_new ();
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[3], 1, 2, 3, 4);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[3]);

  gtk_widget_show_all (hbox);
  response = gtk_dialog_run (GTK_DIALOG (dialog));

  if (response == GTK_RESPONSE_OK)
    {

      {
	PrimBuf prb;
	Elements3D *k3d;
	SetP pk;
	Viewpoint VP1, VP2;
	double koef, xy;
	pr_scale psc;
	FILE *fp;

	sscanf (gtk_entry_get_text (GTK_ENTRY (local_entry[0])),
		"%lf %lf %lf", &VP1.vp.x, &VP1.vp.y, &VP1.vp.z);
	VP1.Z = VP1.vp.z;
	sscanf (gtk_entry_get_text (GTK_ENTRY (local_entry[1])),
		"%lf %lf %lf", &VP2.vp.x, &VP2.vp.y, &VP2.vp.z);
	VP2.Z = VP2.vp.z;
	sscanf (gtk_entry_get_text (GTK_ENTRY (local_entry[2])), "%lf %lf",
		&xy, &koef);
	psc.x = xy;
	psc.y = xy;
	psc.K = koef;
	k3d = get_request (KUPA_3D);
	fp = fopen (gtk_entry_get_text (GTK_ENTRY (local_entry[3])), "w");
	if (fp == NULL)
	  fprintf (stderr, "cant open file: %s\n",
		   gtk_entry_get_text (GTK_ENTRY (local_entry[1])));
	if (k3d && fp)
	  {
	    tensor *tens;

	    tens = get_request (TENS);
	    fprintf (stderr, "running :)\n");
	    pk = project_all (k3d, *tens, 16, &VP1);
	    mask_all (pk);
	    prb = pri_init (2., 2.);
	    prb = plot_projection_all (prb, pk);
	    projection_all_del (pk);
	    if (prb != NULL)
	      prp_step_by_step_ps (fp, psc, prb);
	    free (prb);
	    pk = project_all (k3d, *tens, 8, &VP2);
	    mask_all (pk);
	    prb = pri_init (2., 2.);
	    prb = plot_projection_all (prb, pk);
	    projection_all_del (pk);
	    if (prb != NULL)
	      prp_step_by_step_ps (fp, psc, prb);
	    free (prb);
	  }
	fclose (fp);
      }
    }
  gtk_widget_destroy (dialog);
}

void
svg_to_file (GtkButton * button, gpointer user_data)
{
  GtkWidget *dialog;
  GtkWidget *hbox;
  GtkWidget *stock;
  GtkWidget *table;
  GtkWidget *local_entry[4];
  GtkWidget *label;
  gint response;

  dialog =
    gtk_dialog_new_with_buttons ("Interactive Dialog", GTK_WINDOW (user_data),
				 GTK_DIALOG_MODAL |
				 GTK_DIALOG_DESTROY_WITH_PARENT, GTK_STOCK_OK,
				 GTK_RESPONSE_OK, GTK_STOCK_CANCEL,
				 GTK_RESPONSE_CANCEL, NULL);

  hbox = gtk_hbox_new (FALSE, 8);
  gtk_container_set_border_width (GTK_CONTAINER (hbox), 8);
  gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), hbox, FALSE, FALSE,
		      0);

  stock =
    gtk_image_new_from_stock (GTK_STOCK_DIALOG_QUESTION,
			      GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start (GTK_BOX (hbox), stock, FALSE, FALSE, 0);

  table = gtk_table_new (2, 4, FALSE);
  gtk_table_set_row_spacings (GTK_TABLE (table), 4);
  gtk_table_set_col_spacings (GTK_TABLE (table), 4);
  gtk_box_pack_start (GTK_BOX (hbox), table, TRUE, TRUE, 0);
  label = gtk_label_new_with_mnemonic ("params x y K H");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 0, 1);
  local_entry[0] = gtk_entry_new ();
  gtk_entry_set_text (GTK_ENTRY (local_entry[0]), "2 2 10 40");
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[0], 1, 2, 0, 1);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[0]);

  label = gtk_label_new_with_mnemonic ("viewBox");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 1, 2);
  local_entry[1] = gtk_entry_new ();
  gtk_entry_set_text (GTK_ENTRY (local_entry[1]),
		      "0 0 40 40\" width=\"12in\" height=\"12in");
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[1], 1, 2, 1, 2);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[1]);

  label = gtk_label_new_with_mnemonic ("\"<g>\" params");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 2, 3);

  local_entry[2] = gtk_entry_new ();
  gtk_entry_set_text (GTK_ENTRY (local_entry[2]),
		      "stroke=\"blue\" stroke-width=\".1\" fill=\"none\"");
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[2], 1, 2, 2, 3);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[2]);

  label = gtk_label_new_with_mnemonic ("file name");
  gtk_table_attach_defaults (GTK_TABLE (table), label, 0, 1, 3, 4);
  local_entry[3] = gtk_entry_new ();
  gtk_table_attach_defaults (GTK_TABLE (table), local_entry[3], 1, 2, 3, 4);
  gtk_label_set_mnemonic_widget (GTK_LABEL (label), local_entry[3]);

  gtk_widget_show_all (hbox);
  response = gtk_dialog_run (GTK_DIALOG (dialog));

  if (response == GTK_RESPONSE_OK)
    {
      {
	PrimBuf prb;
	pr_scale psc;
	FILE *fp;

	fp = fopen (gtk_entry_get_text (GTK_ENTRY (local_entry[3])), "w");
	prb = get_request (PRIMITIVES);

	if (fp == NULL)
	  fprintf (stderr, "cant open file: %s\n",
		   gtk_entry_get_text (GTK_ENTRY (local_entry[3])));
	if (prb && fp)
	  {
	    sscanf (gtk_entry_get_text (GTK_ENTRY (local_entry[0])),
		    "%lf %lf %lf %d", &psc.x, &psc.y, &psc.K, &psc.H);
	    fputs ("<?xml version=\"1.0\" standalone=\"no\"?>\n"
		   "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n"
		   "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n",
		   fp);
	    fprintf (fp, "<svg viewBox=\"%s\">\n",
		     gtk_entry_get_text (GTK_ENTRY (local_entry[1])));
	    fprintf (fp, "<g %s>\n",
		     gtk_entry_get_text (GTK_ENTRY (local_entry[2])));
	    if (prb != NULL)
	      prp_step_by_step_svg (fp, psc, prb);
	    fputs ("</g>\n</svg>\n", fp);
	  }
	fclose (fp);
      }
    }
  gtk_widget_destroy (dialog);
}

void
move_box_position (int delta)
{
  SetsContainer *setcon = get_request (KUPAS_3D);
  Elements3D *k3d = get_request (KUPA_3D);
  GtkStatusbar *StatusBar = get_request (STATUSBAR);
  int max, pos;
  gchar *msg;

  max = setcon->sets.nsets;
  pos = setcon->position;
  pos += delta;
  if (pos < 0)
    pos = max - 1;
  if (pos >= max)
    pos = 0;
  setcon->position = pos;
  *k3d = setcon->sets.sets[pos];
  msg = g_strdup_printf ("%d/%d", pos + 1, max);
  gtk_statusbar_push (StatusBar, 0, msg);
  g_free (msg);
  // That is Not all
}

int
event_key_pressed (GtkWidget * widg, GdkEventKey * key, void *data)
{
  tensor *tens;

  tens = get_request (TENS);
  switch (key->hardware_keycode)
    {
    case 100 /* left */ :
      *tens = RotateCoord (0., 5. / 180. * M_PI, *tens);
      break;
    case 104 /* down */ :
      *tens = RotateCoord (5. / 180. * M_PI, 0., *tens);
      break;
    case 102 /* right */ :
      *tens = RotateCoord (0., -5. / 180. * M_PI, *tens);
      break;
    case 98 /* up */ :
      *tens = RotateCoord (-5. / 180. * M_PI, 0., *tens);
      break;
    case 43 /*H*/:
      *tens = RotateCoord (0., 1. / 180. * M_PI, *tens);
      break;
    case 44 /*J*/:
      *tens = RotateCoord (1. / 180. * M_PI, 0., *tens);
      break;
    case 45 /*K*/:
      *tens = RotateCoord (-1. / 180. * M_PI, 0., *tens);
      break;
    case 46 /*L*/:
      *tens = RotateCoord (0., -1. / 180. * M_PI, *tens);
      break;
    case 24 /*Q*/:
    case 9 /* esc */ :
      gtk_main_quit ();
      return 0;
    case 34 /* [ */ :
      move_box_position (-1);
      break;
    case 35 /* ] */ :
      move_box_position (1);
      break;
    default:
      printf ("HWKEY = %d\n", key->hardware_keycode);
      return 0;
    }
  visualization_generator (tens);
  gdk_window_invalidate_rect (widg->window, NULL, TRUE);
  return 0;
}
