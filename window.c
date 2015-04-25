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
#include <stdlib.h>
#include <math.h>
#include <vmath.h>
#include <pgl.h>
#include "probez.h"
#include "interf.h"

typedef struct
{
  int x, y;
} POINT;
static POINT PrevMousePoint;
static pr_scale ext_psc;
static pr_real ext_W = 1.4;
void interactive_dialog_clicked ();
static gboolean
expose_event (GtkWidget * da, GdkEventExpose * event, gpointer * data)
{
  GdkGC *gc;
  GdkColor color;
  PrimBuf prb = NULL;

  gc = gdk_gc_new (da->window);
  color.pixel = 0xFFFFFF;
  gdk_gc_set_foreground (gc, &color);
  gdk_draw_rectangle (da->window, gc, TRUE, 0, 0,
		      da->allocation.width, da->allocation.height);
  prb = get_request (PRIMITIVES);
  if (prb != NULL)
    prp_step_by_step_gdk (da->window, gc, prb);
  g_object_unref (gc);
  return TRUE;
}

static gboolean
configure_event (GtkWidget * widget, GdkEventConfigure * event, gpointer data)
{
  int h, w;

  h = widget->allocation.height;
  w = widget->allocation.width;
  ext_psc.x = -ext_W;
  ext_psc.y = -ext_W;
  ext_psc.H = h;
  ext_psc.K = w < h ? (double) w / 2. / ext_W : (double) h / 2. / ext_W;
  return TRUE;
}

static gboolean
scribble_button_press_event (GtkWidget * widget,
			     GdkEventButton * event, gpointer data)
{
  if (event->button == 1)
    {
      PrevMousePoint.x = (int) event->x;
      PrevMousePoint.y = (int) event->y;
    }
  if (event->button == 3)
    gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button,
		    event->time);

  /* We've handled the event, stop processing */
  return TRUE;
}

void
visualization_generator (const tensor * tens)
{
  PrimBuf prb = NULL;
  Elements3D *k3d = get_request (KUPA_3D);
  Viewpoint *VP = get_request (PT_OF_V);
  static Model mod;
  mod.set = k3d;
  mod.vp = VP;
  mod.rot = (tensor*)tens;
  mod.ncorners = 6;
  prb = image_generator (&mod);
  if (get_request (PRIMITIVES))
    free (get_request (PRIMITIVES));
  set_request (PRIMITIVES, prb);
}

static gboolean
scribble_motion_notify_event (GtkWidget * widget,
			      GdkEventMotion * event, gpointer data)
{
  int x, y;
  tensor *tens;
  GdkModifierType state;

  gdk_window_get_pointer (event->window, &x, &y, &state);
  if (state & GDK_BUTTON1_MASK)
    {
      tens = get_request (TENS);
      *tens =
	RotateCoord (((double) y -
		      (double) PrevMousePoint.y) / (double) ext_psc.H,
		     ((double) PrevMousePoint.x -
		      (double) x) / (double) ext_psc.H, *tens);
      PrevMousePoint.x = x;
      PrevMousePoint.y = y;
      {
	visualization_generator (tens);
	gdk_window_invalidate_rect (widget->window, NULL, TRUE);
      }
    }
  /* We've handled it, stop processing */
  return TRUE;
}

static gint
Delete (GtkWidget * widget, gpointer * data)
{
  gtk_main_quit ();
  return FALSE;
}

static void
fill_menu (GtkWidget * window, GtkWidget * MenuBar)
{
  GtkWidget *FileMenu, *open_item, *file_item,
    *view_item, *point_plane_item, *save2post_item, *ViewMenu, *exit_item,
    *file_ps2_item, *file_save_svg_item;
  FileMenu = gtk_menu_new ();
  ViewMenu = gtk_menu_new ();
  file_item = gtk_menu_item_new_with_label ("File");
  view_item = gtk_menu_item_new_with_label ("View");
  // -------- FILE
  open_item = gtk_menu_item_new_with_label ("Open");

  g_signal_connect (GTK_OBJECT (open_item), "activate",
		    G_CALLBACK (load3d_xml_file), window);
  gtk_menu_shell_append (GTK_MENU_SHELL (FileMenu), open_item);

  save2post_item = gtk_menu_item_new_with_label ("Save to PS");
  g_signal_connect (GTK_OBJECT (save2post_item), "activate",
		    G_CALLBACK (ps_to_stdout), window);
  gtk_menu_shell_append (GTK_MENU_SHELL (FileMenu), save2post_item);
  file_save_svg_item = gtk_menu_item_new_with_label ("Save to SVG");
  g_signal_connect (GTK_OBJECT (file_save_svg_item), "activate",
		    G_CALLBACK (svg_to_file), window);
  gtk_menu_shell_append (GTK_MENU_SHELL (FileMenu), file_save_svg_item);
  file_ps2_item = gtk_menu_item_new_with_label ("Duplet");
  g_signal_connect (GTK_OBJECT (file_ps2_item), "activate",
		    G_CALLBACK (ps_to_file2), window);
  gtk_menu_shell_append (GTK_MENU_SHELL (FileMenu), file_ps2_item);
  exit_item = gtk_menu_item_new_with_label ("Exit");
  g_signal_connect (GTK_OBJECT (exit_item), "activate",
		    G_CALLBACK (gtk_main_quit), NULL);
  gtk_menu_shell_append (GTK_MENU_SHELL (FileMenu), exit_item);
  // --------- VIEW
  point_plane_item = gtk_menu_item_new_with_label ("View point");
  g_signal_connect (GTK_OBJECT (point_plane_item), "activate",
		    G_CALLBACK (dialog_wiewpoint), window);
  gtk_menu_shell_append (GTK_MENU_SHELL (ViewMenu), point_plane_item);
  gtk_menu_item_set_submenu (GTK_MENU_ITEM (file_item), FileMenu);
  gtk_menu_item_set_submenu (GTK_MENU_ITEM (view_item), ViewMenu);
  gtk_menu_shell_append (GTK_MENU_SHELL (MenuBar), file_item);
  gtk_menu_shell_append (GTK_MENU_SHELL (MenuBar), view_item);
}

void
CreateDraw ()
{
  GtkWidget *vBox, *window, *DrawArea, *MenuBar, *StatusBar;

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  MenuBar = gtk_menu_bar_new ();
  fill_menu (window, MenuBar);
  gtk_window_set_title (GTK_WINDOW (window), "ProBez");
  vBox = gtk_vbox_new (FALSE, 0);
  DrawArea = gtk_drawing_area_new ();
  gtk_widget_set_size_request (DrawArea, 400, 300);
  StatusBar = gtk_statusbar_new ();
  gtk_container_add (GTK_CONTAINER (window), vBox);
  gtk_box_pack_start (GTK_BOX (vBox), MenuBar, FALSE, TRUE, 0);
  gtk_box_pack_start (GTK_BOX (vBox), DrawArea, TRUE, TRUE, 0);
  gtk_box_pack_start (GTK_BOX (vBox), StatusBar, FALSE, TRUE, 0);
  set_request (STATUSBAR, StatusBar);
  gtk_signal_connect (GTK_OBJECT (window), "delete_event",
		      GTK_SIGNAL_FUNC (Delete), NULL);
  gtk_signal_connect (GTK_OBJECT (DrawArea), "configure_event",
		      GTK_SIGNAL_FUNC (configure_event), NULL);
  gtk_signal_connect (GTK_OBJECT (DrawArea), "expose_event",
		      GTK_SIGNAL_FUNC (expose_event), NULL);
  g_signal_connect (GTK_OBJECT (DrawArea), "button_press_event",
		    G_CALLBACK (scribble_button_press_event), NULL);
  g_signal_connect (GTK_OBJECT (DrawArea), "motion_notify_event",
		    G_CALLBACK (scribble_motion_notify_event), NULL);

  gtk_widget_show (DrawArea);
  g_signal_connect (GTK_OBJECT (window), "key-press-event",
		    G_CALLBACK (event_key_pressed), "window");

  /* Ask to receive events the drawing area doesn't normally subscribe
   * to */
  gtk_widget_set_events (DrawArea,
			 gtk_widget_get_events (DrawArea) |
			 GDK_LEAVE_NOTIFY_MASK | GDK_BUTTON_PRESS_MASK |
			 GDK_POINTER_MOTION_MASK |
			 GDK_POINTER_MOTION_HINT_MASK | GDK_KEY_PRESS_MASK);
  gtk_widget_show (vBox);
  gtk_widget_show_all (window);	/* do this last */
  set_request (KUPA_3D, calloc (1, sizeof (Elements3D)));
  set_request (KUPAS_3D, calloc (1, sizeof (SetsContainer)));
  {
    Viewpoint *vp;
    tensor *tens;

    vp = malloc (sizeof (Viewpoint));
    vp->vp.x = vp->vp.y = 0.;
    vp->vp.z = 10.;
    vp->Z = 8.;
    set_request (PT_OF_V, vp);
    tens = malloc (sizeof (tensor));
    *tens = InitRotation (0., 0.);
    set_request (TENS, tens);
  }
}

int
main (int argc, char **argv)
{
  gtk_init (&argc, &argv);
  init_requests ();
  CreateDraw ();
  if (argc > 1)
    {
      FILE *fp;
      SetsContainer *sets;

      fp = fopen (argv[1], "r");
      if (fp)
	{
	  fclose (fp);
	  sets = get_request (KUPAS_3D);
	  sets->sets = interpret_sets3d_xml (argv[1]);
	}
    }
  gtk_main ();
  return 0;
}
