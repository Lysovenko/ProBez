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

#ifndef __INTERF_H__
#define __INTERF_H__

#ifdef  __GTK_H__
gboolean load3dfile (GtkWidget * widget, gpointer * data);
gboolean load3d_xml_file (GtkWidget * widget, gpointer * data);
void visualization_generator (const tensor * tens);
void ps_to_stdout (GtkButton * button, gpointer user_data);
void dialog_wiewpoint (GtkButton * button, gpointer user_data);
void ps_to_file2 (GtkButton * button, gpointer user_data);
void svg_to_file (GtkButton * button, gpointer user_data);
int event_key_pressed (GtkWidget * widg, GdkEventKey * key, void *data);
#endif

#ifdef __PROBEZ_H__
typedef struct
{
  Sets3D ks3d;
  int position;
} SetsContainer;
#endif

void init_requests ();
void *get_request (int n);
void set_request (int n, void *what);
void unset_request (int n);

enum Requests
{
  DAT_COMPCFG = 0,
  DAT_CIRCLIN,
  STATUSBAR,
  CameraData,
  PRIMITIVES,
  KUPA_3D,
  PT_OF_V,
  TENS,
  KUPAS_3D,
  G_,
  N_requests
};
#endif
