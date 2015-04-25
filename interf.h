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
