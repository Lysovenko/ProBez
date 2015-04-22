#include <gtk/gtk.h>

/// inter
gboolean load3dfile(GtkWidget *widget,gpointer *data);
gboolean load3d_xml_file(GtkWidget *widget,gpointer *data);	    
void visualization_generator (tensor tens);
 void
ps_to_stdout (GtkButton *button,
			    gpointer   user_data);
 void
dialog_wiewpoint (GtkButton *button,
    			    gpointer   user_data);

void ps_to_file2 (GtkButton *button,
			    gpointer   user_data);
			    
void
svg_to_file (GtkButton *button,
			    gpointer   user_data);
			    
int event_key_pressed(GtkWidget *widg,GdkEventKey*key,void*data);
