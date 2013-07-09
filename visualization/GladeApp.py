import os
import sys

import gtk
import gtk.glade

class GladeApp(dict):
    def __init__(self, glade_filename, main_widget_name=None, domain=None):
        gtk.glade.set_custom_handler(self.custom_handler)
        if os.path.isfile(glade_filename):
            self.glade_path = glade_filename
        else:
            glade_dir = os.path.split(sys.argv[0])[0]
            self.glade_path = os.path.join(glade_dir, glade_filename)
        self.glade = gtk.glade.XML(self.glade_path, main_widget_name, domain)
        if main_widget_name:
            self.main_widget = self.glade.get_widget(main_widget_name)
        else:
            self.main_widget = None
        self.signal_autoconnect()
        self.new()


    def signal_autoconnect(self):
        signals = {}
        for attr_name in dir(self):
            attr = getattr(self, attr_name)
            if callable(attr):
                signals[attr_name] = attr
        self.glade.signal_autoconnect(signals)


    def custom_handler(self, glade, function_name, widget_name,
                       str1, str2, int1, int2):
        if hasattr(self, function_name):
            handler = getattr(self, function_name)
            return handler(str1, str2, int1, int2)


    def __getattr__(self, data_name):
        if data_name in self:
            data = self[data_name]
            return data
        else:
            widget = self.glade.get_widget(data_name)
            if widget:
                self[data_name] = widget
                return widget
            else:
                raise AttributeError, data_name
        

    def __setattr__(self, name, value):
        self[name] = value


    def new(self):
        pass


    def on_keyboard_interrupt(self):
        pass


    def main(self):
        gtk.main()


    def quit(self, *args):
        gtk.main_quit()


    def run(self):
        try:
            self.main()
        except KeyboardInterrupt:
            self.on_keyboard_interrupt()



