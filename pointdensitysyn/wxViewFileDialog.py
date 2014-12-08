import os.path
import wx
import gui

class wxViewFileDialog(gui.ViewFileDialog):
    def __init__(self, parent, fn):
        gui.ViewFileDialog.__init__(self, parent)
        try:
            self.SetTitle(os.path.basename(fn))
            f = open(fn, "r", 0)
            try:
                for s in f.readlines():
                    self.ViewFileTextCtrl.AppendText(s)
            finally:
                f.close()
        except IOError:
            parent.ShowError("Could not open file.")
            self.Close()
            
    def OnClose(self, event):
        self.Destroy()
    
    
