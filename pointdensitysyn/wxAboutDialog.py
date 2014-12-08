import wx
import gui
import version

class wxAboutDialog(gui.AboutDialog):
    def __init__(self, parent):
        gui.AboutDialog.__init__(self, parent)
        self.TitleLabel.SetLabel(version.title)
        self.IconBitmap.SetBitmap(wx.Bitmap(version.icon, wx.BITMAP_TYPE_ANY))
        self.VersionLabel.SetLabel("Version %s" % version.version)
        self.LastModLabel.SetLabel("Last modified %s %s, %s." % version.date)
        self.CopyrightLabel.SetLabel("Copyright 2001-%s %s." % (version.date[2],
                                                                version.author))
        self.LicenseLabel.SetLabel("Released under the terms of the MIT"
                                   " license.")
        self.EmailHyperlink.SetLabel("%s" % version.email)
        self.EmailHyperlink.SetURL("mailto://%s" % version.email)
        self.WebHyperlink.SetLabel("http://%s" % version.homepage)
        self.WebHyperlink.SetURL("http://%s" % version.homepage)
        self.SetIcon(wx.Icon(version.icon, wx.BITMAP_TYPE_ICO))              
        self.Fit()

    def OnClose( self, event ):
        self.Destroy()
    
    
