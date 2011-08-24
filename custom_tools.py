import wx
from matplotlib.patches import Rectangle
from matplotlib.widgets import Lasso
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg

class MyNavToolbar(NavigationToolbar2WxAgg):
    """wx/mpl NavToolbar hack with an additional tools user interaction.
    This class is necessary because simply adding a new togglable tool to the
    toolbar won't (1) radio-toggle between the new tool and the pan/zoom tools.
    (2) disable the pan/zoom tool modes in the associated subplot(s).
    """
    ID_LASSO_TOOL = wx.NewId()
    def __init__(self, canvas):
        super(NavigationToolbar2WxAgg, self).__init__(canvas)

        self.pan_tool  = self.FindById(self._NTB2_PAN)
        self.zoom_tool = self.FindById(self._NTB2_ZOOM)

        self.lasso_tool = self.InsertSimpleTool(5, self.ID_LASSO_TOOL, 
                            wx.ArtProvider.GetBitmap(wx.ART_ADD_BOOKMARK),
                            isToggle=True)
        self.Bind(wx.EVT_TOOL, self.on_toggle_lasso_tool, self.lasso_tool)
        self.Bind(wx.EVT_TOOL, self.on_toggle_pan_zoom, self.zoom_tool)
        self.Bind(wx.EVT_TOOL, self.on_toggle_pan_zoom, self.pan_tool)

    def get_mode(self):
        """Use this rather than navtoolbar.mode
        """
        if self.lasso_tool.IsToggled():
            return 'lasso'
        else:
            return self.mode

    def untoggle_mpl_tools(self):
        """Hack city: Since I can't figure out how to change the way the 
        associated subplot(s) handles mouse events: I generate events to turn
        off whichever tool mode is enabled (if any). 
        This function needs to be called whenever any user-defined tool 
        (eg: lasso) is clicked.
        """
        if self.pan_tool.IsToggled():
            wx.PostEvent(
                self.GetEventHandler(), 
                wx.CommandEvent(wx.EVT_TOOL.typeId, self._NTB2_PAN)
            )
            self.ToggleTool(self._NTB2_PAN, False)
        elif self.zoom_tool.IsToggled():
            wx.PostEvent(
                self.GetEventHandler(),
                wx.CommandEvent(wx.EVT_TOOL.typeId, self._NTB2_ZOOM)
            )
            self.ToggleTool(self._NTB2_ZOOM, False)

    def on_toggle_lasso_tool(self, evt):
        """Lasso tool handler.
        """
        if evt.Checked():
            self.untoggle_mpl_tools()

    def on_toggle_pan_zoom(self, evt):
        """Called when pan or zoom is toggled. 
        We need to manually untoggle user-defined tools.
        """
        if evt.Checked():
            self.ToggleTool(self.ID_LASSO_TOOL, False)
        # Make sure the regular pan/zoom handlers get the event
        evt.Skip()

class ScatterPanel(FigureCanvasWxAgg):
    """Contains the guts for drawing scatter plots.
    """
    def __init__(self, parent, **kwargs):
        self.figure = Figure()
        FigureCanvasWxAgg.__init__(self, parent, -1, self.figure, **kwargs)
        self.canvas = self.figure.canvas
        self.SetMinSize((100,100))
        self.figure.set_facecolor((1,1,1))
        self.figure.set_edgecolor((1,1,1))
        self.canvas.SetBackgroundColour('white')

        self.subplot = self.figure.add_subplot(111)
        self.navtoolbar = None
        self.lasso = None
        self.redraw()

        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)

    def lasso_callback(self, verts):
        pass

    def on_press(self, evt):
        """canvas mousedown handler
        """
        if evt.button == 1:
            if self.canvas.widgetlock.locked(): 
                return
            if evt.inaxes is None: 
                return
            if self.navtoolbar and self.navtoolbar.get_mode() == 'lasso':
                self.lasso = Lasso(evt.inaxes, (evt.xdata, evt.ydata), self.lasso_callback)
                self.canvas.widgetlock(self.lasso)

    def on_release(self, evt):
        """canvas mouseup handler
        """
        # Note: lasso_callback is not called on click without drag so we release
        #   the lock here to handle this case as well.
        if evt.button == 1:
            if self.lasso:
                self.canvas.draw_idle()
                self.canvas.widgetlock.release(self.lasso)
                self.lasso = None
        else:
            self.show_popup_menu((evt.x, self.canvas.GetSize()[1]-evt.y), None)

    def redraw(self):
        self.subplot.clear()
        self.subplot.scatter([1,2,3],[3,1,2])

    def get_toolbar(self):
        if not self.navtoolbar:
            self.navtoolbar = MyNavToolbar(self.canvas)
            self.navtoolbar.Realize()
        return self.navtoolbar

if __name__ == "__main__":
    app = wx.PySimpleApp()
    f = wx.Frame(None, size=(600,600))
    p = ScatterPanel(f)
    f.SetToolBar(p.get_toolbar())            
    f.Show()
    app.MainLoop()