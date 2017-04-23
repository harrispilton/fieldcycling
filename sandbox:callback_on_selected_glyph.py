from bokeh.plotting import figure, output_file, show, ColumnDataSource, hplot
from bokeh.models import HoverTool, Callback, ColumnDataSource
import pandas as pd
output_file("bla.html")

# Mock data
m = np.ones((6,11))
for i in range(2,6):
    for j in range(11):
        m[i,j] = i+j
x = [0,1,2,3,4,5]; y = [0,2,4,6,8,10]
m0 = m.transpose()
m1 = pd.DataFrame(m0, index=['0','1','2','3','4','5','6','7','8','9','10'], columns=[np.arange(0,len(m),1).astype(str)])

#First plot
s1 = ColumnDataSource(data=dict(x=x,y=y))
p1 = figure(tools=["lasso_select"], plot_width=600, plot_height=400)
p1.scatter('x', 'y', fill_color='black', line_color=None, size=10, source=s1)

#Second plot
s2 = ColumnDataSource(data=dict(x=[],y=[],y2=[]))
p2 = figure(plot_width=400, plot_height=400, tools =[])

m1 = ColumnDataSource(m1) #Actual Datasource for the second plot
p2.line(np.arange(0,100,1), 'y' , source=s2) # From original data - series 1
p2.line(np.arange(0,100,1), 'y2' , source=s2) # From original data - series 2

s1.callback = Callback(args=dict(s2=s2, m1=m1), code="""
  var inds = cb_obj.get('selected')['1d'].indices;
  var d1 = m1.get('data'); 
  var d2 = s2.get('data');
  d2['y'] = []
  d2['y2'] = []
  for (i = 0; i < 11; i++) {
    d2['y'].push(d1[inds['0']][i]),
    d2['y2'].push(d1[inds['1']][i])
  }
  s2.trigger('change'); 
""")

layout = hplot(p1, p2)
show(layout)
