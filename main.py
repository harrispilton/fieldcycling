#!/usr/bin/python
# -*- coding: UTF-8 -*-
import os
import sys
from datetime import datetime
import pandas as pd
import numpy as np
import re
import matplotlib as mpl #for the colormaps
from bokeh.charts import Scatter, output_file, show
from bokeh.models.widgets import Panel, Tabs
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.palettes import Spectral5, viridis
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Dropdown, Select, Button
from bokeh.models.callbacks import CustomJS
from tornado.ioloop import IOLoop
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, RangeSlider
from bokeh.server.server import Server
import stelardatafile as sdf
from utils import get_x_axis, model_exp_dec, fun_exp_dec, get_mag_amplitude, magnetization_fit
from scipy.optimize import leastsq

#specify and import data file
path=os.path.join(os.path.curdir,'data')
polymer=sdf.StelarDataFile('297K.sdf',path)
polymer.sdfimport()
nr_experiments = polymer.get_number_of_experiments()

# parameters to dataframe
par_df, columns, discrete, continuous, time, quantileable = polymer.scan_parameters(20)

io_loop = IOLoop.current()
#initially set experiment number ie=1
ie = 1

#TODO: work with more than one tab in the browser: file input, multiple data analysis and manipulation tabs, file output and recent results tab
#TODO: clean modify_doc, it is awfully crowded here...
def modify_doc(doc):
    fid=polymer.getfid(ie) #fid is a dataframe containing fids of experiment ie

    # TODO: more testing on get_x_axis
    tau = get_x_axis(polymer.getparameter(ie))

    #calculate magnetization:
    startpoint=int(0.05*polymer.getparvalue(ie,'BS'))
    endpoint=int(0.1*polymer.getparvalue(ie,'BS')) #TODO: make a range slider to get start- and endpoint interactively
    phi = get_mag_amplitude(fid, startpoint, endpoint,
                            polymer.getparvalue(ie,'NBLK'),
                            polymer.getparvalue(ie,'BS'))

    #prepare magnetization decay curve for fit
    df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
    df['phi_normalized']=(df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] )
    polymer.addparameter(ie,'df_magnetization',df)
    

    
    fit_option = 2
    p0 = [1 , 2 * polymer.getparvalue(ie,'T1MX')**-1, 0]
    df, popt = magnetization_fit(df, p0, fit_option)
    polymer.addparameter(ie,'popt(mono_exp)',popt)
    
    df['fit_phi'] = model_exp_dec(df.tau, *popt)
    

    # convert data to handle in bokeh
    source_fid = ColumnDataSource(data=ColumnDataSource.from_df(fid))
    source_df = ColumnDataSource(data=ColumnDataSource.from_df(df))
    
    # create and plot figures
    p1 = figure(plot_width=300, plot_height=300,
                title='Free Induction Decay', webgl=True)
    p1.line('index', 'im', source=source_fid, color='blue')
    p1.line('index', 'real', source=source_fid, color='green')
    p1.line('index', 'magnitude', source=source_fid, color='red')

    fid_slider = RangeSlider(start=1,end=polymer.getparvalue(ie,'BS'),step=1,callback_policy='mouseup')

    p2 = figure(plot_width=300, plot_height=300,
                title='Magnetization Decay')
    p2.circle_cross('tau', 'phi_normalized', source=source_df, color="navy")
    p2.line('tau', 'fit_phi', source=source_df, color="teal")

    
    # in the plot 4 use followingimpo
    SIZES = list(range(6, 22, 3)) # for some sizes
    COLORS = Spectral5 # for some colors (more colors would be nice somehow)

    def plot_par():
        xs = par_df[x.value].values
        ys = par_df[y.value].values
        x_title = x.value.title()
        y_title = y.value.title()

        kw = dict() #holds optional keyword arguments for figure()
        if x.value in discrete:
            kw['x_range'] = sorted(set(xs))
        if y.value in discrete:
            kw['y_range'] = sorted(set(ys))
        if y.value in time:
            kw['y_axis_type'] = 'datetime'
        if x.value in time:
            kw['x_axis_type'] = 'datetime'
        
        kw['title']="%s vs %s" % (x_title, y_title)


        p4 = figure(plot_height=300, plot_width=600, tools='pan,box_zoom,reset',
                    **kw)

        p4.xaxis.axis_label = x_title
        p4.yaxis.axis_label = y_title


        if x.value in discrete:
            p4.xaxis.major_label_orientation = pd.np.pi / 4 # rotates labels...
                    
        sz = 9
        if size.value != 'None':
            groups = pd.qcut(pd.to_numeric(par_df[size.value].values), len(SIZES))
            sz = [SIZES[xx] for xx in groups.codes]

        c = "#31AADE"
        if color.value != 'None':
            groups = pd.qcut(pd.to_numeric(par_df[color.value]).values, len(COLORS))
            c = [COLORS[xx] for xx in groups.codes]
       
        p4.circle(x=xs, y=ys, color=c, size=sz, line_color="white", alpha=0.6, hover_color='white', hover_alpha=0.5)
        return p4

    def update(attr, old, new):
        tabs.tabs[1].child.children[1] = plot_par()
        #layout_p4.children[1] = plot_par()

    def cb(attr, old, new):
        ## load experiment ie in plot p1 and p2
        ie = new['value'][0]
        fid = polymer.getfid(ie)
        #print(fid)
        #source_fid = ColumnDataSource.from_df(data=fid)
        source_fid.data=ColumnDataSource.from_df(fid)
        #print(source_fid)
        try:
            tau = get_x_axis(polymer.getparameter(ie))
            #print(tau)
            try:
                startpoint=polymer.getparvalue(ie,'fid_range')[0]
                endpoint=polymer.getparvalue(ie,'fid_range')[1]
            except:
                startpoint=int(0.05*polymer.getparvalue(ie,'BS'))
                endpoint = int(0.1*polymer.getparvalue(ie,'BS'))
            phi = get_mag_amplitude(fid, startpoint, endpoint,
                                    polymer.getparvalue(ie,'NBLK'),
                                    polymer.getparvalue(ie,'BS'))
            df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
            df['phi_normalized'] = (df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] )
            polymer.addparameter(ie,'df_magnetization',df)
            fit_option = 2 #mono exponential, 3 parameter fit
            p0=[1.0, polymer.getparvalue(ie,'T1MX')**-1*2, 0]
            df, popt = magnetization_fit(df, p0, fit_option)
            source_df.data = ColumnDataSource.from_df(df)
            
            polymer.addparameter(ie,'popt(mono_exp)',popt)
            print(popt)

            #print(df)
            #print(polymer.getparvalue(ie,'df_magnetization'))
            fid_slider = RangeSlider(start=1,end=polymer.getparvalue(ie,'BS'),range=(startpoint,endpoint),step=1,callback_policy='mouseup')
            tabs.tabs[0].child.children[2]=fid_slider            

        except KeyError:
            print('no relaxation experiment found')
            tau=np.zeros(1)
            phi=np.zeros(1)
            df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
            df['phi_normalized'] = np.zeros(1)
            df['fit_phi'] = np.zeros(1)
            source_df.data = ColumnDataSource.from_df(df)

    
    #this source is only used to communicate to the actual callback (cb)
    source = ColumnDataSource(data=dict(value=[]))
    source.on_change('data',cb)
    
    slider = Slider(start=1, end=nr_experiments, value=1, step=1,callback_policy='mouseup')
    slider.callback=CustomJS(args=dict(source=source),code="""
        source.data = { value: [cb_obj.value] }
    """)#unfortunately this customjs is needed to throttle the callback in current version of bokeh

    def calculate_mag_dec(attr,old,new):
        ie=slider.value
        polymer.addparameter(ie,'fid_range',new['range'])
        print(polymer.getparvalue(ie,'fid_range')) #this works
        start = new['range'][0]
        stop = new['range'][1]
        fid=polymer.getfid(ie)
        tau = polymer.getparvalue(ie,'df_magnetization').tau
        phi = get_mag_amplitude(fid, start, stop,
                                    polymer.getparvalue(ie,'NBLK'),
                                    polymer.getparvalue(ie,'BS'))

        df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
        df['phi_normalized'] = (df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] )

        fit_option = 2 #mono exponential, 3 parameter fit
        p0=polymer.getparvalue(ie,'popt(mono_exp)')
        df, popt = magnetization_fit(df, p0, fit_option)
        source_df.data = ColumnDataSource.from_df(df)
        polymer.addparameter(ie,'df_magnetization',df)
        polymer.addparameter(ie,'popt(mono_exp)',popt)
        pass
    
    source2 = ColumnDataSource(data=dict(range=[], ie=[]))
    source2.on_change('data',calculate_mag_dec)
    fid_slider.callback=CustomJS(args=dict(source=source2),code="""
        source.data = { range: cb_obj.range }
    """)#unfortunately this customjs is needed to throttle the callback in current version of bokeh

    
    def update_parameters():
        par_df, columns, discrete, continuous, time, quantileable = polymer.scan_parameters()
        x.options=columns
        y.options=columns
        size.options=['None']+quantileable
        color.options=['None']+quantileable
    # select boxes for p4
    button_scan = Button(label='Scan Parameters',type="success")
    button_scan.on_click(update_parameters)
    
    x = Select(title='X-Axis', value='ZONE', options=columns)
    x.on_change('value', update)

    y = Select(title='Y-Axis', value='TIME', options=columns)
    y.on_change('value', update)

    size = Select(title='Size', value='None', options=['None'] + quantileable)
    size.on_change('value', update)

    color = Select(title='Color', value='None', options=['None'] + quantileable)
    color.on_change('value', update)




    controls_p4 = widgetbox([button_scan, x,y,color,size], width=150)
    layout_p4 = row(controls_p4,plot_par())

    
    #fitting on all experiments
    p3 = figure(plot_width=300, plot_height=300,
            title='normalized phi vs normalized tau', webgl = True,
                y_axis_type = 'log',
                x_axis_type = 'linear')
    
    #fit magnetization decay for all experiments
    r1=np.zeros(nr_experiments)
    MANY_COLORS = 0
    p3_line_glyph=[]
    for i in range(nr_experiments):
        try:
            par=polymer.getparameter(i)
            fid=polymer.getfid(i)
            tau= get_x_axis(polymer.getparameter(i))
            startpoint=int(0.05*polymer.getparameter(i)['BS'])
            endpoint=int(0.1*polymer.getparameter(i)['BS'])
            phi = get_mag_amplitude(fid, startpoint, endpoint,
                                    polymer.getparameter(i)['NBLK'], polymer.getparameter(i)['BS'])
            df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
            df['phi_normalized']=(df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] )
            polymer.addparameter(i,'df_magnetization',df)

            p0 = [1, 2 * polymer.getparvalue(i,'T1MX'), 0]
            df, popt = magnetization_fit(df,p0, fit_option=2)
            polymer.addparameter(i,'amp',popt[0])
            polymer.addparameter(i,'r1',popt[1])
            polymer.addparameter(i,'noise',popt[2])
            r1[i] = popt[1]
            tau = popt[1]*df.tau
            phi = popt[0]**-1*(df.phi_normalized - popt[2])
            p3_df=pd.DataFrame(data=np.c_[ tau, phi ], columns=['tau', 'phi'])
            source_p3=ColumnDataSource(data=ColumnDataSource.from_df(p3_df))
            p3_line_glyph.append(p3.line('tau', 'phi', source=source_p3))
            MANY_COLORS+=1
        except KeyError:
            print('no relaxation experiment found')
    COLORS=viridis(MANY_COLORS)
    for ic in range(MANY_COLORS):
        p3_line_glyph[ic].glyph.line_color=COLORS[ic]
    par_df['r1']=r1
    layout_p1 = column(slider, p1,fid_slider, p2, p3)
    tab_relaxation = Panel(child = layout_p1, title = 'Relaxation')
    tab_parameters = Panel(child = layout_p4, title = 'Parameters')
    tabs = Tabs(tabs = [tab_relaxation, tab_parameters])
    #print(tabs.tabs[1].child.children[1]) # to access figure in parameters
    doc.add_root(tabs)
    #doc.add_root(layout_p1)
    #doc.add_root(layout_p4)
    doc.add_root(source) # i need to add source to detect changes
    doc.add_root(source2)


bokeh_app = Application(FunctionHandler(modify_doc))

server = Server({'/': bokeh_app}, io_loop=io_loop)
server.start()

if __name__ == '__main__':
    print('Opening Bokeh application on http://localhost:5006/')

    io_loop.add_callback(server.show, "/")
io_loop.start()
