#!/usr/bin/python
# -*- coding: UTF-8 -*-
import os
import pandas as pd
import numpy as np
from bokeh.models.widgets import Panel, Tabs
from bokeh.plotting import figure, ColumnDataSource
from bokeh.palettes import Spectral5, viridis  # @UnresolvedImport
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Select, Button, DataTable, TableColumn, Slider, RangeSlider, TextInput#, ColumnDataSource
from bokeh.models.callbacks import CustomJS
from tornado.ioloop import IOLoop
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.layouts import column, row
#from bokeh.models import ColumnDataSource#, Slider, RangeSlider, TextInput
from bokeh.server.server import Server
import stelardatafile as sdf
from utils import model_exp_dec, fun_exp_dec, get_mag_amplitude, magnetization_fit
from scipy.optimize import leastsq

# in the plot 4 use following impo
SIZES = list(range(6, 22, 3)) # for some sizes
COLORS = Spectral5 # for some colors (more colors would be nice somehow)


def create_plot_1_and_2(source_fid, source_mag_dec):
    # create and plot figure 1
    p1 = figure(plot_width=800, plot_height=500,
                title='Free Induction Decay', webgl=True,
                lod_factor=1000,lod_interval=150,lod_threshold=1000)
    # convert data to handle in bokeh
    # ignore rephased_fid at this point
#     source_fid = ColumnDataSource(data=ColumnDataSource.from_df(rephased_fid))
    p1.line('index', 'im', source=source_fid, color='blue')
    p1.line('index', 'real', source=source_fid, color='green')
    p1.line('index', 'magnitude', source=source_fid, color='red')

    # create and plot figure 2
    p2 = figure(plot_width=400, plot_height=400,
                title='Magnetization Decay')
    p2.circle_cross('tau', 'phi_normalized', source=source_mag_dec, color="navy")
    p2.line('tau', 'fit_phi', source=source_mag_dec, color="teal")
    return p1, p2



def fit_mag_decay_all(polymer, par_df):
    p3 = figure(plot_width=400, plot_height=400,
            title='normalized phi vs normalized tau', webgl = True,
                y_axis_type = 'log',
                x_axis_type = 'linear')

    nr_experiments = polymer.get_number_of_experiments()
    r1=np.zeros(nr_experiments)
    MANY_COLORS = 0
    p3_line_glyph=[]
    for i in range(1, nr_experiments):
        try:
            par=polymer.getparameter(i)
            fid=polymer.getfid(i)
            # TODO: here tau axis is being calculated the 3rd time. Look if it can be simplified
            tau= polymer.get_tau_axis(i)
            try:
                startpoint=polymer.getparvalue(i,'fid_range')[0]
                endpoint=polymer.getparvalue(i,'fid_range')[1]
            except:
                startpoint=int(0.05*polymer.getparvalue(i,'BS'))
                endpoint = int(0.1*polymer.getparvalue(i,'BS'))
                polymer.addparameter(i,'fid_range',(startpoint,endpoint))
            phi = get_mag_amplitude(fid, startpoint, endpoint,
                                    polymer.getparameter(i)['NBLK'], polymer.getparameter(i)['BS'])
            df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
            df['phi_normalized']=(df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] )
            polymer.addparameter(i,'df_magnetization',df)
            
            # initial values for the fit
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
            polymer.addparameter(i,'amp',float('NaN'))
            polymer.addparameter(i,'r1',float('NaN'))
            polymer.addparameter(i,'noise',float('NaN'))
    COLORS=viridis(MANY_COLORS)
    for ic in range(MANY_COLORS):
        p3_line_glyph[ic].glyph.line_color=COLORS[ic]
    par_df['r1']=r1
    return p3

#TODO: work with more than one tab in the browser: file input, multiple data analysis and manipulation tabs, file output and recent results tab
#TODO: clean modify_doc, it is awfully crowded here...
def modify_doc(doc):
    def get_data_frames(ie,):
        fid = polymer.getfid(ie) #read FID or series of FIDs for selected experiment
        try:
            # TODO: merge this with the initial plot. Now its doubled in the code
            tau = polymer.get_tau_axis(ie) #numpy array containing the taus for experiment ie
            try:
                startpoint=fid_slider.range[0] #lower integration bound
                endpoint = fid_slider.range[1] #upper integration bound
            except NameError:
                # fid_slider not initialized for first plot. Use default values:
                startpoint=int(0.05*polymer.getparvalue(ie,'BS'))
                endpoint = int(0.1*polymer.getparvalue(ie,'BS'))
                
            polymer.addparameter(ie,'fid_range',(startpoint,endpoint)) #add integration range to parameters to make it accesible
            phi = get_mag_amplitude(fid, startpoint, endpoint,
                                    polymer.getparvalue(ie,'NBLK'),
                                    polymer.getparvalue(ie,'BS')) # list containing averaged fid amplitudes (which is proportional to a magnetization phi)
            df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi']) # DataFrames are nice
            df['phi_normalized'] = (df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] ) #Normalize magnetization,
            #Note: in the normalized magnetization the magnetization build-up curves and magnetization decay curves look alike
            #Note: this makes it easier for fitting as everything looks like 1 * exp(-R/time) in first order
            polymer.addparameter(ie,'df_magnetization',df) # make the magnetization dataframe accesible as parameter
            fit_option = 2 #mono exponential, 3 parameter fit
            p0=[1.0, polymer.getparvalue(ie,'T1MX')**-1*2, 0] #choose startparameters for fitting an exponential decay
            df, popt = magnetization_fit(df, p0, fit_option) # use leastsq to find optimal parameters
            polymer.addparameter(ie,'popt(mono_exp)',popt) # add fitting parameters for later access
            print(popt) # print the fitting parameters to console (for convenience)
        except KeyError:
            print('no relaxation experiment found')
            tau=np.zeros(1)
            phi=np.zeros(1)
            df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
            df['phi_normalized'] = np.zeros(1)
            df['fit_phi'] = np.zeros(1)
        return fid, df

    def calculate_mag_dec(attr, old, new, start_ie=None):
        ''' Is being call from the callback for the experiment chooser
            loads selected experiment visualize in plot p1 and p2 
            gets experiment number from the slider
            writes source_fid.data from the fid from the polymer object
            writes source_mag_dec.data from the dataframe
            '''
        ie = experiment_slider.value   #get expermient number from the slider
        fid, df = get_data_frames(ie)
        source_fid.data=ColumnDataSource.from_df(fid) #convert fid to bokeh format
        source_mag_dec.data = ColumnDataSource.from_df(df)

    def plot_par():
        ''' Creates plot for the parameters '''
        # TODO: what is x here
        xs = par_df[x.value ].values
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
        ''' Callback for update of figure 1 in parameters tab '''
        tabs.tabs[1].child.children[1] = plot_par()

    def experiment_update(attr, old, new):
        """ Callback for the experiment chooser
        """
        ie = experiment_slider.value
        fid_slider.end = polymer.getparvalue(ie,'BS')
        try:
            fid_slider.range=polymer.getparvalue(ie,'fid_range')
        except:
            startpoint = int(0.05 * polymer.getparvalue(ie,'BS'))
            endpoint = int(0.1 * polymer.getparvalue(ie,'BS'))
            fid_slider.range=(startpoint,endpoint)
        calculate_mag_dec(attr,old,new)
        
    def table_update(attr,old,new):
        ''' callback for loading of data '''
        print('Hello table')
        path=pathbox.value.strip()
        file=filebox.value.strip()
        if file=="*.sdf":
            print('hello IF')
            allsdf=filter(lambda x: x.endswith('.sdf'),os.listdir(path))
            for f in allsdf:
                sdf_list.append(sdf.StelarDataFile(f,path))
        else:
            sdf_list.append(sdf.StelarDataFile(file,path))
        
        filenames=[x.file() for x in sdf_list]
        filenames_df=pd.DataFrame(data=filenames,columns=['file'])
        table_source.data=ColumnDataSource.from_df(filenames_df)
        
    def update_parameters():
        ''' callback for button
            function to call when button is clicked
            for updates parameters of polymer, since they can change during evaluation '''
        par_df, columns, discrete, continuous, time, quantileable = polymer.scan_parameters()
        x.options=columns
        y.options=columns
        size.options=['None']+quantileable
        color.options=['None']+quantileable

    ### This is the start of the script ###
    ### The callbacks are above ###

    #load data:
    polymer = load_data()
    nr_experiments = polymer.get_number_of_experiments()
    # initially set ie = 1
    start_ie = 1
    par_df, columns, discrete, continuous, time, quantileable = polymer.scan_parameters(20)
    # for the initial call get the dataframes without callback
    # they are being updated in following callbacks
    fid, df = get_data_frames(start_ie)
    source_fid = ColumnDataSource(data=ColumnDataSource.from_df(fid))
    source_mag_dec = ColumnDataSource(data=ColumnDataSource.from_df(df))
    # initialy creates the plots p1 and p2
    p1, p2 = create_plot_1_and_2(source_fid, source_mag_dec)
    
    ### initiates widgets, which will call the callback on change ###
    # initiate slider to choose experiment
    experiment_slider = Slider(start=1, end=nr_experiments, value=1, step=1,callback_policy='mouseup', width=800) #select experiment by value
    # initiate slider for the range in which fid shall be calculated
    # select the intervall from which magneitization is calculated from fid
    fid_slider = RangeSlider(start=1,end=polymer.getparvalue(start_ie,'BS'),
                             range=polymer.getparvalue(start_ie,'fid_range'),
                             step=1,callback_policy='mouseup', width=800)

    # initialize empty source for experiment slider
    source = ColumnDataSource(data=dict(value=[]))
    # 'data' is the attribute. it's a field in source, which is a ColumnDataSource
    # initiate experiment_update which is the callback
    source.on_change('data',experiment_update) #source for experiment_slider
    experiment_slider.callback = CustomJS(args=dict(source=source),code="""
        source.data = { value: [cb_obj.value] }
    """)#unfortunately this customjs is needed to throttle the callback in current version of bokeh

    # initialize empty source for fid slider, same as above
    source2 = ColumnDataSource(data=dict(range=[], ie=[]))
    source2.on_change('data',calculate_mag_dec)
    fid_slider.callback=CustomJS(args=dict(source=source2),code="""
        source.data = { range: cb_obj.range }
    """)#unfortunately this customjs is needed to throttle the callback in current version of bokeh

    # same for the update button
    button_scan = Button(label='Scan Parameters',button_type="success")
    button_scan.on_click(update_parameters)
    
    # here comes for callbacks for x, y, size, color
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


    # fit magnetization decay for all experiments
    p3 = fit_mag_decay_all(polymer, par_df)

    ####
    #### TODO: write file input
    #### TODO: select files to import
    #### TODO: discard imported files
    ####

    # load more data:
    table_source=ColumnDataSource(data=dict())
    sdf_list=[polymer]
    filenames=[x.file() for x in sdf_list]
    files_df=pd.DataFrame(data=filenames,columns=['file'])
    table_source.data=ColumnDataSource.from_df(files_df)
    t_columns = [
        TableColumn(field='file', title='Path / Filename'),
        #TableColumn(field='file', title='Filename'),
        ]
    table=DataTable(source=table_source,columns=t_columns)
    pathbox=TextInput(title="Path",value=os.path.curdir)
    filebox=TextInput(title="Filename",value="*.sdf")
    pathbox.on_change('value',table_update)
    filebox.on_change('value',table_update)
    layout_input=column(pathbox,filebox,table)  

    # set the layout of the tabs
    layout_p1 = column(experiment_slider, p1,row(p2, p3), fid_slider)
    tab_relaxation = Panel(child = layout_p1, title = 'Relaxation')
    tab_parameters = Panel(child = layout_p4, title = 'Parameters')
    tab_input = Panel(child = layout_input, title = 'Data In')

    # initialize tabs object with 3 tabs
    tabs = Tabs(tabs = [tab_relaxation, tab_parameters, tab_input])

    doc.add_root(tabs)
    doc.add_root(source) # i need to add source to detect changes
    doc.add_root(source2)

def load_data():
    #specify and import data file
    path=os.path.join(os.path.curdir,'data')
    # initialize StelarDataFile Object
    polymer=sdf.StelarDataFile('glyzerin_d3_300K.sdf',path)
    polymer.sdfimport()
    nr_experiments = polymer.get_number_of_experiments()
    polymer.rephase_fids()
    return polymer

def main():
    bokeh_app = Application(FunctionHandler(modify_doc))
    io_loop = IOLoop.current()
    print('Opening Bokeh application on http://localhost:5006/')
    server = Server({'/': bokeh_app}, io_loop=io_loop)
    server.start()
    io_loop.add_callback(server.show, "/")
    io_loop.start()

if __name__ == '__main__':
    main()

