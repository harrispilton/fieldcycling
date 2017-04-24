class FitModel:
    ###This class should make it easier to add models together so you can fit your data with flexibility
    ###ported from matlab program "NMRD profiler" by Lionel Broche
    def __init__(self): #brauch ich das überhaupt??
        self.modelName='' #string, name of the model (e.g. used for legend)
        self.modelEquation='' #string, put your equation here.
        self.variableName='' #in case you want more than one variable
        self.parameterName=[] #list of strings, name of the parameters appearing in the equation
        self.isFixed=[] #list of booleans, set to 1 if corresponding parameter is fixed, 0 if they are to be optimized
        self.minValue=[] #list of values, minimum values reachable for each parameter, put to -inf for no lower bound
        self.maxValue=[] #list of values, maximum values reachable for each parameter, put to inf for no upper bound
        self.startPoint=[] #list of values, start values for each parameter
        self.bestValue=[] #list of values, estimate from fitting routine
        self.errorBar=[] #list with confidence intervals !!Note: bokeh doesn't support errorBar plots yet. so no need to be specific at the moment #TODO: come up with a good solution
        self.fitobj=object()#fitting object created after the model is used for fitting TODO: is object() OK in this context?
