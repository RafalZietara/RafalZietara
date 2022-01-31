import math
from re import I

import matplotlib.pyplot as plt
import numpy as np
# from numpy.lib import histograms
from scipy import interpolate, optimize
import scipy.integrate as odeInt
from scipy.optimize.minpack import fsolve
import pint as UOM
# from sklearn.cluster import k_means
# from sqlalchemy import false, true
import shellthermo
import os.path
import CoolProp.CoolProp as CP
# http://www.coolprop.org/coolprop/wrappers/Python/index.html

# print("Importing Shellthermo version = %s" % sppts.version)

# path to sptts experimental database
# pathToExpSPPTSDB = None # required for oH2 and pH2
pathToExpSPPTSDB = 'C:\Apps\SPPTS for Python\Pepper\experimental.thermodb' # required for oH2 and pH2


mycol = ['red', 'green', 'blue', 'black','Magenta','Cyan']

class VSCdoc:
    def shortcuts(self):
        '''
        ## folding/unfolding:
        Ctrl + K + 0: fold all levels (namespace, class, method, and block)
        Ctrl + K + 1: namspace
        Ctrl + K + 2: class
        Ctrl + K + 3: methods
        Ctrl + K + 4: blocks
        Ctrl + K + [ or Ctrl + k + ]: current cursor block
        Ctrl + K + j: UnFold  
        ## comment out a section or uncomment
        Ctrl + K + C: comment
        Ctrl + K + U: un-comment
        '''
        pass

# several useful functions
def clip(number,minLim,maxLim):
    # print("Clip: ",minLim, number, maxLim)
    number = np.asarray(number)
    return max(min(number,maxLim),minLim)

def Range(start, end, step =1 ):
    return range(start, end+step, step)

def rError(Exception):
    """Error catching exception. Add own description to it"""
    def __init__(self, value):
        """Initilaise class."""
        self.value = value
    def __str__(self):
        """Create string."""
        return repr(self.value)

def intEx(xval, xref, yref):
    '''funciton to interpolate and extrapolate 
    note: numpy.interp does a bad job when it comes to extrapolation. wrong results even with linear interpolation!
    '''
    f = interpolate.interp1d(xref,yref,fill_value = "extrapolate")
    return f(xval)

class primitive:
    """
    Primitive class for all objects.

    # Info:
    name - name of the matter
    type - class type 
    description - description of the class
    diag - diagnostics level (see below)   
    
    # Verbosity levels (guidance):
    0 - no or little info
    1 - main parameters
    2 - some parameters
    3 - a lot of parameters
    4 - all parameters
    5 - main figures
    7 - all possible figures
    8 - write to file
    9 - write all to file
    """
    def __init__(self, name, type):
        self.name = name
        self.type = type
        self.description = "Primitive class for matter objects to be inherited e.g. fluid or metals"
        self.diag = 0

    def writeDetails(self):
        """Write correlation details into output file."""
        string = 'matter Name={0}; of the type: {1}'.format(self.name, self.type)
        if self.diag>0:
            return string
        else:
            [string + '\nDescription: {0}'.format(self.description)]
            return string

    def setDiag(self,verbosityLevel): # set diagnostics
        self.diag = verbosityLevel

class matter(primitive):
    """
    Base class for matter objects fluids or metals. 
    The class is also used for materials with constant all parameters.
    Derived classes give access to properties either via own correlations or connection to thermo package.
    
    # info:
    name - name of the material
    type - Matter: unknown, metal, fluid, solids (e.g.catalyst)
    source - const, correlation, CoolProp, SPTTS
        None - waiting to be specified
        const - constant data delivered via setProp function
        correlation - derived class taking data from own correlations 
        CoolProp - data from CoolProp
        SPPTS - data from Shell thermo package SPPTS 
    compName - names of components
    # Process variables:
    rho - [kg/m3] density
    Cp - [J/kg/K] heat capacity
    k - [W/m/K] Thermal conductivity of the heat exchanger material
    """
    def __init__(self, name="Material",type="unknown",DataSource=None):
        primitive.__init__(self,name,type)
        self.source = DataSource
        
    def setConstProps(self,Density,ThermalConducitivity,HeatCapacity):
        """if data are delivered as constant values"""
        self.__rho = Density # [kg/m3] density
        self.__Cp = HeatCapacity # [J/kg/K] heat capacity
        self.__k = ThermalConducitivity # [W/m/K] Thermal conductivity of the heat exchanger material
        # self.check()
        self.source = "const"
    
    def setPropDB(self):
        '''placeholder for now'''
        pass

    def check(self):
        if (not self.k) or (not self.rho) or (not self.Cp):
            pass 
            # raise rError('{0} of type {1} with source {2} is not specified correctly, use setConstProps(Density,ThermalConducitivity,HeatCapacity) to define material'.format(self.name,self.type,self.source))
        print('Metal setup check: passed')

    def rho(self,T):
        return self.__rho * np.ones_like(T)
    def Cp(self,T):
        return self.__Cp * np.ones_like(T)
    def k(self,T):
        return self.__k * np.ones_like(T)

class chnl(primitive): 
    """
    Derrived class for channels containing mainly geometry.

    # Geometry
    chnlH # [m] channel height
    chnlW # [m] channel width
    chnlL # [m] channel length
    Area_NS # [m**2] area on NS face 
    Area_WE # [m**2] area on WE face   
    chnlL = 0           # channel Length = length of HX
    chnlW = 0           # channel Width = width of HX
    chnlH = 0           # channel Height
    # specific channel fields will be added in derived classes 

    """
    mat =  []           # material atteched either fluid or 
    def __init__(self,name):
        primitive.__init__(self,name,type="channelGeometry")
        chnlL = None
        chnlW = None
        chnlH = None

    def check(self):
        ''' self check whether all has been set'''
        if (not self.chnlL) or (not self.chnlW) or (not self.chnlH):
            raise rError('Channel missing parameter')
        if (self.chnlL<=0) or (self.chnlW<=0) or (self.chnlH<=0):
            raise rError('Negative dimensions for channel: {0}({1})'.format(self.name,self.type))
    
    def addToHX(self,HXdata):
        ''' update parameters from HX parameters '''
        self.chnlW = HXdata.GeomW
        self.chnlL = HXdata.GeomH
        # below definition is OK for constant mesh data
        self.Area_NS = self.chnlH * self.chnlW
        self.Area_WE = self.chnlL * self.chnlW / HXdata.lyr.nCells
        # self.check()

    def setGeometry(self,channelThickness):
        self.chnlH = channelThickness / 1000      # mm -> m

    def addMaterial(self,F,isHot): 
        pass

class plate(matter, chnl):
    """
    (Derived) Class for solid plates separating channels with hydrogen and refrigerant 
    # info:
    description - description of the material

    # Geometry:
    chnlH # [mm] channel height
    epslon - roughness of the plate
    
    """
    def __init__(self,name):
        matter.__init__(self,name,type="metal")
        self.description = "Solid Metal Plate separating channels with fluid"
        self.epsilon = None

    def check(self):
        # chnl.check(self)
        if (not self.epsilon):
            raise rError('Use: setPlate(Height,epsilon) Roughness of metal')
        elif (self.epsilon<=0):
            raise rError('Material roughness must be >0 for: {0}(specified {1} mm)'.format(self.name,self.epsilon*1000))            
        else:
            pass

    def setPlate(self,Height,epsilon):
        chnl.setGeometry(self,Height)
        self.epsilon = epsilon / 1000   # mm -> m
        self.check()

class H2Ofluid(matter):
    """
    (Derived) Class for accessing fluid propeties
    # Info:
    source - source of process data
    """
    Tref = 400 # reference temperature
    def __init__(self, name="Material",type="fluid"):
        matter.__init__(self, name, type)
        self.description = "Class for fluids with overwritten and constant data handles numpy vectors as input for parameters"
        print("Stream {0} is set up with *fluidH2OConst* package as a {1}".format(name, type))

    def getProps(self,state,props):
        def default(): # submethod defining a wrong dictionary entry
            rError("Wrong property name when calling H2Ofluid")
        
        FlashType = state[0].lower()
        Param1, Param2, Param3 = np.asarray(state[1]), np.asarray(state[2]), np.asarray(state[3])
        nStates, nrProps = np.size(Param1), len(props)
        Result = np.zeros((nrProps,nStates), dtype=np.double)
        for i in range(nStates):
            if Param1.size == 1: # special case when only one state is considered 
                P1, P2, cmp = np.double(Param1), np.double(Param2), Param3
            else:
                P1, P2, cmp = Param1[i], Param2[i], Param3[i]
            try:
                propDict = {
                    'massDens' : self.massDens, # kg/m3
                    'Temp' : self.Temp ,
                    'ThermCond' : self.ThermCond,
                    # 'Pr' : self.PrNr,
                    'EnthMass' : self.Enthalpy
                    # molarWeight # kg/kmol
                    # molar Enthalpy J/mol
                    # specificEnthalpy J/kg
                    }
                for k in range(len(props)):
                    Result[k,i] = propDict.get(props[k],default)(P1,P2,cmp)
                    X=1
            except:
                rError('Invalid fluid specified.')
        return Result


    def massDens(self,T, P, C): # Temperature, Pressure, Composition
        return 1000 #kg/m3

    def PrNr(self,T,P,C): # Temperature, Pressure, Composition
        Result = 0
        mu = self.stream.viscosity()
        cp = self.stream.specificEnthalpy()
        k = self.ThermConduct()
        return mu*Cp/k

    def heatCapacity(self, T, P , C): # Temperature, Pressure, Composition
        return 4180.0 # J/kg/K heat capacity of fluid (Water)

    def Enthalpy(self,T, P, C): # Temperature, Pressure, Composition
        return 4180.0*(T-self.Tref) # J/kg

    def ThermCond(self,T, P, C): # Temperature, Pressure, Composition
        return 0.600 # W/K

    def Temp(self,P,H,C):
        return H/4180.0 + self.Tref
    
    def setDataSource(self, DataSource):
        """set what the source of data:
        const - constant data delivered via setMetal function
        correlation - derived class taking data from own correlations 
        CoolProp - data from CoolProp
        SPPTS - data from Shell thermo package SPPTS 
        """
        self.source = DataSource 

    def setPropDB(self, compNames, model = None):
        compDict = {
            'H2' :'H2',
            'oH2' :'oH2',
            'pH2' :'pH2',
            'He' : 'He',
            'H2O' : 'H2O',
            'N2' : 'N2',
            '1P' : '1P',
            '2P' :'2P',
            '3P' :'3P',
            '4P1': '4P1', # ibutane
            '4P': '4P' # nbutane
        }
        self.compNames = [compDict.get(k) for k in compNames]


class SPTTSFluid(matter):
    """
    Class delivering fluid properties from SPPTS package 
   
    # info:
    name - name of the material
    type - Matter: unknown, metal, fluid, solids (e.g.catalyst)
    source: SPTTS

    # Process variables:

    """
    def __init__(self, name="Material",type="fluid"):
        matter.__init__(self,name,type)
        self.source = "SPTTS {0}".format(shellthermo.version)
        print(self.source) 

    def setPropDB(self,compNames, model = 'CPA/SMIRK'):
        '''placeholder for now
        The density predictions of CPA are generally considered of low quality. That's why we usually use '.withModel(“CPA/SMIRK”)', 
        which combines CPA for the phase equilibrium properties and SMIRK for the density properties. Note that this can lead to 
        issues near the phase boundary, where CPA and SMIRK may not agree on what phases exist.
        As to oH2 and pH2, they are included in the 'experimental' database. 
        To use that instead of the default one, use '.withDatabase(“path/to/experimental.thermodb”)'. 
        The file should be in the Pepper folder in the package I sent you. The component IDs are 'oH2' and 'pH2'. Use them with the SMIRK model, 
        I don't think CPA parameters have been added yet, but SMIRK parameters have.
        '''
        compDict = {
            'H2' :'H2',
            'oH2' :'oH2',
            'pH2' :'pH2',
            'He' : 'He',
            'H2O' : 'H2O',
            'N2' : 'N2',
            '1P' : '1P',
            '2P' :'2P',
            '3P' :'3P',
            '4P1': '4P1', # ibutane
            '4P': '4P' # nbutane
        }
        self.compNames = [compDict.get(k) for k in compNames]
        if pathToExpSPPTSDB != None and os.path.exists(pathToExpSPPTSDB):
            self.solver = shellthermo.EquilibriumSolverBuilder().withDatabase(pathToExpSPPTSDB).withComponents(self.compNames).withModel(model).withPhases("VLL").build()
        else:
            self.solver = shellthermo.EquilibriumSolverBuilder().withComponents(self.compNames).withModel(model).withPhases("VLL").build()

    def check(self):
        if (not self.solver) or (not self.compNames):
            raise rError('Solver and component names not specified correctly for SPTTS, use ???')
        print('Fluid setup check: passed')
   
    def getProps(self,state,props):
        " get more than one input from Shell Thermo"
        def default(): # submethod defining a wrong dictionary entry
            rError("Wrong property name when calling SPTTS")

        FlashType = state[0].lower()
        Param1, Param2, Param3 = np.asarray(state[1]), np.asarray(state[2]), np.asarray(state[3])
        nStates, nrProps = np.size(Param1), len(props)
        Result = np.zeros((nrProps,nStates), dtype=np.double)
        for i in range(nStates):
            if Param1.size == 1: # special case when only one state is considered 
                P1, P2, cmp = np.double(Param1), np.double(Param2), Param3
            else:
                P1, P2, cmp = Param1[i], Param2[i], Param3[i]
            try:
                if FlashType == 'tpflash':
                    self.stream = self.solver.isothermal(P1, P2, cmp)
                elif FlashType == 'phflash': # different types of flashes. see documentation
                    self.stream = self.solver.isenthalpic_T(P1, P2, cmp) # pressure , enthalpy
                else:
                    rError("Flash type {0} is not supported currently".format(FlashType))

                propDict = {
                    'massDens' : self.stream.density, # kg/m3
                    'Temp' : self.stream.T ,
                    'ThermCond' : self.ThermConduct,
                    # 'Pr' : self.PrNr,
                    'EnthMass' : self.stream.enthalpy
                    # molarWeight # kg/kmol
                    # molar Enthalpy J/mol
                    # specificEnthalpy J/kg
                    }       
                for k in range(len(props)):
                    Result[k,i] = propDict.get(props[k],default)()
                    X=1
            except:
                rError('Invalid fluid specified.')
        return Result

    def ThermConduct(self):
        '''Process multiphase output to bulk property
        Note: This function assumes that flash is calculated.'''
        Result = 0
        st = self.stream
        for i in range(st.numberOfPhasesPossible()):
            if st.isPhasePresent(i):
                Result += st.phase(i).fraction()*st.phase(i).thermalConductivity()
         # return Result 
        print("careful: nasty overwrite in SPPTS")
        return 0.6

    def PrNr(self):
        '''This function assumes that flash is calculated'''
        Result = 0
        mu = self.stream.viscosity()
        cp = self.stream.specificEnthalpy()
        k = self.ThermConduct()
        return mu*Cp/k

class CoolPropFluid(matter):
    """
    Class delivering fluid properties from CoolProp package 
   
    # info:
    name - name of the material
    type - Matter: unknown, metal, fluid, solids (e.g.catalyst)
    solver - object with thermo
    source: CoolProp
    # Process variables:
    """

    def __init__(self, name="Material",type="fluid"):
        matter.__init__(self,name,type)
        self.source = "CoolProp ver. {0} (GIT: {1})".format(CP.get_global_param_string("version"), CP.get_global_param_string("gitrevision"))
        print(self.source)
        # TODO: no binary pair for o and p hydrogen so can do with a mixing rule linear or Lorentz-Berthelot
        # CP.apply_simple_mixing_rule('OrthoHydrogen', 'ParaHydrogen','Lorentz-Berthelot')
        # CP.apply_simple_mixing_rule('Hydrogen', 'Hydrogen','Lorentz-Berthelot')


    def setPropDB(self,compNames, model = 'CoolProp'):
        '''placeholder for now'''
        
        # below overwrite because no binary mixing rules.
        compDict = {
            'H2' :'Hydrogen',
            'oH2' :'Hydrogen', # 'OrthoHydrogen',
            'pH2' :'Hydrogen', # 'ParaHydrogen',
            'H2O' : 'Water',
            'N2' : 'N2',
            '1P' : 'Methane',
            '2P' :'Ethane',
            '3P' :'Propane',
            '4P1': 'IsoButane', 
            '4P': 'n-Butane',
            'He': 'Helium'
        }
        # low level initalizaiton of the solver
        # self.CPComp = ""
        # for i in range(len(self.compNames)):
        #     self.CPComp += self.compNames[i] + "&"
        # self.CPComp = self.CPComp[:-1]
        # self.solver = CP.AbstractState('HEOS',self.CPComp)
        self.compNames = [compDict.get(k) for k in compNames]

    def check(self):
        if (not self.solver) or (not self.compNames):
            raise rError('Solver and component names not specified correctly for SPTTS, use ???')
        print('Fluid setup check: passed')

    def getPropsLL(self, state, props): # not working at the moment stopped working when proper compenents were added.
        """
        Low level fluid property call.
        Uses CoolProp low level interface for all supported fluids.

        Inputs:
        state - CoolProp low level syntax [InputPairId,Prop1Value,Prop2Value] e.g.
                state = [CP.PT_INPUTS,20e6,273.15+150].
        props - List of CooProp output properties in level syntax. 
                For example, props = [CP.iPrandtl,CP.iDmass]
        outputProps - Array containing desired output properties.
        Notes: Currently only supports pressure & temperature input pairs.
        """
        def default(): # submethod defining a wrong dictionary entry
            rError("Wrong property name when calling SPTTS")

        FlashType = state[0]
        Param1, Param2, comp = np.asarray(state[1]), np.asarray(state[2]), np.asarray(state[3])
        nStates, nrProps = np.size(Param1), len(props)
        Result = np.zeros((nrProps,nStates), dtype=np.double)
        FT = None
        # Props = [PropIds(k) for k in props]
        if FlashType == 'TPFlash':
            FT = CP.PT_INPUTS
        elif FlashType == 'PHFlash': # different types of flashes. see documentation
            FT = CP.HmassP_INPUTS
        else:
            rError("Flash type {0} is not supported currently".format(FlashType))
        propDict = {
            'massDens' : CP.iDmass ,
            'Temp' : CP.iT,
            'CPMass' : CP.iCpmass, 
            'ThermCond' : CP.iconductivity,
            # 'Pr' : CP.iPrandtl,
            'EnthMass' : CP.iHmass,
            'Visc' : CP.iviscosity 
            }     
        CPprops = [propDict.get(k) for k in props]

        for i in range(nStates):
            if Param1.size == 1: # spcecial case when only one state is considered 
                istate = [FT, np.double(Param1), np.double(Param2)]
            else:
                istate = [FT, np.double(Param1[i]), np.double(Param2[i])]
            try:    
                self.solver.set_mole_fractions([1,0])
                self.solver.update(*istate)
                # it just fails here. I do not know why
                outputProps = [self.solver.keyed_output(k) for k in CPprops]
                if len(outputProps) == 1:
                    Result =  outputProps[0]
                else:
                    Result = outputProps
            except:
                rError('Invalid fluid specified.')
        return Result

    def getProps(self, state, props):
        """
        High level fluid property call.
        Uses CoolProp low level interface for all supported fluids.

        Inputs:
        state - CoolProp low level syntax [InputPairId,Prop1Value,Prop2Value] e.g.
                state = [CP.PT_INPUTS,20e6,273.15+150].
        props - List of CooProp output properties in level syntax. 
                For example, props = [CP.iPrandtl,CP.iDmass]
        outputProps - Array containing desired output properties.
        Notes: Currently only supports pressure & temperature input pairs.
        """
        def default(): # submethod defining a wrong dictionary entry
            rError("Wrong property name when calling SPTTS")

        FlashType = state[0]
        Param1, Param2, comp = np.asarray(state[1]), np.asarray(state[2]), np.asarray(state[3])
        nStates, nrProps = np.size(Param1), len(props)
        Result = np.zeros((nrProps,nStates), dtype=np.double)
        FT = None
 
        propHL = {
            'massDens' : 'D' ,
            'Temp' : 'T',
            'CPMass' : 'C', 
            'ThermCond' : 'L',
            'Pr' : 'Prandtl',
            'EnthMass' : 'Hmass',
            'Visc' : 'V'
            }

        if FlashType == 'TPFlash':
            P = ['T', 'P']
        elif FlashType == 'PHFlash': # different types of flashes. see documentation
            P = ['P', 'H']
        else:
            rError("flash now known {0}".format(FlashType))

        if nStates == 1:
            CompInp = 'HEOS::'
            # CompInp = 'SRK::'
            # CompInp = 'PR::'
            for i in range(len(self.compNames)):
                CompInp += self.compNames[i] + "[" + str(comp[i]) + "]&"
            CompInp = CompInp[:-1]

            Par1 = np.double(Param1)
            Par2 = np.double(Param2)
            for k in range(nrProps):
                CPprops = propHL.get(props[k])
                Result = CP.PropsSI(CPprops,P[0],Par1,P[1],Par2,CompInp)
        else:
            for i in range(nStates):
                cmpnt = comp[i]
                CompInp = 'HEOS::'
                # CompInp = 'PR::'
                # CompInp = 'SRK::'
                for j in range(len(self.compNames)):
                    CompInp += str(self.compNames[j]) + "[" + str(cmpnt[j]) + "]&"
                CompInp = CompInp[:-1]
                Par1 = np.double(Param1[i])
                Par2 = np.double(Param2[i])
                for k in range(nrProps):
                    CPprops = propHL.get(props[k])
                    Result[k,i] = CP.PropsSI(CPprops,P[0],Par1,P[1],Par2,CompInp)

        return Result

class H2fluid:
    '''
    Fluid from literature for Hydrogen

    # info:
    name - name of the material
    type - Matter: unknown, metal, fluid, solids (e.g.catalyst)
    source: Literature
    # Process variables:
    '''
    def __init__(self, name="Hydrogen",type="H2fluid"):
        matter.__init__(self,name,type)
        self.source = "Literature review"

    def setPropDB(self,compNames, model = 'H2fluid'):
        '''placeholder for now'''
        compDict = {
            'H2' :'H',
            'oH2' :'oH2',
            'pH2' :'pH2',
            'H2O' : 'water',
            'N2' : 'N2',
            '1P' : 'Methane',
            '2P' :'Ethane',
            '3P' :'Propane',
            '4P1': 'IsoButane', 
            '4P': 'n-Butane' 
        }
        self.compNames = [compDict.get(k) for k in compNames]

    def getProps(self, state, props):
        """
        Inputs:
        props - List of CooProp output properties in level syntax. 
                For example, props = [CP.iPrandtl,CP.iDmass]
        """
        FlashType = state[0].lower()
        Param1, Param2, Param3 = np.asarray(state[1]), np.asarray(state[2]), np.asarray(state[3])
        nrProps = len(props)
        Result = np.zeros(nrProps, dtype=np.double)
        if Param1.size == 1: # special case when only one state is considered 
            P1, P2, cmp = np.double(Param1), np.double(Param2), Param3
        else:
            P1, P2, cmp = Param1[i], Param2[i], Param3[i]
        try:
            if FlashType == 'tpflash':
                self.stream = self.solver.isothermal(P1, P2, cmp)
            elif FlashType == 'phflash': # different types of flashes. see documentation
                self.stream = self.solver.isenthalpic_T(P1, P2, cmp) # pressure , enthalpy
            else:
                rError("Flash type {0} is not supported currently".format(FlashType))
            propDict = {
                'massDens' : self.Density_Leachman, # kg/m3
                'Temp' : self.Temp ,
                'ThermCond' : self.ThrmlCond,
                'Pr' : self.PrNr,
                'EnthMass' : self.Enthalpy
                # molarWeight # kg/kmol
                # molar Enthalpy J/mol
                # specificEnthalpy J/kg
                }
            for k in range(len(props)):
                Result[k] = propDict.get(props[k])(P1,P2,comp)
        except:
            rError('Invalid fluid specified.')
        return Result
   
    def getParaFrac(self, C):
        posH2, posH2 = self.compNames.index("pH2"), self.compNames.index("oH2")
        # TODO: write here something if the hydrogen was not found
        cH2 = (C[posPH2], C[posOH2])
        result = 0.25 
        if C[0]+C[1] < 0.01: # basically no hydorgen
            # leave with overwrite of para = 0.25
            print(">>> Assume para Fraction = 0.25")
        else:
            result = C[0]/sum(C)
        return result

    def heatCapacity(self, T, P, C):
        # The ideal gas heat capacity at const.P in Btu/lbmole/K;
        # Input: Temperature [K] paraFrac [0-1]
        # cp coefficients for para and ortho hydrogen
        paraH2 = getParaFrac(C)
        cp = np.array([[9.6421, -0.11799, 0.0029401, -0.000021391, 0.000000069883, -0.00000000010698, 6.2306E-14],\
                       [10.769, -0.055323, 0.00052767, -0.0000015156, 0.0000000014248, 0, 0]])
        tempV= np.array([1, T, T**2, T**3, T**4, T**5, T**6]) 
        # cp = (1- para) * cp_ortho + fracPara * cp_para
        Heat_Capacity_Ideal_Gas = (1 - paraH2 ) * sum(cp[1,:] * T) + paraH2 * sum(cp[0,:] * T)
        UOM = Scalar( 1 , 'BTU/lbmole/K')
        UOM.GetValue('kj/kmol/K')
        return Heat_Capacity_Ideal_Gas * UOM

    def HelmholtzEnergy(self,Temperature, Density):
        # TODO: Horrible function I wrote tat needs rewrite.
        # Helmholtz energy [J/mol] returns tuple [ortho, para] 
        # Temperature[K], density [lbm/ft^3] 
        # It is based upon the Leachman Equation of State as found in the masters thesis 
        # by Jacob Leachman at the Univ of Idaho, May 2007 
        # set the ratios to the critical conditions for the Leachman EOS 
        # delta is the ratio of the actual density to the critical density 
        # but first convert density to mol/liter because the critical 
        # density used is 15.538 mol/liter 
        
        emptyTen = [0,0,0,0,0,0,0,0,0,0] # dummy to get lists 10 to whatever
        deltao = Density * 3.2808**3 / 2.2046 / 2.0159 / 15.445 # 15.508 '15.445'
        deltap = Density * 3.2808**3 / 2.2046 / 2.0159 / 15.538 
        # tau is the ratio of the critical temperature to the actual temperature 
        tauo = 33.22 / Temperature # '33.145 / Temperature '33.22 / Temperature'
        taup = 32.938 / Temperature     
        # set the parameter arrays for ortho hydrogen 
        No = [0, -6.83148, 0.01, 2.11505, 4.38353 , 0.211292 , -1.00939 , 0.142086 , -0.87696 , 0.804927, -0.710775 , 0.0639688 , 0.0710858 , -0.087654 , 0.647088]
        Np = [0, -7.33375 , 0.01 , 2.60375 , 4.66279 , 0.68239 , -1.47078 , 0.135801 , -1.05327 , 0.328239 ,  -0.0577833 ,  0.0449743 ,  0.0703464 ,  -0.0401766 ,  0.11951] 
        to =  [0, 0.7333, 1, 1.1372, 0.5136, 0.5638 , 1.6248 , 1.829 , 2.404 , 2.105 , 4.1 , 7.658 , 1.259 , 7.589 , 3.946 ]
        tp = [0, 0.6855 , 1 , 1 , 0.489 , 0.774 , 1.133 , 1.386 , 1.619 , 1.162 ,  3.96 ,  5.276 ,  0.99 ,  6.791 ,  3.19]
        d_lower = [0, 1, 4 , 1 , 1 , 2 , 2 , 3 , 1 , 3 , 2 , 1 , 3 , 1 , 1]
        p  = [0, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 , 0 , 0 , 0 , 0 , 0]   
        phio = emptyTen + [-1.169 ,-0.894 ,-0.04 ,-2.072 ,-1.306]
        phip = emptyTen + [-1.7437 , -0.5516 , -0.0634 , -2.1341 , -1.7777]     
        betao = emptyTen + [-0.4555 , -0.4046 , -0.0869 , -0.4415 , -0.5743]
        betap = emptyTen +[-0.194, -0.2019 , -0.0301 , -0.2383 , -0.3253] 
        gammao = emptyTen + [1.5444 , 0.6627 , 0.763 , 0.6587 , 1.4327]
        gammap = emptyTen + [0.8048,  1.5248 ,  0.6648 ,  0.6832 ,  1.493]     
        D_UPPERo = emptyTen + [0.6366, 0.3876 , 0.9437 , 0.3976 , 0.9626 ]
        D_UPPERp = emptyTen + [1.5487 , 0.1785 , 1.28 , 0.6319 , 1.7104]    
        # the following constants are for ortho hydrogen and were taken from 
        # the revised ideal gas heat capacity function presented by Leachman 
        v_cpo = [0, 2.54151, -2.3661 , 1.00365 , 1.22447 , 0 , 0 , 0] 
        v_cpp = [0, 4.30256 , 13.0289 , -47.7365 , 50.0013 , -18.6261 , 0.993973 , 0.536078] 
        u_cpo = [0, 856 , 1444 , 2194 , 6968 , 999 , 999 , 999] # last 3 numbers = 0 not 999 but entered to avoid division by zero in integral_cp
        u_cpp = [0, 499 , 826.5 , 970.8 , 1166.2 , 1341.1 , 5395 , 10185] 
        # calculate the ideal gas Helmholtz energy parameter alpha_zero as alpha_zero=a0/R/T 
        # a0=h0-RT-Ts0 # ideal gas Helmholtz free energy= ideal gas enthalpy - RT - T * ideal gas  entropy 
        # h0=h0_ref+integral(cp_ideal_gas)dT from T0 to T 
        # s0=s0_ref+integral(cp_ideal_gas / T) dT from T0 to T - R ln (rho T /rho0/T0) 
        # the following uses the function constants from the heat capacity for ideal gas function 
        # assuming that the reference temperature is 300 K 
        # first calculate the integral of cp ideal gas at the current temperature then 
        # subtract the value at the reference temperature 
        # here the gas constant is 3.5748 Btu/lbmole/K 
        TK = Temperature 
        # initialize the value 
        integral_cp_ortho = 3.5748 * 2.5 * TK
        integral_cp_para = integral_cp_ortho
        
        for IJ in Range(1,7):
            integral_cp_ortho = integral_cp_ortho + 3.5748 * v_cpo[IJ] * u_cpo[IJ] / (math.exp(u_cpo[IJ] / TK) - 1)
            integral_cp_para = integral_cp_para + 3.5748 * v_cpp[IJ] * u_cpp[IJ] / (math.exp(u_cpp[IJ] / TK)- 1) 
    
        #now correct for the reference temperature 
        TK = 300 
        integral_cp_ortho = integral_cp_ortho - 3.5748 * 2.5 * TK
        integral_cp_para = integral_cp_para - 3.5748 * 2.5 * TK
    
        for IJ in Range(1,7):
            integral_cp_ortho = integral_cp_ortho - 3.5748 * v_cpo[IJ] * u_cpo[IJ] / (math.exp(u_cpo[IJ] / TK) - 1)
            integral_cp_para = integral_cp_para - 3.5748 * v_cpp[IJ] * u_cpp[IJ] / (math.exp(u_cpp[IJ] / TK) - 1) 
        # add the enthalpy (at 20 psia) of the reference temperature of 300 K in Btu/lbmole 
        # this is the ideal gas enthalpy in Btu/lbmole for the current temperature
        h0o = 3674.4 + integral_cp_ortho
        h0p = 3642.3 + integral_cp_para
    
        # now calculate the integral of cp ideal gas over T at the current temperature then subtract the value at the reference temperature 
        TK = Temperature 
        # 'intitialize the value 
        integral_cp_over_T_ortho = 3.5748 * 2.5 * math.log(TK)
        integral_cp_over_T_para = 3.5748 * 2.5 * math.log(TK)
           
        for IJ in Range(1,7):
            integral_cp_over_T_ortho = integral_cp_over_T_ortho + 3.5748 * u_cpo[IJ] * v_cpo[IJ] / TK * math.exp(u_cpo[IJ] / TK) / (math.exp(u_cpo[IJ] / TK) - 1) - 3.5748 * v_cpo[IJ] * math.log(math.exp(u_cpo[IJ] / TK) - 1)
            integral_cp_over_T_para = integral_cp_over_T_para + 3.5748 * u_cpp[IJ] * v_cpp[IJ] / TK * math.exp(u_cpp[IJ] / TK) / (math.exp(u_cpp[IJ] / TK) - 1) - 3.5748 * v_cpp[IJ] * math.log(math.exp(u_cpp[IJ] / TK) - 1)
    
        # 'now subtract the reference values 
        TK = 300 
        integral_cp_over_T_ortho = integral_cp_over_T_ortho - 3.5748 * 2.5 * math.log(TK)
        integral_cp_over_T_para = integral_cp_over_T_para - 3.5748 * 2.5 * math.log(TK)
        for IJ in Range(1,7):
            integral_cp_over_T_ortho = integral_cp_over_T_ortho - 3.5748 * u_cpo[IJ] * v_cpo[IJ] / TK * math.exp(u_cpo[IJ] / TK) / (math.exp(u_cpo[IJ] / TK) - 1)  + 3.5748 * v_cpo[IJ] * math.log(math.exp(u_cpo[IJ] / TK) - 1)
            integral_cp_over_T_para = integral_cp_over_T_para - 3.5748 * u_cpp[IJ] * v_cpp[IJ] / TK * math.exp(u_cpp[IJ] / TK) / (math.exp(u_cpp[IJ] / TK) - 1) + 3.5748 * v_cpp[IJ] * math.log(math.exp(u_cpp[IJ] / TK) - 1)
            
        # add the entropy of the reference temperature in Btu/lbmole/K and subtract the natural log term 
        # this is the ideal gas entropy in Btu/lbmole/K for the current temperature 
        # here the gas constant is 3.5748 Btu/lbmole/K 
        # here the ideal gas density at 20 psia and 300 K is included 
        s0o = 59.099 + integral_cp_over_T_ortho - 3.5748 * math.log(Density * Temperature / 0.006951 / 300)
        s0p = 55.078 + integral_cp_over_T_para - 3.5748 * math.log(Density * Temperature / 0.006951 / 300) 
        
        # now for the ideal gas Helmholtz energy 
        alpha_zeroo = (h0o - 3.5748 * Temperature - Temperature * s0o) / 3.5748 / Temperature 
        alpha_zerop = (h0p - 3.5748 * Temperature - Temperature * s0p) / 3.5748 / Temperature 
        
        # calculate the residual Helmholtz energy parameter alpha_residual in parts 
        # calculate the first part 
        alpha_resid_oneo, alpha_resid_onep = 0, 0 
    
        for I in Range(1,7):
            alpha_resid_oneo = alpha_resid_oneo + No[I] * deltao**d_lower[I] * tauo**to[I]
            alpha_resid_onep = alpha_resid_onep + Np[I] * deltap**d_lower[I] * taup**tp[I] 
            
        # calculate the second part
        alpha_resid_twoo, alpha_resid_twop = 0, 0 
     
        for I in Range(8,9):
            alpha_resid_twoo = alpha_resid_twoo + No[I] * deltao**d_lower[I] * tauo**to[I] * math.exp(-deltao**p[I])
            alpha_resid_twop = alpha_resid_twop + Np[I] * deltap**d_lower[I] * taup**tp[I] * math.exp(-deltap**p[I])
    
        # calculate the third part 
        alpha_resid_threeo = 0
        for I in Range(10,14):
            if abs(deltao - D_UPPERo[I]) < 0.0000000001 :
                break
            if abs(tauo - gammao[I]) < 0.0000000001:
                break
            alpha_resid_threeo = alpha_resid_threeo + No[I] * deltao**d_lower[I] * tauo**to[I] * math.exp(phio[I] * (deltao - D_UPPERo[I])**2 + betao[I] * (tauo - gammao[I])**2)
            
        # calculate the third part make sure the value is initialized at zero 
        alpha_resid_threep = 0
        for I in Range(10,14):
            if abs(deltap - D_UPPERp[I]) < 0.0000000001 :
                break
            if abs(taup - gammap[I]) < 0.0000000001:
                break
            alpha_resid_threep = alpha_resid_threep + Np[I] * deltap**d_lower[I] * taup**tp[I] * math.exp(phip[I] * (deltap - D_UPPERp[I])**2 + betap[I] * (taup - gammap[I])**2)
            # print(Np[I], deltap, d_lower[I], taup, tp[I], phip[I], D_UPPERp[I], betap[I], gammap[I])
            
        # 'add the three parts of the residual 
        alpha_residualo = alpha_resid_oneo + alpha_resid_twoo + alpha_resid_threeo
        alpha_residualp = alpha_resid_onep + alpha_resid_twop + alpha_resid_threep
        
        # 'add the ideal gas and residual parts and convert to Helmholtz energy 
        # 'which required multiplication by 
        # 'the gas constant (8.314 J/mol/K) and Temperature 
        HelmholtzOrtho = (alpha_zeroo + alpha_residualo) * 8.314 * Temperature
        HelmholtzPara = (alpha_zerop + alpha_residualp) * 8.314 * Temperature 
        return [HelmholtzOrtho,HelmholtzPara]
    
    def HelmholtzEnergyH2(self,Temperature,Density,paraFrac):
        HE = self.HelmholtzEnergy(Temperature,Density)
        return HE[0]*(1-paraFrac)+HE[1]*paraFrac
    
    def HelmholtzNormal(self,Temperature,Density):
        return self.HelmholtzEnergyH2(Temperature,Density,0.25)
    
    def Entropy(self,Temperature, Density, ParaFraction):
        # Entropy[Btu/lbmole/K]; Temperature[K]; Density in [lbm/ft^3]; paraH2 fraction 
        # using the partial derivative S=-(dA/dT)at constant volume set the temperature differential 
        
        dt = 0.001 
        # calculate the individual entropies
        HETminusdt = self.HelmholtzEnergy(Temperature - dt, Density)
        HETplusdt = self.HelmholtzEnergy(Temperature + dt, Density)
        entropy_ortho = -(HETplusdt[0] - HETminusdt[0]) / (2 * dt) 
        entropy_para = -(HETplusdt[1] - HETminusdt[1]) / (2 * dt) 
        #  calculate the mixture entropy and convert J/mol/K to Btu/lbmole/K 
        Entropy = (entropy_para * ParaFraction + entropy_ortho * (1 - ParaFraction)) * 0.43021 
        # Correct for mixing 
        # here the gas constant is 3.5748 Btu/lbmole/K 
        # correction=-R summation(xi*ln(xi)) 
        # limit the para fraction or the log calculation will have trouble 
        if ParaFraction == 1 :
            entropy_correction = -3.5748 * (ParaFraction * math.log(ParaFraction))
        elif ParaFraction == 0:
            entropy_correction = -3.5748 * ((1 - ParaFraction) * math.log(1 - ParaFraction)) 
        else :
            entropy_correction = -3.5748 * (ParaFraction * math.log(ParaFraction) + (1 - ParaFraction) * math.log(1 - ParaFraction))
            print (entropy_correction, ParaFraction, math.log(ParaFraction))
        # # return the entropy in Btu/lbmole/K 
        Entropy = Entropy + entropy_correction 
        return Entropy

    def VaporPressure(self,Temperature):
        # Antoine equation log10(p)=A-B/(T+C)
        # valid 21.01-32.27 K
        # after NIST van Itterbeek, Verbeke 1964
        if (Temperature < 21.01 or Temperature >32.27): 
            print('Warning: vapourPressure function runs with Temperature ouside validity range')  
        coef = [3.54314, 99.395, 7.726] #A,B,C
        return pow(10,coef[0]-coef[1]/(Temperature+coef[2]))

    def Density_Leachman(self,Temperature, Pressure, Comp):
        # Density of normal hydrogen [lbm/ft3] 
        # Temperature [K], Pressure [psia] 
        # the density is not highly sensitive to parafraction 
        # initializing the density guess if one is not made is a significant improvement 
        # as iterations in the liquid region will frequently find the wrong solution of the multiple possible 
        # this section makes a guess if one was not made
        rho = 1 # initial guess
        # use the current guess unless it is below or equal zero
        if Guess<=0:
            # no guess has been entered, make an estimate 
            # for the Temperature >=33.145 K, use the high range correlation 
            if Temperature >= 33.145:
                # set the coefficients based on the current pressure 
                a_3 = -5.3815E-15 * Pressure**2 - 0.000000000050004 * Pressure + 0.000000000061132 
                a_2 = 4.3033E-12 * Pressure**2 + 0.000000048544 * Pressure - 0.000000053475 
                a_1 = -0.00000000098557 * Pressure**2 - 0.000016423 * Pressure + 0.000014879 
                a_0 = 0.000000040261 * Pressure**2 + 0.0022571 * Pressure - 0.0013078 
                rho = a_3 * Temperature**3 + a_2 * Temperature**2 + a_1 * Temperature + a_0
            else:
                if Pressure > 185:
                    # the temperature is less than 33.145 K 
                    # if the pressure is greater than the critical pressure, guess a high value 
                    rho = 5 
                else: 
                    # temperature and pressure are in the possible vapor/liquid zone 
                    # estimate the saturated vapor and liquid densities at this temperature 
                    rho_sat_liq = -0.00029297 * Temperature**4 + 0.02936 * Temperature**3 - 1.099 * Temperature**2 + 18.105 * Temperature - 106.08 
                    rho_sat_vap = 0.00030338 * Temperature**4 - 0.030401 * Temperature**3 + 1.1378 * Temperature**2 - 18.795 * Temperature + 115.51 
                    # estimate the saturated vapor pressure 
                    p_sat = self.VaporPressure(Temperature) 
                    if Pressure > p_sat:
                        # the actual pressure is above the estimated vapor pressure at this temperature 
                        # estimate the density using one close to the saturated liquid density 
                        rho = rho_sat_liq 
                    elif Pressure < p_sat: 
                        # the pressure is less than the saturated vapor pressure 
                        # estimate the density at the saturated vapor value 
                        rho = rho_sat_vap 
                    else:
                        # the pressure is likely equal to the saturated vapor pressure 
                        # estimate the density using a starting guess that is the average of liq and vapor 
                        rho = (rho_sat_vap + rho_sat_liq) / 2 
                  # make a density estimate for the zone where liquid and vapor might coexist 
    
    
        # this iteration scheme will first try to solve it with Newton-Raphson because it is usually fast 
        # but if it takes too long, it will start over with the bisection search 
        # because there are probably two steep slopes that make it very difficult for Newton-Raphson 
        # --------------------------------------------------------------------------------- 
        # set the solution scheme 
        scheme = 2 
        # initialize loop counter 
        counter = 0 
        # set the loop limit 
        counter_limit = 800 
        # set the tolerance for the convergence 
        tolerance = 0.0001 
        # set the finite difference of the density 
        d_rho = 0.00001 
        # set the limits of rho 
        rho_high = 5 
        rho_low = 0.001
        rho_1 = (rho_low+rho_high)/2
        Press_Error_1 = (-(self.HelmholtzNormal(Temperature, rho_1 + d_rho) - self.HelmholtzNormal(Temperature, rho_1 - d_rho)) / (1 / (rho_1 + d_rho) - 1 / (rho_1 - d_rho)) * 1.15248)
        Press_Error_1 = (Press_Error_1 - Pressure) / (Pressure) * 100 
    
        # start the loop that adjusts the assumed density until the pressure calculated 
        # from the partial derivative of the Helmholtz energy is within tolerance of the set pressure
        LoopAgain = True
        while LoopAgain: 
            # increment the counter 
            counter = counter + 1 
            # if scheme=1, use the Newton Raphson method 
            if scheme == 1:
                # --------------------------------------------------------------------------- 
                # save the old values 
                rho_2 = rho_1 
                rho_1 = rho 
                if counter > 2:
                    # make a new guess with the Newton Raphson method
                    print(rho_1, rho_2, Press_Error_2, Press_Error_1)
                    Slope_Press_Error = (Press_Error_2 - Press_Error_1) / (rho_2 - rho_1) 
                    rho = rho - Press_Error / Slope_Press_Error 
                elif counter == 2:
                    # make a new guess if this is the second time around 
                    rho = rho * 0.999  
                # End If
                
                # dont let the density guess be negative - set it to a low number 
                # limit the density to less than 5 lbm/ft3 
                rho = clip(rho, rho_low,rho_high)
    
                # calculate the pressure from the partial derivative 
                # P=-(dA/d 1/rho)at constant temperature 
                # use the parameters for the normal hydrogen calculation as density is not a strong function of the ortho para compotision 
                P_psia = (-(self.HelmholtzNormal(Temperature, rho + d_rho) - self.HelmholtzNormal(Temperature, rho - d_rho)) / (1 / (rho + d_rho) - 1 / (rho - d_rho)) * 1.15248) 
                # calculate the pressure error in percent 
                Press_Error = (P_psia - Pressure) / (Pressure) * 100 
                # save the pressure errors 
                Press_Error_2 = Press_Error_1 
                Press_Error_1 = Press_Error 
            else:
                # --------------------------------------------------------------------------- 
                # now using bisection search (by halves) because the newton raphson method 
                # had trouble with a sharp slope 
                # if this is not the first time through, pick a new rho 
                if counter > 1: 
                    rho = (rho_high + rho_low) / 2
                # dont let the density guess be negative - set it to a low number 
                # limit the density to less than 5 lbm/ft3 
                rho = clip(rho, rho_low,rho_high)
                # calculate the pressure from the partial derivative 
                # P=-(dA/d 1/rho)at constant temperature 
                # use the parameters for the normal hydrogen calculation as density is not a strong function of the ortho para compotision 
                P_psia = (-(self.HelmholtzNormal(Temperature, rho + d_rho) - self.HelmholtzNormal(Temperature, rho - d_rho)) / (1 / (rho + d_rho) - 1 / (rho - d_rho)) * 1.15248) 
                # calculate the pressure error in percent 
                Press_Error = (P_psia - Pressure) / (Pressure) * 100 
                # reset one of the limits based on the results 
                if Press_Error > 0:
                    # reset the upper limit 
                    rho_high = rho 
                else:
                    # reset the lower limit 
                    rho_low = rho  
            # --------------------------------------------------------------------------- 
            # check how many iterations have been run 
            # if it is more than allowed with Newton-Raphson, 
            # restart with the bisection method 
            if counter > 10:
                scheme = 2 
            # loop until the counter gets too high or the error is small enough 
            LoopAgain = (counter < counter_limit and math.fabs(Press_Error) > tolerance)
        # check if it is beyond the tolerance
        #if math.abs(Press_Error) > tolerance:
        #    break 
        # save the value 
        Density_Leachman = rho
        return rho
        
    def Enthalpy(self,Temperature, Pressure, ParaFraction):
        # Enthalpy[J/mol] 
        # from H = PV - T^2*(d(A/T)/dT)constant V, N 
        # Temperature [K]; Pressure [psia]; para hydrogen molar fraction 
        # find the stream density in lbm/ft3 
        
        rho = self.Density_Leachman(Temperature, Pressure, 0) 
        # set the Temperature increment in Kelvin 
        dt = 0.001 
        # find the ortho and para hydrogen [J/mol] and combine them
        HEplusdt =  self.HelmholtzEnergy(Temperature + dt, rho)
        HEminusdt= self.HelmholtzEnergy(Temperature - dt, rho)
        enthalpy_ortho = -Temperature**2 * (HEplusdt[0] / (Temperature + dt) - HEminusdt[0] / (Temperature - dt)) / ((Temperature + dt) - (Temperature - dt)) 
        enthalpy_para = -Temperature**2 * (HEplusdt[1] / (Temperature + dt) - HEminusdt[1] / (Temperature - dt)) / ((Temperature + dt) - (Temperature - dt)) 
        Enthalpy = enthalpy_ortho * (1 - ParaFraction) + enthalpy_para * ParaFraction 
        # convert J/mol to Btu/lbmole 
        Enthalpy = Enthalpy * 0.43021 
        # now add the PV term in Btu/lbmole 
        Enthalpy = Enthalpy + Pressure / rho * 144 * 2.0159 / 777.65 
        Enthalpy /= 0.43021 
        return Enthalpy

    def EnthNP(self,T, P, pH2):
        # Enthalpy[J/mol] 
        # from H = PV - T^2*(d(A/T)/dT)constant V, N 
        # Temperature [K]; Pressure [psia]; para hydrogen molar fraction 
        # find the stream density in lbm/ft3 
        # TODO: rewrite this as it should be 
        if (T.size() == P.size()) and (T.size() == pH2.size()):
            Enth = T * 0
        else:
            rError("EnthNP function expects all input params to have the same size")

        for i in range(len(T)):
            rho = self.Density_Leachman(T[i], P[i], 0) 
            # set the Temperature increment in Kelvin 
            dt = 0.001 
            # find the ortho and para hydrogen [J/mol] and combine them
            HEplusdt =  self.HelmholtzEnergy(T[i] + dt, rho)
            HEminusdt= self.HelmholtzEnergy(T[i] - dt, rho)

            enthalpy_ortho = -T[i]**2 * (HEplusdt[0] / (T[i] + dt) - HEminusdt[0] / (T[i] - dt)) / ((T[i] + dt) - (T[i] - dt)) 
            enthalpy_para = -T[i]**2 * (HEplusdt[1] / (T[i] + dt) - HEminusdt[1] / (T[i] - dt)) / ((T[i] + dt) - (T[i] - dt)) 
            Enthalpy = enthalpy_ortho * (1 - pH[2]) + enthalpy_para * pH2[i] 
            #  convert J/mol to Btu/lbmole 
            Enthalpy = Enthalpy * 0.43021 
            # now add the PV term in Btu/lbmole 
            Enthalpy = Enthalpy + Pressure / rho * 144 * 2.0159 / 777.65 
            Enthalpy /= 0.43021
            Enth[i] = Enthalpy 
        return Ent

    def EnthHelp(self,T, P, pH2, H):
        return H - self.Enthalpy(T,P,pH2)

    def getTfromHP(self,Told,P,H,pH2):
        # I hate this function this should be from EOS. Not this way !!
        # bi-section way
        # maybe interpolation would be better
        PARM = (P, H, pH2,H)
        return optimize.fsolve(EnthHelp, Told, args=PARM, xtol=1e-8)
   
    def Thermal_Conductivity(self,Temperature, Pressure, ParaFraction, Density=0): 

        if (Density):
            rho = Density
        else: 
            rho = self.Density_Leachman(Temperature, Pressure, 0) 
    
        # set the reduced temperature with the critical temperature for para 
        T_red = Temperature / 32.938 
        # set the thermal conductivity for the para hydrogen part 
        # this is a five part equation so set the constant for each part as a function of pressure 
        # these were determined from a curve fit of thermal conductivity divided by the density 
        para_4 = -4.7484 * (1 / Pressure)**3 + 0.39601 * (1 / Pressure)**2 + 0.060825 * (1 / Pressure) - 0.0000169 
        para_3 = 138.2 * (1 / Pressure)**3 - 11.591 * (1 / Pressure)**2 - 1.6886 * (1 / Pressure) + 0.00030286 
        para_2 = -1341.9 * (1 / Pressure)**3 + 113.59 * (1 / Pressure)**2 + 19.183 * (1 / Pressure) + 0.00091879 
        para_1 = 4695.7 * (1 / Pressure)**3 - 404.31 * (1 / Pressure)**2 - 17.669 * (1 / Pressure) - 0.011962 
        para_0 = -3655.5 * (1 / Pressure)**3 + 334.47 * (1 / Pressure)**2 + 3.3007 * (1 / Pressure) + 0.029922 
    
    
        # set the reduced temperature with the critical temperature for ortho 
        T_red = Temperature / 33.22
        # set the thermal conductivity for the ortho hydrogen part 
        # this is a five part equation so set the constant for each part as a function of pressure 
        # these were determined from a curve fit of thermal conductivity divided by the density 
        ortho_4 = 0.194631 * (1 / Pressure)**3 + 0.010058 * (1 / Pressure)**2 - 0.0221291 * (1 / Pressure) - 0.0000364344 
        ortho_3 = -1.55936 * (1 / Pressure)**3 - 0.662601 * (1 / Pressure)**2 + 0.37491 * (1 / Pressure) + 0.000820967 
        ortho_2 = -27.3931 * (1 / Pressure)**3 + 10.6882 * (1 / Pressure)**2 + 4.38153 * (1 / Pressure) - 0.0034411 
        ortho_1 = 297.937 * (1 / Pressure)**3 - 60.1667 * (1 / Pressure)**2 + 7.63341 * (1 / Pressure) + 0.000249204 
        ortho_0 = -633.093 * (1 / Pressure)**3 + 101.419 * (1 / Pressure)**2 - 7.25464 * (1 / Pressure) + 0.0198012 
    
        # find the thermal conductivity of para and ortho hydrogen 
        para_therm_cond = (para_4 * T_red**4 + para_3 * T_red**3 + para_2 * T_red**2 + para_1 * T_red**1 + para_0) * rho
        ortho_therm_cond = (ortho_4 * T_red**4 + ortho_3 * T_red**3 + ortho_2 * T_red**2 + ortho_1 * T_red**1 + ortho_0) * rho
    
        # limit the range but it also depends on temperature because the graph resembles a check mark 
        if Temperature > 100:
            para_therm_cond = clip(para_therm_cond,0.05,0.3)
            ortho_therm_cond = clip(ortho_therm_cond,0.05,0.3)
        else:
            para_therm_cond = clip(para_therm_cond,0.01,0.12)
            ortho_therm_cond = clip(ortho_therm_cond,0.01,0.12)
            
        # add the two parts and return Btu/hr/ft/K 
        Thermal_Conductivity = para_therm_cond * ParaFraction + ortho_therm_cond * (1 - ParaFraction) 
        return Thermal_Conductivity

    def ThrmlCondvect(self,Temperature, Pressure, ParaFraction, Density=0):
        return np.vectorize(self.Thermal_Conductivity(Temperature, Pressure, ParaFraction, Density))
            

    def ThrmlCond(self,Temperature, Pressure, ParaFraction, Density=0): 
        # Thermal conductivity [Btu/hr/ft/K] 
        # Temperature [Kelvin], pressure [psia] 
        # fraction of para hydrogen [0-1] 
        # density [lbm/ft3]

        T = np.array(Temperature,dtype='float64')
        P = np.array(Temperature,dtype='float64')
        pH2 = np.array(Temperature,dtype='float64')
   
        if (Density):
            rho = Density
        else: 
            rho = self.Density_Leachman(Temperature, Pressure, 0) 
    
        # set the reduced temperature with the critical temperature for para 
        T_red = Temperature / 32.938 
        # set the thermal conductivity for the para hydrogen part 
        # this is a five part equation so set the constant for each part as a function of pressure 
        # these were determined from a curve fit of thermal conductivity divided by the density 
        para_4 = -4.7484 * (1 / Pressure)**3 + 0.39601 * (1 / Pressure)**2 + 0.060825 * (1 / Pressure) - 0.0000169 
        para_3 = 138.2 * (1 / Pressure)**3 - 11.591 * (1 / Pressure)**2 - 1.6886 * (1 / Pressure) + 0.00030286 
        para_2 = -1341.9 * (1 / Pressure)**3 + 113.59 * (1 / Pressure)**2 + 19.183 * (1 / Pressure) + 0.00091879 
        para_1 = 4695.7 * (1 / Pressure)**3 - 404.31 * (1 / Pressure)**2 - 17.669 * (1 / Pressure) - 0.011962 
        para_0 = -3655.5 * (1 / Pressure)**3 + 334.47 * (1 / Pressure)**2 + 3.3007 * (1 / Pressure) + 0.029922 
    
    
        # set the reduced temperature with the critical temperature for ortho 
        T_red = Temperature / 33.22
        # set the thermal conductivity for the ortho hydrogen part 
        # this is a five part equation so set the constant for each part as a function of pressure 
        # these were determined from a curve fit of thermal conductivity divided by the density 
        ortho_4 = 0.194631 * (1 / Pressure)**3 + 0.010058 * (1 / Pressure)**2 - 0.0221291 * (1 / Pressure) - 0.0000364344 
        ortho_3 = -1.55936 * (1 / Pressure)**3 - 0.662601 * (1 / Pressure)**2 + 0.37491 * (1 / Pressure) + 0.000820967 
        ortho_2 = -27.3931 * (1 / Pressure)**3 + 10.6882 * (1 / Pressure)**2 + 4.38153 * (1 / Pressure) - 0.0034411 
        ortho_1 = 297.937 * (1 / Pressure)**3 - 60.1667 * (1 / Pressure)**2 + 7.63341 * (1 / Pressure) + 0.000249204 
        ortho_0 = -633.093 * (1 / Pressure)**3 + 101.419 * (1 / Pressure)**2 - 7.25464 * (1 / Pressure) + 0.0198012 
    
        # find the thermal conductivity of para and ortho hydrogen 
        para_therm_cond = (para_4 * T_red**4 + para_3 * T_red**3 + para_2 * T_red**2 + para_1 * T_red**1 + para_0) * rho
        ortho_therm_cond = (ortho_4 * T_red**4 + ortho_3 * T_red**3 + ortho_2 * T_red**2 + ortho_1 * T_red**1 + ortho_0) * rho
    
        # limit the range but it also depends on temperature because the graph resembles a check mark 
        if Temperature > 100:
            para_therm_cond = clip(para_therm_cond,0.05,0.3)
            ortho_therm_cond = clip(ortho_therm_cond,0.05,0.3)
        else:
            para_therm_cond = clip(para_therm_cond,0.01,0.12)
            ortho_therm_cond = clip(ortho_therm_cond,0.01,0.12)
            
        # add the two parts and return Btu/hr/ft/K 
        Thermal_Conductivity = para_therm_cond * ParaFraction + ortho_therm_cond * (1 - ParaFraction) 
        return Thermal_Conductivity


    def Viscosity(self,Temperature, Pressure): 
        """ Viscosity [centiPoise] is not senstivie to ortho-para composition 
        Temperature [K], Pressure [psia] 
         it consists of a high temperature section, a knee section and a low temperature section and these temperature transitions depend upon the pressure"""
        # High and low temperature zone [K]
        temperature_H = 0.057 * Pressure + 22.6 
        temperature_L = 2.5 * math.log(Pressure) + 14.2
        # vH and vL are coefficients for high and low temperature zone respecive
        vH =np.array([ (-0.0000000030096 * Pressure**2 + 0.0000044929 * Pressure + 0.00048368), \
              (0.00000000002744 * Pressure**2 - 0.000000036277 * Pressure + 0.000038261),\
              (-5.3637E-14 * Pressure**2 + 0.00000000006793 * Pressure - 0.000000031886),\
              0])
        vL = np.array([(-0.014037 * math.log(Pressure) + 0.16369),\
              (0.0023049 * math.log(Pressure) - 0.019595),\
              (-0.00011352 * math.log(Pressure) + 0.00084578),\
              (0.0000018175 * math.log(Pressure) - 0.000012611)])
        Viscosity, ViscL, ViscH = 0, 0, 0  
        if Temperature >= temperature_H:
            for i in Range(0,3):
                Viscosity = Viscosity + vH[i] * pow(Temperature,i)
        elif Temperature <= temperature_L:
            for i in Range(0,3):
                Viscosity = Viscosity + vL[i] * pow(Temperature,i)     
        else:
            # use correlation for knee zone 
            # use a linear interpolation between the high an low values at their ends 
            # first the high temperature part
            for i in Range(0,3):
                ViscL = ViscL + vL[i] * pow(temperature_L,i)
                ViscH = ViscH + vH[i] * pow(temperature_H,i)     
            Viscosity = ViscH + (ViscL - ViscH) / (temperature_H - temperature_L) * (temperature_H - Temperature) 
        return Viscosity
    
    def ViscositySimple(self,T, P): 
        """ Viscosity [centiPoise] (not senstivie to ortho-para composition) 
        Temperature [K], Pressure [psia] 
        it consists of: a high temperature zone, a knee section and a low temperature section and these temperature transitions depend upon the pressure"""
        # limits for high and low temperature zone [K]
        T_H = 0.057 * P + 22.6 
        T_L = 2.5 * math.log(P) + 14.2
        # vH and vL are coefficients for high and low temperature zone respecive
        TV= np.array([1, T, T**2, T**3]) 
        vH =np.array([ (-0.0000000030096 * P**2 + 0.0000044929 * P + 0.00048368), \
              (0.00000000002744 * P**2 - 0.000000036277 * P + 0.000038261),\
              (-5.3637E-14 * P**2 + 0.00000000006793 * P - 0.000000031886),\
              0])
        vL = np.array([(-0.014037 * math.log(P) + 0.16369),\
              (0.0023049 * math.log(P) - 0.019595),\
              (-0.00011352 * math.log(P) + 0.00084578),\
              (0.0000018175 * math.log(P) - 0.000012611)])
        Viscosity, ViscL, ViscH = 0, 0, 0  
        if T >= T_L:
            ViscH = sum(vH * TV)
        if T <= T_H:
            ViscL = sum(vL * TV)
        # Combine by linear interpolation between the high an low values at their ends 
        Viscosity = ViscH + (ViscL - ViscH) / (T_H - T_L) * (T_H - T) 
        return Viscosity
    
    def PrNr(self):
        '''This function assumes that flash is calculated'''
        Result = 0
        mu = self.stream.viscosity()
        cp = self.stream.specificEnthalpy()
        k = self.ThermConduct()
        return mu*Cp/k





class catalyst:
    """
    Base class for catalysts. It contains parameters of the catalyst but also
    caculates exquilibrium o-p composition for a given temperature. 
   
    # info:
    name - name of the material
    type - catalyst

    Geometry:
    Morphology of the catalyst

    """
    def __init__(self, name="Generic Catalyst",type="catalyst"):
        primitive.__init__(self,name,type)
        mesh_size = 0.012 # 0.012-0.022 inch
        reactionConstants = 1
    

    def Equilibrium(self, T):
        '''Flyn[1] valid between 15K and 280K return
        values for limits for points outside'''
        Temp = np.asarray(T)
        Temp = clip(Temp,15.0,280.0)
        pH2eq = 0.1/(math.exp(-175/Temp)+0.1) - 7.06e-9*Temp**3 +3.42e-6*Temp**2 -6.2e-5*Temp - 0.00227
        return pH2eq

class stream(matter):
    """
    Class to handle process streams. It has methods for thermo properties from fluid and process parmeters for streams
     
    # Process Parameters:
    T_in - [K] Temperature inlet
    P_in - [Pa] Pressure inlet
    pH2 - [-] parahydrogen fraction inlet
    m - [kg/s] mass flow inlet
    P_out - [Pa] Pressure outlet
    h_in - [J/kg] Enthalpy inlet 
    fluid - refrence to the fluid
    
    """

    def __init__(self,name):
        """Initialise class."""
        matter.__init__(self,name,type="stream")
        # self.pH2 = None # fraction of parahydrogen        
        self.T_in = None # Inlet Temperature [K]      
        self.P_in = None # Pressure in [Pa]
        self.P_out = None # Pressure out [Pa]
        self.m = None # mass flow [kg/s]

        self.fluid = None # set fluid for channel
        self.h_in = None # mass Enthalpy [J/kg/K]
        self.dP_poly = None # Pressure drop correlation TODO:Rafal: Not sure if it is here or on the channel
        # compositions only for hydrogen because thermo props do not support
        self.compName = None # ['oH2','pH2'] # ,'N2','C1','C2','C3','iC4','nC4'] 
        self.comp = None 

    def setFluid(self, compNames,thermoDB, Thermo = 'None'):
        # TODO: here a switch between Thermo packages
        # self.fluid = fluidX("Fluid")

        dict = { # TODO: construct this list as a global at the class creation
            'h2fluid' : H2fluid,
            'h2ofluid': H2Ofluid,
            'sptts' : SPTTSFluid,
            'coolprop': CoolPropFluid
            }
        self.mat = dict.get(thermoDB.lower(), rError('Thermodynamic package {0} is not recognized for strem {1}'.format(thermoDB,self.name)))()
        # translation for compnents that are not in the database

        self.mat.setPropDB(compNames, Thermo)

    def setStream(self, compFract, m, Tin, Pin, Pout):
    # def setStream(self, compFract, m, Tin, Pin, Pout):
        cF = np.asarray(compFract,dtype=np.double) # composition vector
        self.comp = cF/sum(cF) # sum to 1
        self.T_in = np.double(Tin)  # Inlet Temperature [K]        
        self.m = np.double(m)  # mass flow [kg/s]
        self.P_in = np.double(Pin * 1000)  # Pressure in [Pa]
        self.P_out = np.double(Pout * 1000)  # Pressure out [Pa]

        # TODO: checks here that all is OK
        # you can request here to have Pout or DP specified.
        # if missing specify pressure drop correlation.
        self.h_in = self.mat.getProps(['TPFlash', self.T_in, self.P_in, self.comp],['EnthMass'])


    def check(self):
        """Function to check that class has been defined correctly."""
        # fluid type
        # self.fluid = fluidX()
        if not self.mat:
            raise rError('fluid not specified')
        if (not self.T_in) and (not self.h_in):
            raise rError('T_in (inlet Tempr) or h_in (inlet ?) not specified correctly')
        if not self.m:
            raise rError('m (mass flow) not specified')
        if not self.P_in:
            raise rError('P_in (inlet pressure) not specified')
        if not self.P_out:
            self.P_out = self.P_in
            print('PH_out (outlet pressure) not specified and set to equal PH_in ie. DP=0')
            
        if self.dP_poly != None and M.f_CorrelationH != 0:
            raise rError("dP_poly= {0} has been specified. This is not supported in combination with f_Correlation = {1}. To use pressure drop prescribed by polynominal, set M.f_CorrelationH = 0 \n Bailing Out!".format(self.dP_poly, M.f_CorrelationH))

class chnlfluid(stream,chnl):
    """
    (Derived) Class for fluid channels
    
    Proc. Param:
    chnlm  # [kg/s] flowrate through the individual channel (!) rather than stream
    - stream - refrence to the stream
    - isHot - [bool] is the channel hot ? (or cold = false)
    - Cat - refrence to catalyst data
    - hasCat - (bool) catalyst present?

    # Geometry (when finned):
    - finH  # [m] fin height
    - self.finTh # [m] fin thickness
    - finSp # [m] fin spacing  
    - finL # [m] length of the fins
    """

    def __init__(self,name):
        stream.__init__(self,name)
        self.type = "fluid"
        self.description = "Channel with a fluid flowing inside. The channel has no catalyst or solid inserts"
        self.hasCat = False

    def check(self):
        pass

    def connectStream(self,F,isHot):
        self.stream = F
        self.isHot = isHot

    def addCatalyst():
        self.Cat = catalyst()
        self.hasCat = True

    def setGeometry(self,chnlThickness):
        chnl.setGeometry(self,chnlThickness)

    def addFinGeometry(self, finHeight, finThickness, finSpacing, finLength):
        self.finH  = finHeight/1000     # mm -> m
        self.finTh = finThickness/1000  # mm -> m
        self.finSp  = finSpacing/1000   # mm -> m
        self.finL    = finLength/1000   # mm -> m
        self.chnlH = self.finH + self.finTh

    def addToHX(self,HXdata):
        ''' update channel parameters from HX parameters '''
        chnl.addToHX(self,HXdata)
        self.chnlm = self.stream.m / HXdata.nLyr # flowrate through the individual channel 
        # TODO: NOTE1: here I assume that a stream is only once in a layer which may not be true.
        # you may have ABCB arrangement where B (refrigerant) is twice 

        # NOTE2: when fins are considered we may want to think about it 
        # self.chnlH = self.finH + self.finTh
    

class htc:
    """Class template for heat transfer."""

    def __init__(self):
        """Initialise class."""
        self.type = []  # set what type of flow is supported.
        self.name = 'template'  # name by which correlation will be called
        self.description = ''  # provide a description of the correlation, include validity range
        self.reference = ''  # provide a reference for the correlation

    def getHTC(self,T, P): # define later
        """Calculate heat transfer."""
        # add formula to calculate  Nusselt number.
        HTC = None
        return [HTC, HTC_flag]

    def isValid(self, T, P):
        """Check validaity of using ."""
        # add code to calculate if function is within validaity range
        self.validity = 'not checked'

    def info(self):
        string = 'Correlation Name={0}; Types of flow supported: {1}; Validity check:'.format(self.name, self.type, self.validity)
        if long is False:
            return string
        else:
            [string + '\nDescription: {0} \nReference: {1}'.format(self.description, self.reference)]
            return string

class htc_pipes:
    """
    classical heat transfer correlation for pipes (dittus boelter)
    Nu = 0.023 * Re**0.8 * Pr**n
    """

    def __init__(self):
        """Initialise class."""
        self.name = 'DittusBoetler'  # name by which correlation will be called
        self.description = 'Dittus Boelter Correlation for Heat Transfer'  # provide a description of the correlation, include validity range
        self.reference = 'see ref XXX'  # provide a reference for the correlation
        self.ValidityRange = ''

    def getHTC(self, Hcc, Pcc, Twall, m, Area, CharLength):
        """Calculate heat transfer."""
        # Compute Bulk properties
        # Calculate bulk temperature
        
        state = ['PHFlash', Pcc, Hcc]
        prtis = ['Pr', 'massDens', 'Visc', 'Temp', 'ThermCond']
        props = getProps(state, prtis)
        Pr, rhocc, mucc, Tcc, condcc = props[0], props[1], props[2], props[3], props[4]

        # Calculate Reynolds number
        U = abs(m/(rhocc*Area))
        Re = rhocc*U*CharLength/mucc

        # Dittus Boelter Equation
        if Twall > Tcc:  # heating of fluid
            n = 0.3
        else:       # cooling of fluid
            n = 0.4
        Nu = 0.023 * Re**0.8 * Pr**n

        HTC = Nu * condcc / CharLength
        return [HTC, 'htc']

    def isValid(self, T, P):
        """Check validaity of using ."""
        # add code to calculate if function is within validaity range
        self.validity = 'not checked'

    def info(self, long=False):
        string = 'Correlation Name={0}; Types of flow supported: {1}; Validity check:'.format(self.name, self.type, self.validity)
        if long is False:
            return string
        else:
            [string + '\nDescription: {0} \nReference: {1}'.format(self.description, self.reference)]
            return string

class htc_Joshi_Webb1987: # TODO: finish based on the paper
    """Class template for heat transfer."""
    def __init__(self):
        """Initialise class."""     
        self.name = 'Joshi_Webb1987'  # name by which correlation will be called
        self.description = 'Corrlation for Off-set strip fins.'  # provide a description of the correlation
        self.reference = 'Joshi Webb 1987'  # provide a reference for the correlation
class channelFinned: # serrated channel
    """
    Info:
    name = "stream"
    Type = "fluid"

    Proc. Param:
    - chnlm  # [kg/s] flowrate through the individual channel 
    - stream - refrence to the stream
    - isHot - [bool] is the channel hot ? (or cold = false)
    - Cat - refrence to catalyst data
        
    Geometry:
    - finH  # [m] fin height
    - self.finTh # [m] fin thickness
    - finSp # [m] fin spacing  
    - finL # [m] length of the fins

    Correlations:
    - corr - object with correlations for calculation of pressure drop and HTC from Nusselt number
    """
    name = "stream"
    Type = "fluid"
    # CH = 0
    # chnlH = 0
    chnlL = 0 # TODO: can we remove since provided 
    chnlW = 0 # TODO: can we remove since provided 
        
    def __init__(self):
        # Geometry # reactor     # exchanger 
        finH  = []   # a 0.264 inch# 0.185 fin height [mm]
        finTh = []    # t 0.016    # 0.016 fin Thickness [mm]
        finSp  = []  # S 0.046     # 0.073 fin spacing [mm]
        finL    = []  # l 0.25      # 0.5 serration length [mm]
        # I think not here
        plateThickness = [] 
        chnlL = 1 # length of the channel
        cat = []

    #hasCatalyst
    def check(self):
        """Check what needs to be defined to start calculations"""

        # if not self.finH: # and (self.finThickness) and (self.finSpacing, finLenth
        #    raise rError('Channel geometry not defined. Please provide\
        #    finHeight, finThickness, finSpacing, finLenth')
        print("Channel check: pass")
        
    def setGeometry(self, finHeight, finThickness, finSpacing, finLength):
        self.finH  = finHeight/1000     # mm -> m
        self.finTh = finThickness/1000  # mm -> m
        self.finSp  = finSpacing/1000   # mm -> m
        self.finL    = finLength/1000   # mm -> m
        self.chnlH = self.finH + self.finTh

    def addToHX(self,HXdata):
        ''' update parameters from HX parameters '''
        self.chnlH = self.finH + self.finTh
        self.chnlW = HXdata.GeomW
        self.chnlL = HXdata.GeomH
        self.Area_NS = self.chnlH * self.chnlW
        self.Area_WE = self.chnlL * self.chnlW / HXdata.lyr.nCells        
        self.chnlm = self.stream.m / HXdata.nLyr # flowrate through the individual channel 
        # TODO: here I assume that a stream is only once in a layer which may not be true.
        # you may have ABCB arrangement where B (refrigerant) is twice 

    def setCatalyst(self, C):
        self.Cat = C

    def connectStream(self,F,isHot):
        self.stream = F
        self.isHot = isHot

    def setHTC(self, Name=None):
        # get directory and pick the correct one
        self.corrHTC = None
        pass

    def setDP(self, Name=None):
        # get directory and pick the correct one
        self.corrDP = None
        pass

class layer:
    CH = [] # channels
    LH = 0 # [mm] layer height
    nCells = None # [-] number of discretization cells
    # xCell = [] # position of the middle of the cells along the lenght of the channel
    # xFace = [] # position of the face (cell wall) along length of the channel 
    # Predefined fields to make live graphs during solver execution fast ie. update fields rather than redraw on the canvas. 
    # TODO: this is really ugly but required by the solver
    # Face values:
    # T, P, pH2, H = [], [], [] # tempearture, pressure, paraHydrogen, enthalpy content field on the faces of cells
    # mCh, sCh = [], [] # list of metal and streams for channel iterations
    # CH - reference to Channels
            
    def __init__(self,nCells):
        CH = [] # list of channels, metals
        self.nCells = nCells
        self.methodSS = 'hybr' # default
        self.methodDyn = 'RK23' # default
        self.tSpan = (0,2000) # by default run for 2000s
    def check(self):
        """Check what needs to be defined to start calculations"""
        if len(self.CH) < 4:
            raise rError('Specified {0} channels. At least 2 fluid \
            and 2 metal channels required'.format(len(self.CH)))
        
    def addChannel(self,CH,M):
        self.CH.append(CH)
        self.CH.append(M)
        print(" Layers added: {0}".format(len(self.CH)))
        self.lyrH = CH.chnlH + M.chnlH 

    def addGrid(self,Length,nCells):
        # make computational grid
        # equaly spaced
        xFace = np.zeros(nCells+1)
        dist = Length/nCells
        for i in range(nCells + 1):
            xFace[i] = i * dist
        xCell = np.zeros(nCells)
        for i in range(nCells):
            xCell[i] = (xFace[i]+xFace[i+1])/2
        self.xFace = xFace
        self.xCell = xCell 

        # for now take, this part of code is in a wrong module. this should be calculated once only and refered 
        # TODO: change it for adaptive mesh. most likely in the future needed for tracking phase transition
        # DXNS = np.ones(nCells) # distances - the same for all channels
        # DXNS *= self.CH[0].chnlL / nCells /2.0 # channel length by number of cells and by 2.
        # DXWE = np.zeros(nChnls)

        # for i in range(len(self.CH)): # distances to walls characteristic to channels
        #    DXWE[i] = self.CH[i].chnlH / 2.0       
        self.nCells = nCells

    def addToHX(self,HXdata):
        for L in range(len(self.CH)):
            self.CH[L].addToHX(HXdata)
        self.addGrid(HXdata.GeomH,self.nCells)

    def setSS(self, method):
        if method=='RK45' or method=='RK23' or method=='DOP853':
            # method=='Radau' or method=='BDF' or method=='LSODA'
            self.methodSS = method
        else:
            rError("{0} is an unknown solver for time dependent simulation".format(method))

    def setDyn(self, method, tSpan):
        '''set dynamic solver
        method - solver
        tSpan - (Tstart Tstop) integration time span
        '''
        if method=="hybr" or method=='lm' or method=='lmN':
            # method=='Radau' or method=='BDF' or method=='LSODA'
            self.methodDyn = method
            self.tSpan = tSpan
        else:
            rError("{0} is an unknown solver for time dependent simulation".format(method))



    def print(self):
        print("Number of layers defined:",len(self.CH))
        for ind in range(len(self.CH)):
            print("Channel: {0}-> {1}".format(ind,self.CH[ind].name))
    
    def setBC(self,T,i):
        if (self.CH[i].isHot): # for hot stream
            T[0] = self.CH[i].stream.T_in
        else: # for cold stream
            T[-1] = self.CH[i].stream.T_in
        return T

    def drawTempPlots(self,info):
        plt.clf()
        markers = ['o','o','x','+','','o','x','+','']
        
        for i in range(len(self.T)):
            plt.plot(self.xFace,self.T[i], marker=markers[i], color=mycol[i],  label=str(i))
   
        plt.title('Water water system {0}'.format(info))
        plt.ylabel('Temperature')
        plt.xlabel('Distance [m]')
        plt.legend()

        plt.draw()
        plt.pause(0.00001)


    def initVariables(self):
        """ make initialization of variables 
        - in the future possibly read in from a file
        """
      
        apprT = 2.0        
        nCells = self.nCells
        nFaces = nCells+1
        nChnls = len(self.CH)
        self.comp = [None] * nChnls

        FM = np.ones([nChnls,nCells+1]) # face matrix for initialization
        h = FM.copy() # enthalpy
        # setup iteration lists for materials and streams
        if (not hasattr(self, 'mCh')):      
            """ write tables with indexes metal or stream channels
            if metal no temperature set
            if stream hot: T set at 0. do P interpolation
            if stream is not hot at nCells+1 do P interpolation reverse 
            """
            # TODO: change channel structure to have name and type type should be metal and stream
            self.mCh, self.sCh = [],[] # metal and stream channel list for iteration        
            for i in range(nChnls):
                if self.CH[i].type=="metal":
                    self.mCh.append(i)
                elif (self.CH[i].type=="fluid"):
                    self.sCh.append(i)
                else:
                    print("stream {0} of type {1}".format(i,self.CH[i].name))
                    raise rError('Unknown Channel type ')
           
           
           # create a field of uniform composition for a channel
        for i in self.sCh:
            cmp = self.CH[i].stream.comp
            tmp = np.zeros((nFaces,cmp.size))
            for j in range(nFaces):
                for k in range(cmp.size):
                    tmp[j,k] = cmp[k]
            self.comp[i] = tmp.copy() 

        # note: Field at Cell Center is average between oposite cell faces since differential size is assumed
        if (not hasattr(self, 'T')): # is it defined
            # TODO: for now init normal pH2, change once reaction is in place 
            self.T, self.P,  self.H = FM.copy(),  FM.copy() * 0.25, FM.copy()
            # self.pH2 = FM.copy()
        # initialization of tempeartues for: 
        # hot stream: find hotest out of cold streams and add appr temperature of 5 K
        # cold stream: find coldest out of hot streams and take away appr temperature of 5 K
        # then interpolate between inlet and outlet
        THot,TCold = 0, 1000, 
        for i in self.sCh:
            if (self.CH[i].isHot): # for hot stream
                TCold = min(TCold,self.CH[i].stream.T_in) 
            else:
                THot = max(THot,self.CH[i].stream.T_in)
        THot += apprT 
        TCold -= apprT
        
        DQexchInit = np.ones(nChnls)*1e10 # note: large values init
        # Calculate enthalpy change for inlet and outlet
        # TODO: once thermodynamics in place add init for pH2 (equilbrium at every point) and check whether it makes any difference 
        # init T, P, calculate enthalpy change per pass.
        for i in self.sCh:
            if (self.CH[i].isHot): # for hot stream 
                self.T[i] = np.interp(self.xFace,[self.xFace[0],self.xFace[-1]],[self.CH[i].stream.T_in, THot])
                # T[i,0] = self.CH[i].stream.T_in # boundary condition
                self.P[i] = np.interp(self.xFace,[self.xFace[0],self.xFace[-1]],[self.CH[i].stream.P_in, self.CH[i].stream.P_out])
                # TODO: no pressure drop correlation here, simple DP implementd.
            else: # for cold stream reverse
                # self.T[i] = np.linspace(TCold, self.CH[i].stream.T_in, nCells)
                self.T[i] = np.interp(self.xFace,[self.xFace[0],self.xFace[-1]],[TCold, self.CH[i].stream.T_in])
                # T[i,-1] = self.CH[i].stream.T_in # boundary condition
                self.P[i] = np.interp(self.xFace,[self.xFace[0],self.xFace[-1]],[self.CH[i].stream.P_out, self.CH[i].stream.P_in])
                # TODO: no pressure drop correlation here. simple DP implementd
            # enthalpy change with initial initialization in channel "i"

            Ctmp = self.comp[i]
            state = [ 'TPFlash', [self.T[i,0], self.T[i,-1]], 
                    [self.P[i,0], self.P[i,-1]] , [Ctmp[0], Ctmp[-1] ] ]
            getProps = self.CH[i].stream.mat.getProps(state,['EnthMass'])
            self.H[i,0], self.H[i,-1] = getProps[0,0], getProps[0,1]     

            DQexchInit[i] = (self.H[i,-1]-self.H[i,0]) * self.CH[i].chnlm

        self.drawTempPlots("init")

        QExch = min(abs(DQexchInit)) # Heat Flows
        DQ = np.zeros(nFaces)
        DQ = np.interp(self.xFace,[self.xFace[0],self.xFace[-1]],[0, QExch])

        for i in self.sCh:
            if (self.CH[i].isHot): # for hot stream
                self.H[i] = self.H[i,0] - DQ / self.CH[i].chnlm
            else:
                self.H[i] = self.H[i,-1] + (QExch-DQ) / self.CH[i].chnlm

            self.T[i] = self.CH[i].stream.mat.getProps(['PHFlash', self.P[i], self.H[i], self.comp[i]],['Temp'])
            x= 1
            # old definition
            # self.T[i] = self.CH[i].stream.mat.getTfromHP(self.T[i],self.P[i],self.H[i],self.pH2[i,0])
        # TODO: now I should recalculate temperatures of all streams.

        # initialize global variables for metals
        for i in self.mCh:
            # fill in init for the metal part
            # TODO: average is incorrect but for now has to do.
            j = i + 1
            if j >= nChnls:
                j = 0 
            self.T[i] = (self.T[j]+self.T[i-1])/2.0
            self.H[i] = self.T[i] * self.CH[i].Cp(self.T[i])

        self.drawTempPlots("start")

        self.count =0  
        s =1 


    def dHdtSS(self, Hcc, time=None):
        # unified function for steady state and dynamic 
        # in SS - dHdt should be 0 and the value is treated as error that needs to be minimized.

        # for convenience create N,S,W,E subscripts North, South etc. on the faces of the cell (and make values average) 
        # useful declarations:
        DynSim = (time is not None)

        nCells = self.nCells
        nChnls = len(self.CH)
        FM0, FM1= np.zeros([nChnls,nCells+1]), np.ones([nChnls,nCells+1])
        CM0, CM1 = np.zeros([nChnls,nCells]), np.ones([nChnls,nCells])
        # create distances mesh
        DX = self.xFace[1:]-self.xFace[:-1]

        # declare values at the cell center 
        rho, QR, dHdt = CM1.copy(), CM0.copy(), CM0.copy() # average density, heat of reaction, Pressure
        TE, TW = CM0.copy(), CM0.copy()
        alpha, k, h = CM0.copy(), CM0.copy(), CM0.copy()

        Hcc = np.reshape(Hcc,(-1,nCells))


        compcc = [None] * nChnls

        self.H = intEx(self.xFace, self.xCell, Hcc)
        # restore Boundary Conditions
        for i in self.sCh:
            ch = self.CH[i]
            if ch.isHot: # Restore boundary conditions, TODO(hold): later check whether need to restore other params: P, T, H, pH2
                self.H[i,0] = ch.stream.h_in
            else:                
                self.H[i,-1] = ch.stream.h_in
            # self.T[i] = ch.stream.mat.getTfromHP(self.T[i],self.P[i],self.H[i],self.pH2[i])
            self.T[i] = ch.stream.mat.getProps(['PHFlash', self.P[i], self.H[i], self.comp[i]],['Temp'])
        self.drawTempPlots("mid")


        for i in self.mCh:
            self.T[i] = self.H[i]/self.CH[i].Cp(self.T[i])
        
        self.drawTempPlots("mid")
        
        TN, TS = self.T[:,:-1], self.T[:,1:] 
        hN, hS = self.H[:,:-1], self.H[:,1:]

        # create cell center values
        Tcc = (TN + TS) * 0.5
        Pcc = (self.P[:,:-1]+self.P[:,1:])*0.5
        # pH2cc = (self.P[:,:-1]+self.P[:,1:])*0.5

        for i in self.sCh:
            # need to iterate because metal spaces are empty
            cmp = self.comp[i]
            compcc[i] = (cmp[:-1,:]+cmp[1:,:])*0.5

        # calculate properties required for derivative
        for i in self.sCh:
            state = ['TPFlash', Tcc[i], Pcc[i], compcc[i]]
            props = ['EnthMass', 'massDens','ThermCond']
            getProps = self.CH[i].stream.mat.getProps(state,props)
            h[i], rho[i], k[i]= getProps[0], getProps[1], getProps[2]
        for i in self.mCh:
            k[i] = self.CH[i].k(Tcc[i])
            rho[i] = self.CH[i].rho(Tcc[i]) # later make it dependent on temperature and other

        ######################
        # calculate HTC (heat transfer coefficient) fluid -> metal
        # TODO: introduce Nuselt number
        # Calculate midpoint property values - (conductivity)      
        # Nu = hL/k k - conductivity, h convective heat transfer
        # htc now calculated with correlations.
        # htcH = NuH * kH / G.Lc_H
        for i in self.sCh:
            ch = self.CH[i]
            #DP[i] = ch.corr.dp(ch)
            # alpha = ch.corr.alpha(ch)
        alpha = CM1.copy()*25 # for water W/m2/K
        ######################
        
        # wall temperature WE        
        for i in self.sCh:
            TW[i] = Tcc[i] + k[i-1]/(k[i-1]+alpha[i]*self.CH[i-1].chnlH/2) * (Tcc[i-1]-Tcc[i])
            TE[i-1] = TW[i]
            TE[i] = Tcc[i] + k[i+1]/(k[i+1]+alpha[i]*self.CH[i+1].chnlH/2) * (Tcc[i+1]-Tcc[i])
            TW[i+1] = TE[i]
        
        if DynSim:
            # TODO: write here something not to display every step
            self.drawTempPlots(" at {0} min".format(time))
            pass 
        else:
            self.drawTempPlots(" {0}".format(" steady state"))
            pass

        # formulate derivative:
        for i in self.sCh: # for fluids
            ch = self.CH[i]
            chnlm = ch.chnlm * (1-ch.isHot*2) # mass flow with sign (cold flow in oposite direction)
            dHdt[i] = (chnlm * (hS[i]-hN[i]) + QR[i] \
                + alpha[i] * (TW[i]-Tcc[i]) * ch.Area_WE + alpha[i] * (TE[i] - Tcc[i]) * ch.Area_WE ) \
                / (rho[i] * DX * ch.Area_NS)
        for i in self.mCh: # for metal
            ch = self.CH[i]
            dHdt[i] = (k[i] * (TN[i]-Tcc[i])/DX * 2 * ch.Area_NS + k[i] * (TS[i]-Tcc[i])/DX * 2 * ch.Area_NS \
                +  k[i] * (TW[i]-Tcc[i])/ch.chnlH * 2  * ch.Area_WE +   k[i] * (TE[i]-Tcc[i])/ch.chnlH * 2 * ch.Area_WE ) \
                / (rho[i] * DX * ch.Area_NS)

        self.count +=1         
        # name = "Hcc_"+ str(self.count ) + ".txt"
        # np.savetxt(name, np.c_[Hcc[0]], header='Hcc file', fmt='%d',delimiter='\t')

        dHdt = dHdt.flatten('C')
        return dHdt

    def dHdtDyn(self,time,Hcc):
        # ode solver requires different argument order which is different than for SS solvers.
        return self.dHdtSS(Hcc, time)
        
    def dHdtSSmin(self,time,Hcc):
        '''sum square for certain solvers'''
        error = self.GenEquNew(Hcc, True, time)
        return sum(error**2)

    def solve(self,runDyn): # SS or Dyn
        '''
        main solver by 
        '''
        # def dHdtDyn(self,time,Hcc):
        #     ''' ode solver requires different argument order which is different than for SS solvers '''
        #     return self.dHdtSS(Hcc, time)

        # def dHdtSSmin(self,time,Hcc):
        #     '''sum square for certain solvers'''
        #     error = self.GenEquNew(Hcc, True, time)
        #     return sum(error**2)

        sol, status, mesg, success = None, None, None, None
        self.initVariables()
        Hcc = (self.H[:,1:] + self.H[:,:-1]) * 0.5
        Hcc = Hcc.flatten('C')

        if runDyn: # run in Dynamics
            tSpace = np.linspace(self.tSpan[0], self.tSpan[1], 50)
            method = self.methodDyn
            # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
            if method=='RK45' or method=='RK23' or method=='DOP853' :
                # Explicit Runge-Kutta methods ('RK23', 'RK45’, 'DOP853’) should be used for non-stiff problems and implicit methods 
                # 'DOP853’ is recommended for solving with high precision (low values of rtol and atol).
                options={'atol': 1.e-6, 'atol': 1.e-12} # Default values are 1e-3 for rtol and 1e-6 for atol.
                sol = odeInt.solve_ivp(self.dHdtDyn, self.tSpan, Hcc, method=method, t_eval=tSpace, dense_output=False, events=None, vectorized=True, options={'atol': 1.e-4, 'atol': 1.e-8})
            elif method=='Radau' or method=='BDF' or method=='LSODA' :
                # 'Radau’, 'BDF’ for stiff problems
                # first try to run 'RK45’. If fails 'Radau’ or 'BDF’. 'LSODA’ can also be a good universal choice
                # 'Radau', 'BDF', 'LSODA' require jacobian 
                raise rError("Method {0} requires Jackobian set. Jacobian construction is not implemented yet.", method)
            else:
                raise rError("Solver note set or unknown or typo")
            self.H = sol.y
        else: # Run in steady state
            method = self.methodSS
            # method = 'lmN'
            if method == 'hybr':
                sol = optimize.root(self.dHdtSS, Hcc, method=method, options={'xtol': 1.e-12})
            elif method == 'lm':
                sol = optimize.root(self.dHdtSS, Hcc, method='lm', options={'eps': 1.e-3, 'xtol': 1.e-12, 'ftol': 1e-12})
            elif method == 'lmN': # 'Newton-CG': # what is the difference with the below?
                # sol = optimize.root(self.dHdtSS, Hcc, method='lm', options={'eps': 1.e-3, 'xtol': 1.e-12})
                sol = optimize.root(self.dHdtSS, Hcc, method='lm', options={'eps': 1.e-6, 'xtol': 1.e-14})
            elif method == 'df-sane':
                sol = optimize.root(self.dHdtSS, Hcc, method=method, options={'ftol': 1.e-12})
            elif method == 'L-BFGS-B' or method == 'TNC':
                Hlim = [min(Hcc), max(Hcc)]
                bounds = []
                for i in range(len(Hcc)):
                    bounds.append((Hlim[0], Hlim[1]))
                sol = optimize.minimize(dHdtSSmin,Hcc,args=(runDyn),method=method,bounds=bounds,options={'disp':True})
            elif method =='fsolve':
                sol, infodict, status, mesg = optimize.fsolve(self.dHdtSS, Hcc, full_output=1)
            else:
                raise rError("solver note set.")
            if status == None: # all methods except for fsolve
                self.H = sol.x
                success = sol.success
                status = sol.status
                mesg = sol.message            
                # params for for verbose
            else: 
                self.H = sol

        if runDyn:
            self.drawTempPlots(" at {0} min".format(round((self.tSpan[1])/60)))
        else:
            self.drawTempPlots(" - steady state (Final)")
        runDyn = None
 
class HXcore:
    # L is None # Layers class
    def __init__(self):
        # TODO define later
        self.lyr = [] # Layers class
    
    def check(self):
        """Check what needs to be defined to start calculations"""
        if not self.lyr:
            raise rError('No layers specified')
        if self.lyr.count() < 4:
            raise rError('Specified {0} layers. At least 2 fluid \
            and 2 metal required'.format(self.lyr.count()))
        # TODO: loop thorugh all layers to check whether properly defined.
        # TODO: check whether solver defined.
        
    def addlayers(self,L,nLayers,GeomH,GeomW):
        # Add layers only if they were not added before
        if not self.lyr:
            self.lyr = L
            self.nLyr = nLayers       # [-] number of layers in the core
            self.GeomH = GeomH/1000 # [mm]->m Height of the exchanger = length of channel
            self.GeomW = GeomW/1000 # [mm]->m Width of the exchagner
            self.GeomL = self.lyr.lyrH * self.nLyr
            self.lyr.addToHX(self)
            print("Heat exchanger dimensions WxLxH: {0:3.1f}m x {1:3.1f}m x {2:3.1f}m".format(self.GeomW,self.GeomL,self.GeomH))
        else:
            rError('HXcore has layers added already')

    def solve(self,runDyn=False):
        # check whether all is set up good
        self.lyr.solve(runDyn)
   

    def plot(self):
        """plot solution temperatures"""
        # TODO: plot temperature profiles mark which are cold and hot by checking isHot flag.  
        pass



###################################################################
###################################################################
###################################################################
# test Themro

# F1 = stream("HotH2")
# compNames = ["H2O", "1P"] # use Shell Thermo names, I translate comp names for every fluid.
# compFracs = [1.0, 0.0]
# F1.setFluid(compNames,'SPTTS', 'CPA/SMIRK')
# F1.setStream(compFracs,0.03,362.0, 253.0,250.0) # setStream(self,compFracs, m, Tin, Pin, Pout):
# F1.check()

# F2 = stream("HotH2")
# compNames = ["H2O", "1P"] # use Shell Thermo names, I translate comp names for every fluid.
# compFracs = [1.0, 0.0]
# F2.setFluid(compNames,'h2Ofluid', 'CPA/SMIRK')
# F2.setStream(compFracs,0.03,362.0, 253.0,250.0) # setStream(self,compFracs, m, Tin, Pin, Pout):
# F2.check()

# T = [250.0, 290.0, 340.0]
# P = [1.0e5, 1.0e5, 1.0e5]
# C = [compFracs, compFracs, compFracs]

# state = ['TPFlash', T, P, C]
# prtis = ['EnthMass']
# K = F1.mat.getProps(state,prtis)
# L = F2.mat.getProps(state,prtis)
# print(K)
# print(L)


### setup: metals
# aluminium 
# https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
# https://www.enggcyclopedia.com/2011/09/absolute-roughness/

Al = plate("Alumninium")
Al.setConstProps(2710.0,237.0,887.0) # density, thermal conductivity, Heat Capacity
Thickness = 1.6 # [mm] = 0.064 inches thick
Al.setPlate(Thickness,0.0015) # channel thickness, roughness
SS = plate("Stainless steel type 304")
SS.setConstProps(8000.0,14.4,468)
SS.setPlate(Thickness,0.015)

water = True

if water:
    # setup: Streams
    Feed = stream("HotH2")
    compNames = ["H2O", "1P"] # use Shell Thermo names, I translate comp names for every fluid.
    compFracs = [1.0, 0.0]

    # compNames, compFracs= ['oH2','pH2'], [0.75,0.25]
    # compNames = ['H2O','N2','1P','2P','3P','4P1','4P']
    # compFracs = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # TODO: remove pH2 from setStream

    # Feed.setFluid(compNames,'CoolProp', True)
    # Feed.setFluid(compNames,'SPTTS', 'CPA/SMIRK')
    Feed.setFluid(compNames,'h2Ofluid')
    Feed.setStream(compFracs,0.03,362.0, 253.0,250.0) # setStream(self,compFracs, m, Tin, Pin, Pout):
    Feed.check()

    compNames = ['H2O','1P']
    compFracs = [1.0, 0.0]

    Refr = stream("Refr")
    # Refr.setFluid(compNames,'CoolProp', True)
    # Refr.setFluid(compNames,'SPTTS', 'CPA/SMIRK')
    Refr.setFluid(compNames,'h2Ofluid')
    Refr.setStream(compFracs,0.05,295.0, 253.0,250.0)
    Refr.check()

    # setup: channel
    finHeight  = 6.7    # [mm] a 0.264 inch# 0.185
    finThickness = 0.4  # [mm] t 0.016    # 0.016
    finSpacing  = 1.16  # [mm] S 0.046     # 0.073
    finLenth    = 6.35  # [mm] l 0.25      # 0.5 serration length

    CH1 = chnlfluid("H2 stream")
    # CH1.setGeometry(finHeight, finThickness, finSpacing, finLenth)
    CH1.connectStream(Feed,True)
    CH1.setGeometry(10) # mm
    # CH1.check()

    CH3 = chnlfluid("Refrigerant")
    # CH3.setGeometry(finHeight, finThickness, finSpacing, finLenth)
    CH3.connectStream(Refr,False)
    CH3.setGeometry(10) # mm
    # CH3.check()

    # setup: Streams
    Feed3 = stream("HotH2")
    # Feed3.setFluid(compNames,'CoolProp', True)
    # Feed3.setFluid(compNames,'SPTTS', 'SMIRK')
    Feed3.setFluid(compNames,'h2Ofluid')

    Feed3.setStream(compFracs,0.01,342.0, 253.0,250.0) # setStream(self,compFracs, m, Tin, Pin, Pout):

    # second hydrogen stream
    CH5 = chnlfluid("H2 stream small")
    # CH1.setGeometry(finHeight, finThickness, finSpacing, finLenth)
    CH5.connectStream(Feed3,True)
    CH5.setGeometry(10) # mm
    # setup: Layers
    L = layer(50) # nr of cells 
    L.addChannel(CH1,Al)
    L.addChannel(CH3,Al)
    # L.addChannel(CH5,Al)

    L.print()

    # setup: Heat Exchangers
    HEX = HXcore()
    # HEX.addlayers(L,500,7600,2000)
    HEX.addlayers(L,50,760,200)
    # HEX.addlayers(L,500,1000,200) # L,nLayers,GeomH,GeomW):

    # HEX.addSolver(S)
    HEX.solve(True) # run in dynamic? default: False
    # HEX.plotTemp()
    wait =1 
else:
    
    # setup: Streams
    # SPTTS has no oHydrogen and pHydrogen defined at least this version
    # CoolProp has oHydrogen and pHydrogen defined but is missing binary interaction parameters
    # CoolProp prevents me to even enter default rules for them interaction parameters defined because they have the same CAS in the system (!) 
    
    Feed = stream("HotH2")
    
    # compNames, compFracs= ['oH2','pH2'], [0.75,0.25]
    compNames, compFracs= ['H2','He'], [0.75,0.25]
    # Feed.setFluid(compNames,'CoolProp', True)
    Feed.setFluid(compNames,'SPTTS', 'CPA')
    Feed.setStream(compFracs,0.03,90.0, 280.0,250.0) # setStream(self,compFracs, m, Tin, Pin, Pout):
    Feed.check()

    compNames, compFracs= ['oH2','pH2'], [0.75,0.25]
    Refr = stream("Refr")
    # Refr.setFluid(compNames,'CoolProp', True)
    # Refr.setFluid(compNames,'SPTTS', 'CPA')
    Refr.setFluid(compNames,'CoolProp', True)
    Refr.setStream(compFracs,0.05,50.0, 280.0,250.0)
    Refr.check()

    CH1 = chnlfluid("H2 stream")
    CH1.connectStream(Feed,True)
    CH1.setGeometry(10) # mm

    # setup: channel
    finHeight  = 6.7    # [mm] a 0.264 inch# 0.185
    finThickness = 0.4  # [mm] t 0.016    # 0.016
    finSpacing  = 1.16  # [mm] S 0.046     # 0.073
    finLenth    = 6.35  # [mm] l 0.25      # 0.5 serration length
    # CH1.setGeometry(finHeight, finThickness, finSpacing, finLenth)
    # CH1.check()

    CH3 = chnlfluid("Refrigerant")
    CH3.connectStream(Refr,False)
    CH3.setGeometry(10) # mm

    # setup: Layers
    L = layer(50) # nr of cells 
    L.addChannel(CH1,Al)
    L.addChannel(CH3,Al)

    L.setSS('RK23') # 'RK45' 'RK23' 'DOP853'
    L.setDyn('lm',(0,1000)) # 0-1000s

    L.print()

    # setup: Heat Exchangers
    HEX = HXcore() # create object 
    HEX.addlayers(L,500,7600,2000) # L,nLayers,GeomH,GeomW):

    # HEX.addSolver(S)
    HEX.solve(False) # run in dynamic? default: False