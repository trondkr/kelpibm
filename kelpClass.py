# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import os
import numpy as np
from datetime import datetime, timedelta

from opendrift.models.opendrift3D import OpenDrift3DSimulation, Lagrangian3DArray
from opendrift.elements import LagrangianArray

# Defining the oil element properties
class PelagicPlankton(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs
    """

    variables = LagrangianArray.add_variables([
        ('area', {'dtype': np.float32,
                         'units': 'm2',
                         'default': 0.}),
        ('length', {'dtype': np.float32,
                         'units': 'm',
                         'default': 0.}),
        ('diameter', {'dtype': np.float32,
                         'units': 'm',
                         'default': 0.}),
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0.}),
        ('plantpart', {'dtype': np.int32,
                     'units': '',
                     'default': 0}),             
        ('weight', {'dtype': np.float32,
                         'units': 'kg',
                         'default': 0.}),
        ('volume', {'dtype': np.float32,
                         'units': 'm3',
                         'default': 0.}),
        ('density', {'dtype': np.float32,
                         'units': 'kg/m3',
                         'default': 0.})])


  # def updateEggDevelopment(self):
  #      # Update percentage of egg stage completed
  #      amb_duration = np.exp(3.65 - 0.145*self.environment.sea_water_temperature) #Total egg development time (days) according to ambient temperature (Ellertsen et al. 1988)
  #      days_in_timestep = self.time_step.total_seconds()/(60*60*24)  #The fraction of a day completed in one time step
  #      amb_fraction = days_in_timestep/(amb_duration) #Fraction of development time completed during present time step 
  #      self.elements.stage_fraction += amb_fraction #Add fraction completed during present timestep to cumulative fraction completed
  #      self.elements.hatched[self.elements.stage_fraction>=1] = 1 #Eggs with total development time completed are hatched (1)


class PelagicPlanktonDrift(OpenDrift3DSimulation,PelagicPlankton):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

        Under construction.
    """

    ElementType = PelagicPlankton

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_ice_area_fraction', 'land_binary_mask',
                          'sea_floor_depth_below_sea_level',
                          'ocean_vertical_diffusivity',
                          'sea_water_temperature',
                          'sea_water_salinity',
                          'turbulent_kinetic_energy',
                          'turbulent_generic_length_scale',
                          'upward_sea_water_velocity',
                          'x_wind',
                          'y_wind'
                          ]

    # Vertical profiles of the following parameters will be available in
    # dictionary self.environment.vertical_profiles
    # E.g. self.environment_profiles['x_sea_water_velocity']
    # will be an array of size [vertical_levels, num_elements]
    # The vertical levels are available as
    # self.environment_profiles['z'] or
    # self.environment_profiles['sigma'] (not yet implemented)
    required_profiles = ['sea_water_temperature',
                         'sea_water_salinity',
                         'ocean_vertical_diffusivity']
    required_profiles_z_range = [-150, 0]  # The depth range (in m) which
                                          # profiles shall cover

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_floor_depth_below_sea_level': 100,
                       'ocean_vertical_diffusivity': 0.02,  # m2s-1
                       'sea_water_temperature': 10.,
                       'sea_water_salinity': 34.,
                       'surface_downward_x_stress': 0,
                       'surface_downward_y_stress': 0,
                       'turbulent_kinetic_energy': 0,
                       'turbulent_generic_length_scale': 0,
                       'upward_sea_water_velocity': 0
                       }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta'}


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(PelagicPlanktonDrift, self).__init__(*args, **kwargs)

    def update_terminal_velocity(self, Tprofiles=None, Sprofiles=None, z_index=None):
        """Calculate terminal velocity for kelp

        Created by Trond Kristiansen, 10.04.2017 

        """
        g = 9.81  # ms-2

        # Pelagic kelp properties that determine buoyancy
        kelparea = self.elements.area 
        kelpdiameter=self.elements.diameter  # assuming spoherical downward facing area
        kelpvolume = self.elements.volume 
        kelplength = self.elements.length  
        kelpmass = self.elements.weight 
        kelpdensity = self.elements.density 

        # prepare interpolation of temp, salt
        if not (Tprofiles==None and Sprofiles==None):
            if z_index==None:
                z_i = range(Tprofiles.shape[0]) # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],z_i,bounds_error=False) # evtl. move out of loop
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Tprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into this function, if not, use reader by calling self.environment
        if Tprofiles==None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * weight_upper + Tprofiles[lower, range(Tprofiles.shape[1])] * (1-weight_upper) 
        if Sprofiles==None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * weight_upper + Sprofiles[lower, range(Sprofiles.shape[1])] * (1-weight_upper) 

        # The density difference bettwen a pelagic egg and the ambient water
        # is regulated by their salinity difference through the
        # equation of state for sea water.

        DENSw = self.sea_water_density(T=T0, S=S0)
        dr = kelpdensity - DENSw  # density difference
       
        # water dynamic viscosity
        # Cylinder falling through water of density. Forces accounted fro is weight, drag, and buoyancy
        # https://www.physicsforums.com/threads/terminal-velocity-of-a-cylinder-freefalling-through-a-fluid.871446/
        my_w = 0.001*(1.7915 - 0.0538*T0 + 0.007*(T0**(2.0)) - 0.0023*S0)
        # ~0.0014 kg m-1 s-1
        Cd=0.82 # long cylinder
        W = np.sqrt(((kelpdiameter**2) * kelplength*np.pi*abs(dr)*g)/DENSw*2*Cd*kelparea)
        W = np.where(dr<0, W, -W)
       
        self.elements.terminal_velocity = W


    def update(self):
        """Update positions and properties of buoyant particles."""

        # Update element age
        self.elements.age_seconds += self.time_step.total_seconds()

        # Turbulent Mixing
        self.update_terminal_velocity()

        if self.get_config('processes:turbulentmixing') is True:
            self.vertical_mixing()

        # Horizontal advection
        self.advect_ocean_current()
       
        # Vertical advection
        if self.get_config('processes:verticaladvection') is True:
            self.vertical_advection()

        # Deactivate elements hitting land
        # self.deactivate_stranded_elements()
