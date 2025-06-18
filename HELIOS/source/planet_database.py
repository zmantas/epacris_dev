# ==============================================================================
# Module with a database of planetary system parameters used in HELIOS
# Copyright (C) 2020 - 2022 Matej Malik
# ==============================================================================
# This file is part of HELIOS.
#
#     HELIOS is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     HELIOS is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You find a copy of the GNU General Public License in the main
#     HELIOS directory under <license.txt>. If not, see
#     <http://www.gnu.org/licenses/>.
# ==============================================================================

try:
    from HELIOS_vulcan.source import phys_const as pc
    
except ImportError:
    try:
        from source import phys_const as pc
    except ImportError:
        import phys_const as pc



class Planet(object):

    def __init__(self, R_p, g_p, a, T_star, R_star, g_star, metal_star, R_p_unit,):

        self.R_p = R_p
        if R_p_unit == "R_Earth":
            self.R_p *= pc.R_EARTH / pc.R_JUP

        self.pl_radE = R_p
        self.g_p = g_p
        self.a = a
        self.T_star = T_star
        self.R_star = R_star
        self.g_star = g_star
        self.metal_star = metal_star


planet_lib = {}

# Units are [R_p]=R_Earth or R_Jup, [g]=cm s^-2 or [g]=log(cm s^-2), [a]=AU, [R_star]=R_Sun, [T_star]=K

planet_lib["GJ_1214b"] = Planet(R_p=2.85, R_p_unit="R_Earth",
                                g_p=760,
                                a=0.01411,
                                T_star=3026,
                                R_star=0.216,
                                g_star=4.944,
                                metal_star=0.39
                                )  # references: Harpsoe et al. (2013)

planet_lib["wasp-39b"] = Planet(R_p=14.23532, R_p_unit="R_Earth",
                                g_p=430,
                                a=0.0486,
                                T_star=5400.,
                                R_star=0.895,
                                g_star=4.3893300,
                                metal_star=-0.030
                                )  # references: 

planet_lib["GJ1214_0"] = Planet(R_p=2.0, R_p_unit="R_Earth",
                                g_p=2454,
                                a=0.000897,
                                T_star=3026,
                                R_star=0.216,
                                g_star=4.944,
                                metal_star=0.39
                                )  # references: Harpsoe et al. (2013)

planet_lib["HD85512_0"] = Planet(R_p=2.0, R_p_unit="R_Earth",
                                g_p=2454,
                                a=0.00580,
                                T_star=4305,
                                R_star=0.69,
                                g_star=4.596,
                                metal_star=0.39
                                )  # references: Harpsoe et al. (2013)

planet_lib["K381_0"] = Planet(R_p=2.0, R_p_unit="R_Earth",
                                g_p=2454,
                                a=0.0273,
                                T_star=6407,
                                R_star=1.466,
                                g_star=4.215,
                                metal_star=0.39
                                )  # references: Harpsoe et al. (2013)

planet_lib["HD_209458b"] = Planet(R_p=1.380, R_p_unit="R_Jupiter",
                               g_p=930,
                               a=0.04747,
                               T_star=6117,
                               R_star=1.162,
                               g_star=4.368,
                               metal_star=0.02
                               )  # references: Southworth (2010)

planet_lib["55Cnce"] = Planet(R_p=1.947, R_p_unit="R_Earth",
                               g_p=2222,
                               a=0.016,
                               T_star=5214,
                               R_star=0.98,
                               g_star=4.458,
                               metal_star=0.35
                               )  # references: Southworth (2010)

planet_lib["TOI561b"] = Planet(R_p=1.425, R_p_unit="R_Earth",
                               g_p=967.14,
                               a=0.01055,
                               T_star=5326,
                               R_star=0.843,
                               g_star=4.50,
                               metal_star=-0.40
                               )  # 

planet_lib["K2-141b"] = Planet(R_p=1.51, R_p_unit="R_Earth",
                               g_p=2187.76,
                               a=0.007478,
                               T_star=4599,
                               R_star=0.68,
                               g_star=4.62,
                               metal_star=0.
                               )  #

planet_lib["TOI1807b"] = Planet(R_p=1.37, R_p_unit="R_Earth",
                               g_p=1344.568,
                               a=0.0120,
                               T_star=4757,
                               R_star=0.69,
                               g_star=4.55,
                               metal_star=-0.04
                               )  # references: Southworth (2010)

planet_lib["TOI-2431b"] = Planet(R_p=1.406, R_p_unit="R_Earth",
                               g_p=1788.,
                               a=0.00633,
                               T_star=4099,
                               R_star=0.63,
                               g_star=4.59,
                               metal_star=0
                               )  # references: all over
							  
planet_lib["HD213885b"] = Planet(R_p=1.75, R_p_unit="R_Earth",
                               g_p=2754.,
                               a=0.02012,
                               T_star=5978,
                               R_star=1.101,
                               g_star=4.3827,
                               metal_star=-0.04
                               )  # references: all over	

planet_lib["TOI2260b"] = Planet(R_p=1.62, R_p_unit="R_Earth",
                               g_p=1308.438,
                               a=0.0097,
                               T_star=5534,
                               R_star=0.94,
                               g_star=4.62,
                               metal_star=0.22
                               )  # references: all over								   

planet_lib["HD219134b"] = Planet(R_p=1.5, R_p_unit="R_Earth",
                               g_p=1859.38,
                               a=0.037,
                               T_star=4699.0,
                               R_star=0.778,
                               g_star=4.567,
                               metal_star=0.11
                               )  # exoplanets eu/wiki		

planet_lib["TOI1634b"] = Planet(R_p=1.789, R_p_unit="R_Earth",
                               g_p=2359.9,
                               a=0.01545,
                               T_star=3550.0,
                               R_star=0.45,
                               g_star=4.833,
                               metal_star=0.23
                               )  # exoplanets eu/wiki	

planet_lib["TOI2427b"] = Planet(R_p=1.80, R_p_unit="R_Earth",
                               g_p=2722.2,
                               a=0.0202,
                               T_star=4072.0,
                               R_star=0.65,
                               g_star=4.62,
                               metal_star=0.0
                               )  # exoplanets eu/wiki	

planet_lib["ISAAC"] = Planet(R_p=1.0, R_p_unit="R_Jupiter",
                               g_p=1000,
                               a=17.505,
                               T_star=6000,
                               R_star=1,
                               g_star=4.43,
                               metal_star=0.0
                               )  # exoplanets eu/wiki	

planet_lib["HEarth"] = Planet(R_p=1.0, R_p_unit="R_Earth",
                               g_p=980,
                               a=1.0,
                               T_star=5778,
                               R_star=1,
                               g_star=4.44,
                               metal_star=0.0
                               )  # Earth parameters


planet_lib["Alpha_Ab"] = Planet(R_p=1, R_p_unit="R_Jupiter",
                               g_p=1000,
                               a=1.0,
                               T_star=5790,
                               R_star=1.2175,
                               g_star=4.32,
                               metal_star=0.2
                               )  # 
							   
planet_lib["Alpha_Ac"] = Planet(R_p=1, R_p_unit="R_Jupiter",
                               g_p=1000,
                               a=2.0,
                               T_star=5790,
                               R_star=1.2175,
                               g_star=4.32,
                               metal_star=0.2
                               )  #

planet_lib["Alpha_Ae"] = Planet(R_p=1, R_p_unit="R_Jupiter",
                               g_p=1000,
                               a=1.9,
                               T_star=5790,
                               R_star=1.2175,
                               g_star=4.32,
                               metal_star=0.2
                               )  #  							   

planet_lib["Alpha_Ad"] = Planet(R_p=0.9, R_p_unit="R_Jupiter",
                               g_p=1530,
                               a=1.9,
                               T_star=5790,
                               R_star=1.2175,
                               g_star=4.32,
                               metal_star=0.2
                               )  # A puffy version

planet_lib["Alpha_Ax"] = Planet(R_p=1., R_p_unit="R_Jupiter",
                               g_p=2000,
                               a=2.0,
                               T_star=5790,
                               R_star=1.2175,
                               g_star=4.32,
                               metal_star=0.2
                               )  # A puffy version
							   
planet_lib["Alpha_Af"] = Planet(R_p=1.1, R_p_unit="R_Jupiter",
                               g_p=1000,
                               a=1.9,
                               T_star=5790,
                               R_star=1.2175,
                               g_star=4.32,
                               metal_star=0.2
                               )  # A puffy version							  

planet_lib["K2-18b"] = Planet(R_p=2.61, R_p_unit="R_Earth",
                               g_p=1244,
                               a=0.1591,
                               T_star=3457,
                               R_star=0.44,
                               g_star=4.79,
                               metal_star=0.12
                               )  # A puffy version		

if __name__ == "__main__":
    print("This module stores information about planetary systems. No guarantee that anything here is remotely correct.")