.. _Cyclus: http://www.fuelcycle.org/
.. _Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page

Welcome New Bright-lite User!
=============================
Bright-lite is a collection of modules for the Cyclus_ Fuel Cycle Simulator 
that allow for medium fidelity reactor modeling. There are currently 
three modules in the Bright-lite suite. 

- Bright-lite Reactor Facility
- Bright-lite Fuel Fabrication Facility
- Bright-lite Reprocessing Facility

**Bright-lite Reactor Facility**: A reactor modeling software that uses burnup, criticality, and 
transmutations matrix curves to determine input and output isotopic compositions. The reactor
can operate outside of the Bright-lite suite as a reactor however it will no long have start
up and blending capabilities. 

**Bright-lite Fuel Fabrication Facility**: A module that communicates directly with the Bright-lite
Reactor Facility that allows for the reactor facility to access blending functions that will 
determine the isotopic composition of fuel that the reactor will be using in order to match a
specific constraint. 

**Bright-lite Reprocessing Facility**: A module that seperates specific isotopes out of a material
that is sent to it. This currently requires an external text file to indicate which isotopes
should be seperated into which streams. The structure of this text file can be seen here_.

Bright-lite works in conjuction with the Cyclus_ Fuel Cycle Simulator. 

Bright-lite is currently only being actively supported for Ubuntu.

------------
Something here about how it works.
------------

Installation
------------
To use Bright-lite first you need to install it. Currently Bright-lite has
the following dependencies. 

- Cyclus_ Fuel Cycle Simulator

The following dependencies will be required in the future

- Eigen_

To install Bright-lite please follow these instructions.

1) Clone the Bright-lite repository from github.
2) Change directory into the Bright-lite directory using the following
   command. 
    cd Bright-lite
   	
3) Use the following command inside the Bright-lite directory.
    python install.py
   	
This will add the Bright-lite module to the cyclus environment, and allow
you to use Bright-lite in Cyclus simulations. 

------------
Using Bright-lite
------------
Bright-lite requires at least 6 inputs from the users to operate fully. While
there are several other inputs associated with the Bright-lite module all of 
these other inputs come with a default value. 

The six required inputs are

- **in_commods**: This field is a one or more than indicates that possible sources of 
  fuel for the reactor. The values in this field should be commodities that exist 
  inside of the simulation.  
- **out_commod**: This field should be filled out with the cyclus commodity that will
  connect the reactor facility to the facility that will be directly handling the 
  waste.
- **libraries**: This is a one or more field that indicates the Bright-lite library 
  the reactor will be using. Note: Adding additionally libraries to this list
  will enable the 'library interpolation'_ capabilities in Bright-lite but also
  requires that the user input parameters and values to be interpolated upon. The
  interpolation feature is intended for advanced users. 
- **target_burnup**: This field indicates to the reactor what the target burnup for the 
  reactor will be. If this is set to 0, the reactor will operate in forward mode. If 
  this value is not set to zero the Bright-lite reactor must be connected to a
  Bright-lite fuel fabrication facility.
- **core_mass**: This field indicates the total mass of fuel inside of the core. This mass
  does not include structural components, it is only the mass of fuel to be burned.
- **generated_power**: This indicates the total thermal generating power of the core. 
  The electrical generated_power will be this value times the effiency of the reactor
  (a input set to default at 33% but is user adjustable).
  
------------
Something something results
------------

------------
Library Interpolation
------------
.. _'library interpolation'

The libraries used in Bright-lite are often associated with several parameters. For example
an LWR reactor library might have parameters for burnup, and enrichment. If as a user, you
require a different value for these parameters there are two possible methods for obtaining it
First, a new library can be generated externally from Bright-lite using tools available (XSGEN
for example). It is also possible to create a dynamic library that matches your desired parameters
using Bright-lite's built in library interpolation tool.

This tool is used using two key components in the Bright-lite input schema.

- **libraries** 
- **interpolation_pairs**

For example, there may be two LWR libraries that fit into an LWR library suite. 

- Reactor 1
 - Burnup: 50 MWd/kgIHM
 - Enrichment: 5% U235
- Reactor 2
 - Burnup: 30 MWd/kgIHM
 - Enrichment: 3.3% U235
 
If a new library with the following parameters is desired

- Dynamic Reactor
 - Burnup: 40 MWd/kgIHM
 - Enrichment: 4% U235


------------
Format of Reprocessing Plant Text File
------------
.. _here:
 
	BEGIN
	 isotope1n fraction1n
	 
	 isotope2n fraction2n
	 
	 ...
	 
	 isotopeN fractionN
	 
	END
	
	BEGIN
	
	 isotope1k fraction1k
	 
	 isotope2k fraction2k
	 
	 ...
	 
	 isotopeK fractionK
	 
	END

------------

