.. _Cyclus: http://www.fuelcycle.org/
.. _Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page

*****************************
Welcome New Bright-lite User!
*****************************
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

============
Installation
============
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

============
Using Bright-lite
============
Bright-lite requires at least 6 inputs from the users to operate fully. While
there are several other inputs associated with the Bright-lite module all of 
these other inputs come with a default value. 

-----------
The six required inputs are
-----------

- **in_commods**: This field is a one or more than indicates that possible sources of 
  fuel for the reactor. The values in this field should be commodities that exist 
  inside of the simulation.  
- **out_commod**: This field should be filled out with the cyclus commodity that will
  connect the reactor facility to the facility that will be directly handling the 
  waste.
- **libraries**: This is a one or more field that indicates the Bright-lite library 
  the reactor will be using. Note: Adding additionally libraries to this list
  will enable the library interpolation_ capabilities in Bright-lite but also
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

-----------
Operational Modes
-----------
Bright-lite has two operational modes. The mode is indicated using inputs to the 
Bright-lite reactor module. Forward mode is chosen by setting the reactor *target_burnup*
mode to be equal to 0. Blending mode requires the *target_burnup* field to be a non negative
value. Additionally, blending mode requires the Bright-lite ReactorFacility to be connected to a
Bright-lite FuelfabFacility. 
 
^^^^^^^^^^^^
Forward Mode
^^^^^^^^^^^^
In forward mode Bright-lite accepts a fuel composition and burns it foward in time to match
a target. It does this by advancing the fluence of each batch in the core until the target is met.

Currently forward mode works only with *criticality* and *burnup* targets

^^^^^^^^^^^^
Blending Mode
^^^^^^^^^^^^
As stated above using the blending function in Bright-lite requires connecting a Bright-lite
ReactorFacility to a Bright-lite FuelfabFacility. The *in_commods* field of the Bright-lite
reactor should include all of the fuel fabrication facilities that the reactor can be connected to. 

Currently there are only two blending modes available in Bright-lite. These modes are described by
a target-constraint pairing. The two available pairs currently are:
 
1) Burnup - Criticality: The blender will create a fuel that meets a target burnup when criticality is equal to the given constraint. This set of constraints only requires a non negative number to be entered into the *target_burnup* field. 
2) Burnup - Conversion Ratio: The blender will create a fuel that meets a target burnup when conversion ratio is equal to the given constraint. This is achieved by setting the *CR_target* input field of the reactor to 
be equal to a number greater than 0 (note that there is no upper bound limit in the code for this this but physically it should not exceed 2). Additionally the *target_burnup* field must be a non negative value for this to work. 

------------
Something something results
------------

------------
Library Interpolation
------------
.. _interpolation:

The libraries_ used in Bright-lite are often associated with several parameters. For example
an LWR reactor library might have parameters for burnup, and enrichment. If as a user, you
require a different value for these parameters there are two possible methods for obtaining it
First, a new library can be generated externally from Bright-lite using tools available (XSGEN
for example). It is also possible to create a dynamic library that matches your desired parameters
using Bright-lite's built in library interpolation tool.

This tool is used using two key components in the Bright-lite input schema.

**libraries** 
- To enable library interpolation here simple add more than one library to the field. This is done
  simply by adding another val to the input field. That is...::
  
  <val>extLWR</val>
  represents a reactor library using just the *extLWR* library. However by adding another library::
  
  <val>extLWR</val>
  <val>lowLWR</val>
  Bright-lite will make a new library based on the interpolation pairs and the values inside of 
  these two libraries. 
**interpolation_pairs**
- Once two or more libraries have been selected at least one interpolation pair will need to be added. 
  An interpolation pair is a <"Parameter", Value> pair. The parameter represents a common parameter 
  shared by the libraries, and the value is the target value for the new dynamic library in that 
  parameter. 

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

The following xml should be added to the reactor archetype.
::

 <libraries>
  <val>Reactor 1</val>
  <val>Reactor 2</val>
 </libraries>
 <interpolation_pairs>
  <key>BURNUP</key>
  <val>40</val>
  <key>ENRICHMENT</key>
  <val>4</val>
 </interpolation_pairs>

------------
Available Libraries
------------
.. _libraries:

Recommended Libraries

- lowLWR - A standard PWR library.
 - Enrichment: 2.2 %U235
 - Burnup: 20 MWd/kgIHM
 - PNL: 0.903
 - Batches: 3

- standLWR
 - Enrichment: 3.3 %U235
 - Burnup: 33 MWd/kgIHM
 - PNL:0.911 
 - Batches: 3

- extLWR
 - Enrichment: 5% U235
 - Burnup: 50 MWd/kgIHM
 - PNL: 0.883
 - Batches: 3

- BWRMOX
 - Burnup:
 - PNL:
 - Batches: 
 
- PWRMOX
 - Burnup:
 - PNL:
 - Batches:

- DUPIC
 - Burnup:
 - PNL:
 - Batches:
 
- FR25
 - Burnup:
 - Conversion Ratio: 0.25
 - PNL:
 - Batches:
 
- FR25MOX
 - Burnup:
 - Conversion Ratio: 0.25:
 - PNL:
 - Batches:
 
- FR50
 - Burnup:
 - Conversion Ratio: 0.5:
 - PNL:
 - Batches:
 
- MOXMA
 - Burnup:
 - PNL:
 - Batches:
 
Additional Libraries

- E5_50
 - Enrichment
 - Burnup
 - PNL
 - Batches
 
- E5_60
 - Enrichment
 - Burnup
 - PNL
 - Batches
 
- E7_100
 - Enrichment
 - Burnup
 - PNL
 - Batches
 
- E9_100
 - Enrichment
 - Burnup
 - PNL
 - Batches
 
------------
Format of Reprocessing Plant Text File
------------
.. _here:
::

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

