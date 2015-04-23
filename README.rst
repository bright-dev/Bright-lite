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

Bright-lite Reactor Facility is a reactor modeling software that uses burnup, criticality, and 
transmutations matrix curves to determine input and output isotopic compositions. The reactor
can operate outside of the Bright-lite suite as a reactor however it will no long have start
up and blending capabilities. 

Bright-lite Fuel Fabrication Facility: A module that communicates directly with the Bright-lite
Reactor Facility that allows for the reactor facility to access blending functions that will 
determine the isotopic composition of fuel that the reactor will be using in order to match a
specific constraint. 

Bright-lite Reprocessing Facility: A module that seperates specific isotopes out of a material
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

 - in_commods: This field is a one or more than indicates that possible sources of 
   fuel for the reactor. The values in this field should be commodities that exist 
   inside of the simulation.  
------------
Something something results
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

