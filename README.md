
## Survey Simulator 2.0  Documentation


The Survey Simulator (sometimes compiled as `Driver` or `SSim` or `SurveySimulator`) takes a model of 
the orbits and physical parameters of objects in the outer solar system and determines which model objects a 
particular survey would detect and track. 

The user provides the model and files that describe the characterization of the survey of interest.
For a description of the how to generate a model see the `F95` or `F77` or `python` subdirectoris, 
depending on your preferred language.  Each of these contains and `example` subdirectory that 
provides examples / clues have how to use the Survey Simulator in that language.

For descriptions of the formats of the characterization files see the `Surveys` sub-directory.

---
### Licence
This software is released under the terms of the European Union Public
Licence v1.1 (EUPL v.1.1). See files eupl1.1.-en_0_0.pdf (Preamble) and
eupl1.1.-licence-en_0.pdf (detailed description of the licence).

This source code is provided as is, with no warranty of any kind.  The user takes full responsiblity for any damage to system, and
for any scientific conclusion drawn.

### Contact
The primary contact for the Survey Simulator code is:  
* Jean-Marc Petit: Jean-Marc.Petit@normalesup.org

### Acknowledgement

Cite **Petit, J.-M., et al., AJ, Vol 142 ID 131 (2011)** if you make 
use of the SurveySimulator, or the CFEPS L7SyntheticModel-v09 Kuiper belt model.

If you make use of the Survey characterizations and detections please cite
the appropriate survey paper:  
* CFEPS:
   * _Petit, J.-M., et al., AJ, Vol 142 ID 131 (2011)_
* OSSOS:
   * _Bannister et al (2016) AJ, 152, 70_  
   * _Bannister et al (2018) ApJS, 236, 18_
* HiLat:
   * _Petit et al (2017), AJ, 153, 236_
* MA Survey:
   * _Alexandersen et al (2016), AJ, 152, 111_

---
## Overview
This package of provides programmes and data aiming at simulating large-scale,
well-calibrated KBO surveys like CFEPS and OSSOS. The ultimate goal is to 
compare:

A. An orbital and size distribution model of the outer Solar System  
**to**  
B. Those same distributions as observed by a survey.  

The SurveySimulator provides the information need to allow a quantitative 
comparison of the parameter distributions provided by a set of
detected Kuiper belt objects that were actually found in a survey to what 
the survey would have found give a model of the Kuiper belt. 

The Survey Simulator itself provides enables the 'biasing' of a model of
the Kuiper Belt in a way that mimics the observational biases in a survey.

### Quick Start
Look in `Simulator/F{77|95}/fortran/exmaple/README.parametric` for an example of how to run Simulator.

####  Survey Characterization
The parameters that describe the observational biases are provided by the
individual surveys. The more precise the characterization of a given survey's
observational biases, the more accurate the results of the SurveySimulator
will be. Characterizations include the observing dates and pointing
locations of the discovery observations along with an accurate
estimate of the detection efficiency as a function of apparent
magnitude in the filter of the survey and an object's rate of motion
on the sky at the time of each potentail observation.

### The Simulation Process
The OSSOS SurveySimulator _draws_ objects from a _model_
and computes if the model object would have been detected by the
survey(s) defined by the survey configuration files. 
Using the survey configuration settings the Simulator computes 
if the model object is inside one of the imaged fields, bright enough to be detected, 
and within the detectable range of rates-of-motion on the sky.
The _model_ objects that would have been detected are call simulated
_detections_.  The SurveySimulator further classifies, based on the 
input characterization, which of those simulated _detections_ would have been
_tracked_ to high-precision orbits.  The SurveySimulator provides to the 
user this list of _detected_ and _tracked_ _model_ objects. 

The user should then compare the _tracked_ model objects to a given survey's
detections via some statistical method (see 
_Lawler et al., Frontiers in Astronomy and Space Sciences, Volume 5, id.14, 2018_). 
This last step is not part of the simulator itself, leaving users to 
decide on their own statistical approach.

Note that list of real detections from a survey are not required to 
RUN the SurveySimulator.  The true detections from the survey do not 
influence the output of the SurveySimulator. The true detections can be used
after execution to measure the validity of the model via comparison 
between the orbital and H-mag distribution of the _tracked_ 
detection, given the model, and the real detections from the Survey whose 
characterization was used as input into the SurveySimulator.

---
## Package Contents  
The SurveySimulator-2.0 release consists of a Driver program, subroutines that
define trans-neptunian or other outer Solar System objects (either from a 
parametric model or a lookup table), data that describe the CFEPS, OSSOS, MA and some other survey 
characterizations, and a list of the real classified objects discovered during 
the CFEPS project (to be compared to the output of the Survey Simulator).

The common source codes for the survey simulator, available as F77 or F95.
The F77/fortan and F95/fortran directories contain the Fortran code that is the simulator.  
There are also F77/python and F95/python that provide examples of building Python callable
modules from the two fortran branches.
The F95 and F77 versions are identitical, just different language implementation.
 See the README files for details.  

The architecture of the SurveySimulator is to provide a subroutine named GiMeObj
(for examples see the `InnerHotModel.f{95}` or `ReadModelFromFile.f{95}` source code)
when compiling the Driver program.  The simulator is run by calling *Driver*.

See Simulator/F95/fortran/example or Simulator/F77/fortan/example for examples of how
to compile (e.g. `make InnerHotModel`) a survey simulator.

### Simulator/

#### Simulator/fortran/{F95|F77}/SurevySubs.f
This contains the source code for the Survey Simulator.  in particular Detos1 which determines, based on
the given orbit and the survey characterization area, which sources are detected.

The Detos1 subrouting is described in the source code but reproduced here to make you aware of the inputs and outputs.

```fortran

subroutine Detos1 (o_m, jday, hx, color, gb, ph, period, amp, surnam, seed, &
          debug, &
       flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, mt, &
       jdayp, ic, surna, h_rand, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine determines if a given object is seen by the survey
! described in the directory \verb|surnam|.
! An object is described by its ecliptic (J2000) barycentric osculating
! elements given at time \verb|jday|.
! This version uses polygons to describe the footprint of the block on
! the sky.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J.-M. Petit  Observatoire de Besancon
! Version 1 : January 2006
! Version 2 : May 2013
! Version 3 : March 2016
! Version 4 : May 2016
!             Changed API to remove size of arrays, added parameter
!             statement to define array sizes (in include file).
!             Continue looping on pointings until object is detected,
!             characterized and tracked. Don't stop at first detection.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     o_m   : orbital elements of object (t_orb_m)
!     jday  : Time of elements [JD] (R8)
!     hx    : Absolute magnitude of object in 'x' band, what ever this is (R8)
!     color : Array of colors (128*R8)
!                IACHAR(FILTER)+1 : FILTER - X
!     gb    : opposition surge factor G, Bowell formalism (R8)
!     ph    : phase of lightcurve at epoch jday [rad] (R8)
!     period: period of lightcurve [day] (R8)
!     amp   : amplitude of lightcurve [mag] (R8)
!     surnam: Survey directory name (CH)
!
! OUTPUT
!     seed  : Random number generator seed (I4)
!     flag  : Return flag (I4): 
!                0: not found 
!                1: found, but not tracked
!                2: found and tracked
!                3: characterized, but not tracked
!                4: characterized and tracked
!     ra    : Right ascension at detection [rad] (R8)
!     dec   : Declination at detection [rad] (R8)
!     d_ra  : Right ascension rate [rad/day] (R8)
!     d_dec : Declination rate [rad/day] (R8)
!     r     : Sun-object distance [AU] (R8)
!     delta : Earth-object distance [AU] (R8)
!     m_int : Intrinsic apparent magnitude, in x-band (R8)
!     m_rand: Averaged randomized magnitude, in detection filter (R8)
!     eff   : Actual efficiency of detection (R8)
!     isur  : Identification number of survey the object was in (I4)
!     mt    : Mean anomaly at discovery [rad] (R8)
!     jdayp : Time of discovery [JD] (R8)
!     ic    : Index of color used for survey (I4)
!     surna : Detection survey name (CH10)
!     h_rand: Absolute randomized magnitude, in detection filter (R8)
!     ierr  : error flag
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
'''
  
```
#### Simulator/{F95|F77}/fortran/ReadModelFromFile.f 
This contains the source code for a "GiMeObj" routine that reads 
an (orbital+size+colour+lightcurve) model from a file (lookup table).
Use `make ReadModelFromFile` to build a `Driver` program that can be used 
to run simulation of observing the Kuiper belt model described in `ReadModelFromFile.f`. 

#### Simulator/{F95|F77}/fortran/InnerHotModel.f 
Contains the source code for a GiMeObj routine that generates
objects according to some parametric prescription. 
Use `make InnerHotModel` to build a `Driver` program that can be used 
to run simulation of observing the Kuiper belt model described in `InnerHotModel.f`. 
See the `example` directories to understand how to use the Simulator.

#### Simulator/F{77|95}/fortran/example
Examples of how to run the Simulator.

### Simulator/{F95|F77}/python
Contains Makefiles that build a python module that can be used to Drive the 
survey simulator.  We provide a F77 and F95  module building capacity. 
See `README.python` for details.

### Python
This directory contains the source code for a GiMeObj routine that generates
objects according to some parametric presciption, along with a Makefile to
generate the executable. This example uses the Driver.py Python program as
driver instead of the usual Fortran program Driver.f.

### Models/
#### Models/L7model-3.0-9.0  
Contains the CFEPS L7 model of the debiased Kuiper Belt) as an example input file for 
the survey simulator. See `ReadModelFromFile.f{95}` for example code that uses 
this file as input.
#### Models/InnerHot.in 
Contains the parameters needed by the `InnerHotModel.f{95}` version of the `GiMeObj` 
subroutine.

### Surveys
The configuration files that configure the Simulator for various Survey characterizations. 
#### Surveys/cfeps
Pointing history and efficiency functions for 
each cfeps block, including the cfeps presurvey 'block'; it also contains 
the list of real detected objects (cfeps.detections).  Other calibrated 
surveys could be substituted (or added) once their characterization is 
specified in the correct format.

#### Surveys/OSSOS
Pointing history and efficiency functions for
each OSSOS block; it also contains the list of real detected objects
including their dynamical class.
 
#### Surveys/All_Surveys
Union of all surveys: CFEPS, HiLat, Alexandersen and
OSSOS, but SEE IMPORTANT INFORMATION in that subdir's README.allsurveys file
about how care must be taking when combining the surveys.

#### Surveys/All_r_Surveys
Similar to _All_Surveys_ but restricted to objects detected with the
MegaPrime r filter (OSSOS, Alexandersen, HiLat and L3h block from CFEPS).

---
