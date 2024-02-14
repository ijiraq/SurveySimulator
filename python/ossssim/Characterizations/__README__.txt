
ObsSummary README for the OSSOSv11 release 

-------------------------------------------------------------------
Contents

Directories and files in ObsSummary directory

All_Surveys     subdir  The full sample, including all detections from OSSOS,
                        CFEPS, HiLat and MA.
                        WARNING: 6 objects from CFEPS have been (blindly, not
                        targeted) re-observed by OSSOS. These objects appear
                        only once in All_Surveys.detections file. We have made
                        the choice to list them as OSSOS detections, with their
                        detection circumstances (in particular the filter - r -
                        and the magnitude) at time of OSSOS. The detection
                        circumstances at time of CFEPS can be found in the
                        CFEPS/CFEPS.detections file.
All_r_Surveys   subdir  All detections from OSSOS, CFEPS, HiLat and MA acquired
                        in MegaPrime r filter. This comprise all OSSOS, MA,
                        HiLat and the L3h block from CFEPS. In this sub-sample,
                        there is no duplicate object (object detected in 2
                        different surveys).
CFEPS           subdir  All detections from CFEPS. This includes the objects
                        detected in CFEPS and re-detected in OSSOS. These
                        objects are listed here with their CFEPS detection
                        circumstances.
OSSOS           subdir  All detections from OSSOS and only OSSOS.
OSSOS-MA        subdir  All detections from OSSOS and MA.

Correspondence.list     Matches MPC names to Surveys.  Human readable format
OSSOSv11ScatterPlots.pdf Plot of uncertainties and orbital elements
README_ObsS_v11.txt     this file
Uncharact.detections    Orbital elements and classes for Uncharacterised 
                        objects that were tracked outside of the discovery
                        dark run. Same format as All_Surveys_v11.detections
Uncharact-nt.detections Orbital elements and classes for Uncharacterised and
                        non-tracked objects. Same format as above.

Detection files for the various characterized samples can be found in the above
directories. File names are <sample>/<sample>[_v11].detections and
<sample>/<sample>[_v11].CDS. As an example, the full sample files are in
All_Surveys as:

All_Surveys_v11.detections Orbital elements and classes. This is a
                        space-separated value file, with name of the column
                        given on the last comment line (starting with #) at the
                        start of the file. See header of CDS file for a
                        description of the columns.
			Contains OSSOS, CFEPS, HiLat and MA.
All_Surveys_v11.CDS        Same as All_Surveys_v11.detections in CDS format

-------------------------------------------------------------------
File Formats

-The .detections file gives an averaged magnitude (and band, and uncertainty)
    for each object DURING THE DISCOVERY NIGHT'S TRIPLE (this is the only
    magnitude that matters for the characterization).  The 'surmised' absolute
    magnitude (that is, the H_r magnitude surmised using the average m_r and
    the discovery geometry) is also given (it has the same error). The dynamical
    classification of the objects is given in the first several columns; see the
    wiki page and file headers for details.

-------------------------------------------------------------------
Caveats/clarifications for the v11 release

- There are 49/37/84/67/147/87/67/55/105/146 characterized objects in
    E/O/L/H/P/M/S/T/C/D block, respectively.  
    The unclassified uo3eXX, uo3oXX, uo3lXX, uo4hXX, uo5pXXX, uo5mXX, uo5sXX,
    uo5tXX, uo5cXXX, uo5dXXX objects are not part of the release but can be
    obtained by request to coreossos@taos.asiaa.sinica.edu.tw  although
    they are listed in the file in the top level directory
- Orbital classifications are 'secured' according to SSBN08 nomenclature. 
  The detections file gives an 'S'=secure or 'I'=insecure rating
- Orbital classifications for less than 10% of the objects is still insecure.
- Note that usage of the CFEPS survey will require the user make choices related
    to what they believe the g-r colour distributions are...
- All products for CFEPS, Hilat and MA (abbrevation for the two survey blocks
    mal and mah done in Mike Alexandersen's thesis) were previously published,
    and are being provided in this release to the OSSOS team for convenience. 
