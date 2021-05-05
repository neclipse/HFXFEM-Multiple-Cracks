# HFXFEM-Verified-Singlecrack

# Author info
Chang Huang : chuan25@lsu.edu, huangchang73@gmail.com

# Affiliation
Geomechanics Group leaded by Dr. Shengli Chen (shenglichen@lsu.edu)
Department of Civil and Environmental Engineering
Louisiana State University

# Funds support
The development of this framework is supported by the ACS Petroleum Research Fund, American Chemical Society (PRF# 56743-DNI9), Industrial Ties Research Subprogram, Board of Regents, Louisiana [LEQSF(2016-19)-RD-B-02], and Economic Development Assistantships, Louisiana State University (Award No. 000408).

# Code development description
The project was developed for poroelastic material with discontinuities. It is designed to study hydraulic fracturing in porous rock formation, with/without existing natural fractures. The underlying method for discontinuity modeling is extended finite element method (XFEM) (Moës et al. 1999). The master branch is verified for single crack against latest analytical solution (Dontsov 2017).

This Matlab package is developed using the Object-oriented paradiagm. The core structure is inherited from my other [repository](https://github.com/neclipse/FEA-in-Matlab-NSMOOM).Please consider to use the two algorithm flowcharts and one class aggregation map for quick understanding of the package.
 
More detailed introduction of the code can be seen in Chapter 4 of my dissertation at LSU Digital Commons, *Modelling Hydraulic Fracturing Initiation and Propagation in Porous Rock Formations*. The theoretical details of XFEM formulation and application examples are also available in our latest [publication](https://www.onepetro.org/journal-paper/SPE-204476-PA).

# References

- Dontsov, E. V. 2017. An approximate solution for a plane strain hydraulic fracture that accounts for fracture toughness, fluid viscosity, and leak-off. International Journal of Fracture 205 (2): 221-237. http://doi.org/10.1007/s10704-017-0192-4.
- Moës, Nicolas, Dolbow, John, and Belytschko, Ted. 1999. A finite element method for crack growth without remeshing. International journal for numerical methods in engineering 46 (1): 131-150. http://doi.org/10.1002/(SICI)1097-0207(19990910)46:1<131::AID-NME726>3.0.CO;2-J.
