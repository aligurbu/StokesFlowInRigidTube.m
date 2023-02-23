[![DOI](https://zenodo.org/badge/413033129.svg)](https://zenodo.org/badge/latestdoi/413033129)

# StokesFlowInRigidTube.m

- High-order weakly singular boundary element method code is presented for analyzing pressure-driven viscous flow in the rigid vessel with varying cross-sections.
- This code was developed for part of [my dissertation](https://www.researchgate.net/publication/355033649_Simulations_of_Red_Blood_Cell_Flow_by_Boundary_Integral_Methods) to simulate the red blood cell flow using boundary integral methods.
- This repository contains the code for the concepts and examples presented in Chapter 5.

## Numerical examples

|Straight vessel| Constricted vessel | Long constricted vessel |
| :-: | :-: | :-: |
|<img src="Results/ShortMicrocapillary_16El/ShortMicrocapillary_16El.png">|<img src="Results/RefinedConstrictedVessel_16El/RefinedConstrictedVessel_16El.png">|<img src="Results/LongConstrictedVessel_16El/LongConstrictedVessel_16El.png">|
| | Speed profile | |
|<img src="Results/ShortMicrocapillary_16El/SpeedProfileInsideVessel.png">|<img src="Results/RefinedConstrictedVessel_16El/SpeedProfileInsideVessel.png">|<img src="Results/LongConstrictedVessel_16El/SpeedProfileInsideVessel.png">|
|<img src="Results/ShortMicrocapillary_16El/SpeedProfileInsideVesselSideView.png">|<img src="Results/RefinedConstrictedVessel_16El/SpeedProfileInsideVesselSideView.png">|<img src="Results/LongConstrictedVessel_16El/SpeedProfileInsideVesselSideView.png">|
|<img src="Results/ShortMicrocapillary_16El/SpeedCenterlineVessel.png">|<img src="Results/RefinedConstrictedVessel_16El/SpeedCenterlineVessel.png">|<img src="Results/LongConstrictedVessel_16El/SpeedCenterlineVessel.png">|
| | Velocity vector fields | |
|<img src="Results/ShortMicrocapillary_16El/InletOutletVelocityProfile.png">|<img src="Results/RefinedConstrictedVessel_16El/InletOutletVelocityProfile.png">|<img src="Results/LongConstrictedVessel_16El/InletOutletVelocityProfile.png">|
|<img src="Results/ShortMicrocapillary_16El/VelocityProfileInsideVesselSideView.png">|<img src="Results/RefinedConstrictedVessel_16El/VelocityProfileInsideVesselSideView.png">|<img src="Results/LongConstrictedVessel_16El/VelocityProfileInsideVesselSideView.png">|
| | Norm of traction | |
|<img src="Results/ShortMicrocapillary_16El/NormTractionProfile.png">|<img src="Results/RefinedConstrictedVessel_16El/NormTractionProfile.png">|<img src="Results/LongConstrictedVessel_16El/NormTractionProfile.png">|
|<img src="Results/ShortMicrocapillary_16El/NormTractionLine.png">|<img src="Results/RefinedConstrictedVessel_16El/NormTractionLine.png">|<img src="Results/LongConstrictedVessel_16El/NormTractionLine.png">|
| | Pressure | |
|<img src="Results/ShortMicrocapillary_16El/NormPressureProfile.png">|<img src="Results/RefinedConstrictedVessel_16El/NormPressureProfile.png">|<img src="Results/LongConstrictedVessel_16El/NormPressureProfile.png">|
|<img src="Results/ShortMicrocapillary_16El/NormPresureLine.png">|<img src="Results/RefinedConstrictedVessel_16El/NormPresureLine.png">|<img src="Results/LongConstrictedVessel_16El/NormPresureLine.png">|
| | Norm of shear | |
|<img src="Results/ShortMicrocapillary_16El/NormShearProfile.png">|<img src="Results/RefinedConstrictedVessel_16El/NormShearProfile.png">|<img src="Results/LongConstrictedVessel_16El/NormShearProfile.png">|
|<img src="Results/ShortMicrocapillary_16El/NormShearLine.png">|<img src="Results/RefinedConstrictedVessel_16El/NormShearLine.png">|<img src="Results/LongConstrictedVessel_16El/NormShearLine.png">|

## Citation

    @phdthesis{gurbuz2021Thesis,
    title={Simulations of Red Blood Cell Flow by Boundary Integral Methods},
    author={G\"urb\"uz, Ali},
    year={2021},
    school={State University of New York at Buffalo}
    }
