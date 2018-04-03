/* mpcci_quantities.h
 * 
 * This file was automatically generated by the Perl script "qc"
 *    Source file: "mpcci_quantities.def"
 *    Compiled by: pbayrasy@catalan, Wed May  6 17:47:43 2015
 * 
 * Quantity id (QID) definitions: The QID consists of the ...
 *       - dimension                (DIM) [1...64]
 *       - interpolation            (IOK) [0... 1]
 *       - extrapolation            (XOK) [0... 1]
 *       - field interpolation type (INT) [1...15]
 *       - physical meaning         (PHY) [1...15]
 *       - 4 validity bits          (VAL) : 0=global, quantity valid for Volume/Face/Line/Point coupling
 *       - tensor order             (ORD) [0... 7] 0=scalar/1=vector/2=tensor/3=3dr order tensor...
 *       - ramping                  (ROK) [0... 1]
 *       - signature                (SIG) [0..255]
 *    and is specified within a single unsigned.
 * 
 *    nibbles: |0000|0000|0|000|0000|0000|0000|00|00|0000|
 *             |---SIG---|R|ORD|-PHY|-VAL|-INT|XI|---DIM-|
 * 
 *    $Id: mpcci_quantities.h 1568 2015-05-06 15:51:33Z pbayrasy $
 */

#define MPCCI_QID_INT_SWITCH               0x01000b01   /* Control switch(Int) */
#define MPCCI_QID_REAL_SWITCH              0x02000b01   /* Control switch(Real) */
#define MPCCI_QID_PHYSICAL_TIME            0x03000a01   /* Physical time */
#define MPCCI_QID_TIMESTEP_SIZE            0x04000a01   /* Time step size */
#define MPCCI_QID_TIMESTEP_COUNT           0x05000b01   /* Time step number */
#define MPCCI_QID_ITERATION_COUNT          0x06000b01   /* Iteration number */
#define MPCCI_QID_GLOBAL_RESIDUAL          0x07000b01   /* Global residual */
#define MPCCI_QID_REF_PRESSURE             0x0b000b01   /* Reference pressure */
#define MPCCI_QID_VOLTAGE1                 0x15000b01   /* Electric voltage - phase 1 */
#define MPCCI_QID_VOLTAGE2                 0x16000b01   /* Electric voltage - phase 2 */
#define MPCCI_QID_VOLTAGE3                 0x17000b01   /* Electric voltage - phase 3 */
#define MPCCI_QID_VOLTAGE4                 0x18000b01   /* Electric voltage - phase 4 */
#define MPCCI_QID_CURRENT1                 0x19000b01   /* Electric current - phase 1 */
#define MPCCI_QID_CURRENT2                 0x1a000b01   /* Electric current - phase 2 */
#define MPCCI_QID_CURRENT3                 0x1b000b01   /* Electric current - phase 3 */
#define MPCCI_QID_CURRENT4                 0x1c000b01   /* Electric current - phase 4 */
#define MPCCI_QID_MO_POSITION              0x1f170b03   /* Moving obstacle CG position */
#define MPCCI_QID_MO_ANGLE                 0x20170b03   /* Moving obstacle CG angle */
#define MPCCI_QID_MO_VELOCITY              0x21150b03   /* Moving obstacle CG velocity */
#define MPCCI_QID_MO_OMEGA                 0x22150b03   /* Moving obstacle CG angular velocity */
#define MPCCI_QID_UGS_00                   0x28000b01   /* Global scalar 00 */
#define MPCCI_QID_UGS_01                   0x29000b01   /* Global scalar 01 */
#define MPCCI_QID_UGS_02                   0x2a000b01   /* Global scalar 02 */
#define MPCCI_QID_UGS_03                   0x2b000b01   /* Global scalar 03 */
#define MPCCI_QID_UGS_04                   0x2c000b01   /* Global scalar 04 */
#define MPCCI_QID_UGS_05                   0x2d000b01   /* Global scalar 05 */
#define MPCCI_QID_UGS_06                   0x2e000b01   /* Global scalar 06 */
#define MPCCI_QID_UGS_07                   0x2f000b01   /* Global scalar 07 */
#define MPCCI_QID_UGV_00                   0x28100b03   /* Global vector 00 */
#define MPCCI_QID_UGV_01                   0x29100b03   /* Global vector 01 */
#define MPCCI_QID_UGV_02                   0x2a100b03   /* Global vector 02 */
#define MPCCI_QID_UGV_03                   0x2b100b03   /* Global vector 03 */
#define MPCCI_QID_UGV_04                   0x2c100b03   /* Global vector 04 */
#define MPCCI_QID_UGV_05                   0x2d100b03   /* Global vector 05 */
#define MPCCI_QID_UGV_06                   0x2e100b03   /* Global vector 06 */
#define MPCCI_QID_UGV_07                   0x2f100b03   /* Global vector 07 */
#define MPCCI_QID_WALLTEMPERATURE          0x018562c1   /* Boundary temperature */
#define MPCCI_QID_WALLHEATFLUX             0x028664c1   /* Boundary normal heat flux density */
#define MPCCI_QID_WALLHTCOEFF              0x038562c1   /* Boundary heat transfer coefficient */
#define MPCCI_QID_FILMTEMPERATURE          0x048562c1   /* Film temperature */
#define MPCCI_QID_HEATRATE                 0x0585f3c1   /* Heat rate */
#define MPCCI_QID_WALLFORCE                0x0a9563c3   /* Boundary absolute force vector */
#define MPCCI_QID_RELWALLFORCE             0x0b9563c3   /* Boundary relative force vector */
#define MPCCI_QID_TEMPERATURE              0x1585f2c1   /* Temperature */
#define MPCCI_QID_TOTALTEMP                0x1685e2c1   /* Total Temperature */
#define MPCCI_QID_ABSPRESSURE              0x1f85e4c1   /* Absolute pressure */
#define MPCCI_QID_OVERPRESSURE             0x2085e4c1   /* Relative pressure */
#define MPCCI_QID_POREPRESSURE             0x2185e4c1   /* Pore pressure */
#define MPCCI_QID_ACSTPRESSURE             0x2285e4c1   /* Acoustic pressure */
#define MPCCI_QID_TOTALPRESSURE            0x2385f4c1   /* Total pressure */
#define MPCCI_QID_DYNPRESSURE              0x2485e4c1   /* Dynamic pressure */
#define MPCCI_QID_STATICPRESSURE           0x2585f4c1   /* Static pressure */
#define MPCCI_QID_VELOCITY                 0x2995f2c3   /* Velocity vector */
#define MPCCI_QID_ACCELERATION             0x2a95f2c3   /* Acceleration vector */
#define MPCCI_QID_VELOCITYMAG              0x2b85f2c1   /* Velocity magnitude */
#define MPCCI_QID_ANGULARCOORD             0x2c971143   /* Angular coordinate */
#define MPCCI_QID_ANGULARVELOCITY          0x2d951243   /* Angular velocity */
#define MPCCI_QID_ANGULARACCELERATION      0x2e951243   /* Angular acceleration */
#define MPCCI_QID_POREFLOW                 0x3106c4c1   /* Pore fluid flow */
#define MPCCI_QID_BODYFORCE                0x339284c3   /* General body force density vector */
#define MPCCI_QID_LORENTZFORCE             0x349284c3   /* Lorentz force density vector */
#define MPCCI_QID_FORCE                    0x3592f3c3   /* Force */
#define MPCCI_QID_TORQUE                   0x3692f3c3   /* Torque */
#define MPCCI_QID_HEATFLUX                 0x3d1684c3   /* Heat flux density vector */
#define MPCCI_QID_HEATSOURCE               0x3e0384c1   /* General heat source density */
#define MPCCI_QID_ENTHALPY                 0x3f0084c1   /* Enthalpy density */
#define MPCCI_QID_JOULEHEAT                0x400384c1   /* Joule heat density */
#define MPCCI_QID_JOULEHEATLIN             0x410384c1   /* Joule heat linearization */
#define MPCCI_QID_VOLFLOWVECT              0x471183c3   /* Volume flow vector */
#define MPCCI_QID_VOLFLOWRATE              0x480153c1   /* Volume flow rate */
#define MPCCI_QID_MASSFLOWVECT             0x491183c3   /* Mass flow vector */
#define MPCCI_QID_MASSFLOWRATE             0x4a0153c1   /* Mass flow rate */
#define MPCCI_QID_MASSFLUXVECT             0x4b1184c3   /* Mass flux vector */
#define MPCCI_QID_MASSFLUXRATE             0x4c0154c1   /* Mass flux rate */
#define MPCCI_QID_NPOSITION                0x5297e143   /* Nodal position */
#define MPCCI_QID_POINTPOSITION            0x54971143   /* Point position */
#define MPCCI_QID_CURRENTDENSITY           0x5b16e4c3   /* Electric current density vector */
#define MPCCI_QID_MAGNETICFLUX             0x5c16e4c3   /* Magnetic flux density vector */
#define MPCCI_QID_MAGNETICFIELD            0x5d16e2c3   /* Magnetic field vector */
#define MPCCI_QID_ELECTRICFLUX             0x5e16e4c3   /* Electric flux vector */
#define MPCCI_QID_ELECTRICFIELD            0x5f16e2c3   /* Electric field vector */
#define MPCCI_QID_CHARGEDENSITY            0x6000e2c1   /* Charge density */
#define MPCCI_QID_ELECTRICPOT              0x6105e2c1   /* Electric Potential */
#define MPCCI_QID_DENSITY                  0x650482c1   /* Density */
#define MPCCI_QID_SPECHEAT                 0x660482c1   /* Specific heat */
#define MPCCI_QID_THERMCONDX               0x6f04e2c1   /* Thermal conductivity - x */
#define MPCCI_QID_THERMCONDY               0x7004e2c1   /* Thermal conductivity - y */
#define MPCCI_QID_THERMCONDZ               0x7104e2c1   /* Thermal conductivity - z */
#define MPCCI_QID_THERMCOND1               0x7204e2c1   /* Thermal conductivity - xyz */
#define MPCCI_QID_THERMCOND3               0x7314e2c3   /* Thermal conductivity -(x,y,z) */
#define MPCCI_QID_ELECTCONDX               0x7904e2c1   /* Electric conductivity - x */
#define MPCCI_QID_ELECTCONDY               0x7a04e2c1   /* Electric conductivity - y */
#define MPCCI_QID_ELECTCONDZ               0x7b04e2c1   /* Electric conductivity - z */
#define MPCCI_QID_ELECTCOND1               0x7c04e2c1   /* Electric conductivity - xyz */
#define MPCCI_QID_ELECTCOND3               0x7d14e2c3   /* Electric conductivity - (x,y,z) */
#define MPCCI_QID_ELECTRESVX               0x8304e2c1   /* Electric resistivity - x */
#define MPCCI_QID_ELECTRESVY               0x8404e2c1   /* Electric resistivity - y */
#define MPCCI_QID_ELECTRESVZ               0x8504e2c1   /* Electric resistivity - z */
#define MPCCI_QID_ELECTRESV1               0x8604e2c1   /* Electric resistivity - xyz */
#define MPCCI_QID_ELECTRESV3               0x8714e2c3   /* Electric resistivity - (x,y,z) */
#define MPCCI_QID_YI_00                    0x0a08e2c1   /* Species mass fraction 00 */
#define MPCCI_QID_YI_01                    0x0b08e2c1   /* Species mass fraction 01 */
#define MPCCI_QID_YI_02                    0x0c08e2c1   /* Species mass fraction 02 */
#define MPCCI_QID_YI_03                    0x0d08e2c1   /* Species mass fraction 03 */
#define MPCCI_QID_YI_04                    0x0e08e2c1   /* Species mass fraction 04 */
#define MPCCI_QID_YI_05                    0x0f08e2c1   /* Species mass fraction 05 */
#define MPCCI_QID_YI_06                    0x1008e2c1   /* Species mass fraction 06 */
#define MPCCI_QID_YI_07                    0x1108e2c1   /* Species mass fraction 07 */
#define MPCCI_QID_YI_08                    0x1208e2c1   /* Species mass fraction 08 */
#define MPCCI_QID_YI_09                    0x1308e2c1   /* Species mass fraction 09 */

/* additional defines */
#define MPCCI_QID_UGS_COUNT   8
#define MPCCI_QID_UGV_COUNT   8
#define MPCCI_QID_YI_COUNT  10


#define MPCCI_QNM_INT_SWITCH               "IntFlag"
#define MPCCI_QNM_REAL_SWITCH              "RealFlag"
#define MPCCI_QNM_PHYSICAL_TIME            "PhysicalTime"
#define MPCCI_QNM_TIMESTEP_SIZE            "DeltaTime"
#define MPCCI_QNM_TIMESTEP_COUNT           "TimeStepNo"
#define MPCCI_QNM_ITERATION_COUNT          "IterationNo"
#define MPCCI_QNM_GLOBAL_RESIDUAL          "Residual"
#define MPCCI_QNM_REF_PRESSURE             "RefPressure"
#define MPCCI_QNM_VOLTAGE1                 "Voltage1"
#define MPCCI_QNM_VOLTAGE2                 "Voltage2"
#define MPCCI_QNM_VOLTAGE3                 "Voltage3"
#define MPCCI_QNM_VOLTAGE4                 "Voltage4"
#define MPCCI_QNM_CURRENT1                 "Current1"
#define MPCCI_QNM_CURRENT2                 "Current2"
#define MPCCI_QNM_CURRENT3                 "Current3"
#define MPCCI_QNM_CURRENT4                 "Current4"
#define MPCCI_QNM_MO_POSITION              "CGPosition"
#define MPCCI_QNM_MO_ANGLE                 "CGAngle"
#define MPCCI_QNM_MO_VELOCITY              "CGVelocity"
#define MPCCI_QNM_MO_OMEGA                 "CGOmega"
#define MPCCI_QNM_UGS_00                   "gs00"
#define MPCCI_QNM_UGS_01                   "gs01"
#define MPCCI_QNM_UGS_02                   "gs02"
#define MPCCI_QNM_UGS_03                   "gs03"
#define MPCCI_QNM_UGS_04                   "gs04"
#define MPCCI_QNM_UGS_05                   "gs05"
#define MPCCI_QNM_UGS_06                   "gs06"
#define MPCCI_QNM_UGS_07                   "gs07"
#define MPCCI_QNM_UGV_00                   "gv00"
#define MPCCI_QNM_UGV_01                   "gv01"
#define MPCCI_QNM_UGV_02                   "gv02"
#define MPCCI_QNM_UGV_03                   "gv03"
#define MPCCI_QNM_UGV_04                   "gv04"
#define MPCCI_QNM_UGV_05                   "gv05"
#define MPCCI_QNM_UGV_06                   "gv06"
#define MPCCI_QNM_UGV_07                   "gv07"
#define MPCCI_QNM_WALLTEMPERATURE          "WallTemp"
#define MPCCI_QNM_WALLHEATFLUX             "WallHeatFlux"
#define MPCCI_QNM_WALLHTCOEFF              "WallHTCoeff"
#define MPCCI_QNM_FILMTEMPERATURE          "FilmTemp"
#define MPCCI_QNM_HEATRATE                 "HeatRate"
#define MPCCI_QNM_WALLFORCE                "WallForce"
#define MPCCI_QNM_RELWALLFORCE             "RelWallForce"
#define MPCCI_QNM_TEMPERATURE              "Temperature"
#define MPCCI_QNM_TOTALTEMP                "TotalTemp"
#define MPCCI_QNM_ABSPRESSURE              "AbsPressure"
#define MPCCI_QNM_OVERPRESSURE             "OverPressure"
#define MPCCI_QNM_POREPRESSURE             "PorePressure"
#define MPCCI_QNM_ACSTPRESSURE             "AcstPressure"
#define MPCCI_QNM_TOTALPRESSURE            "TotalPressure"
#define MPCCI_QNM_DYNPRESSURE              "DynPressure"
#define MPCCI_QNM_STATICPRESSURE           "StaticPressure"
#define MPCCI_QNM_VELOCITY                 "Velocity"
#define MPCCI_QNM_ACCELERATION             "Acceleration"
#define MPCCI_QNM_VELOCITYMAG              "VelocityMagnitude"
#define MPCCI_QNM_ANGULARCOORD             "AngularCoordinate"
#define MPCCI_QNM_ANGULARVELOCITY          "AngularVelocity"
#define MPCCI_QNM_ANGULARACCELERATION      "AngularAcceleration"
#define MPCCI_QNM_POREFLOW                 "PorousFlow"
#define MPCCI_QNM_BODYFORCE                "BodyForce"
#define MPCCI_QNM_LORENTZFORCE             "LorentzForce"
#define MPCCI_QNM_FORCE                    "Force"
#define MPCCI_QNM_TORQUE                   "Torque"
#define MPCCI_QNM_HEATFLUX                 "HeatFlux"
#define MPCCI_QNM_HEATSOURCE               "HeatSource"
#define MPCCI_QNM_ENTHALPY                 "Enthalpy"
#define MPCCI_QNM_JOULEHEAT                "JouleHeat"
#define MPCCI_QNM_JOULEHEATLIN             "JouleHeatLin"
#define MPCCI_QNM_VOLFLOWVECT              "VolumeFlow"
#define MPCCI_QNM_VOLFLOWRATE              "VolumeFlowRate"
#define MPCCI_QNM_MASSFLOWVECT             "MassFlowVect"
#define MPCCI_QNM_MASSFLOWRATE             "MassFlowRate"
#define MPCCI_QNM_MASSFLUXVECT             "MassFluxVect"
#define MPCCI_QNM_MASSFLUXRATE             "MassFluxRate"
#define MPCCI_QNM_NPOSITION                "NPosition"
#define MPCCI_QNM_POINTPOSITION            "PointPosition"
#define MPCCI_QNM_CURRENTDENSITY           "CurrentDensity"
#define MPCCI_QNM_MAGNETICFLUX             "MagneticFlux"
#define MPCCI_QNM_MAGNETICFIELD            "MagneticField"
#define MPCCI_QNM_ELECTRICFLUX             "ElectricFlux"
#define MPCCI_QNM_ELECTRICFIELD            "ElectricField"
#define MPCCI_QNM_CHARGEDENSITY            "ChargeDensity"
#define MPCCI_QNM_ELECTRICPOT              "ElectricPot"
#define MPCCI_QNM_DENSITY                  "Density"
#define MPCCI_QNM_SPECHEAT                 "SpecificHeat"
#define MPCCI_QNM_THERMCONDX               "ThermCondX"
#define MPCCI_QNM_THERMCONDY               "ThermCondY"
#define MPCCI_QNM_THERMCONDZ               "ThermCondZ"
#define MPCCI_QNM_THERMCOND1               "ThermCond1"
#define MPCCI_QNM_THERMCOND3               "ThermCond3"
#define MPCCI_QNM_ELECTCONDX               "ElectrCondX"
#define MPCCI_QNM_ELECTCONDY               "ElectrCondY"
#define MPCCI_QNM_ELECTCONDZ               "ElectrCondZ"
#define MPCCI_QNM_ELECTCOND1               "ElectrCond1"
#define MPCCI_QNM_ELECTCOND3               "ElectrCond3"
#define MPCCI_QNM_ELECTRESVX               "ElectrResX"
#define MPCCI_QNM_ELECTRESVY               "ElectrResY"
#define MPCCI_QNM_ELECTRESVZ               "ElectrResZ"
#define MPCCI_QNM_ELECTRESV1               "ElectrRes1"
#define MPCCI_QNM_ELECTRESV3               "ElectrRes3"
#define MPCCI_QNM_YI_00                    "YI00"
#define MPCCI_QNM_YI_01                    "YI01"
#define MPCCI_QNM_YI_02                    "YI02"
#define MPCCI_QNM_YI_03                    "YI03"
#define MPCCI_QNM_YI_04                    "YI04"
#define MPCCI_QNM_YI_05                    "YI05"
#define MPCCI_QNM_YI_06                    "YI06"
#define MPCCI_QNM_YI_07                    "YI07"
#define MPCCI_QNM_YI_08                    "YI08"
#define MPCCI_QNM_YI_09                    "YI09"
