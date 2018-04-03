#ifndef __MPCCI4_ELEMENTS_HEADER_INCLUDED
#define __MPCCI4_ELEMENTS_HEADER_INCLUDED
/*****************************************************************************************
 * mpcci_elements.h
 *
 * MpCCI header file with basic element type definitions
 *
 * Element type (ETYP) definitions: The ETYP consists of the ...
 *       - element dimension (0=point, 1=line, 2=face, 3=volume) (DIM)
 *       - element unique signature for a unique shape (0-f)     (SIG)
 *       - element no. of geometrical vertices <= nnodes         (NVX)
 *       - element basic shape                                   (SHP)
 *       - element no. of nodes in element/cell                  (NND)
 *    and is specified within a single unsigned.
 *
 *    nibbles: |0000|0000|0000|0000|0000|0000|0000|0000|
 *             |-DIM|-SIG|-NVX|----------NND-----------|
 *             |------SHP-----|----------NND-----------|
 *
 * $Id: mpcci_elements.h 104 2013-08-16 11:56:59Z pbayrasy $
 *
 ****************************************************************************************/

#include "mpcci_defs.h"

/*
 * This is just an upper limit for the no. of suppported different element types
 * used as a limit for array dimensions
 */
#define MPCCI_ETYPES_MAX               64

#define MPCCI_ETYP(_shp,_nne)        ( (_shp) | (_nne) )

#define MPCCI_ETYP_DIM(_etyp)          ((MPCCI_UCAST(_etyp)) >> 28)          /* get element dimension 0..3 */
#define MPCCI_ETYP_SHP(_etyp)          ((MPCCI_UCAST(_etyp)) & 0xfff00000)   /* mask the full shape */
#define MPCCI_ETYP_SIG(_etyp)          ((MPCCI_UCAST(_etyp)) & 0x0f000000)   /* mask the signature only */
#define MPCCI_ETYP_NNE(_etyp)          ((MPCCI_UCAST(_etyp)) & 0x000fffff)   /* mask the no. of nodes */
#define MPCCI_ETYP_NVX_BASE(_etyp)   ( ((MPCCI_UCAST(_etyp)) >> 20) & 0x0f ) /* get no. of vertices == minimal no. of nodes */


/* point elements with and without coordinates */
#define MPCCI_ESHP_IPNT                0x00000000 /* integration point of a network code (no coordinates required) */
#define MPCCI_ESHP_CPNT                0x00100000 /* point with coordinates */

/* beam elements */
#define MPCCI_ESHP_LINE                0x10200000 /* 102: line with at least 2 nodes */

/* face type elements */
#define MPCCI_ESHP_TRIA                0x20300000 /* 203: face with at least 3 nodes */
#define MPCCI_ESHP_QUAD                0x20400000 /* 204: face with at least 4 nodes */
#define MPCCI_ESHP_MEMB                0x21400000 /* 214: sig=1: FEA thin membrane */
#define MPCCI_ESHP_SHELL               0x22400000 /* 224: sig=2: FEA thick shell element */
#define MPCCI_ESHP_POLY                0x2f500000 /* 2f5: sig=0xf: polygonal face with > 4 vertices */

/* volume elements */
#define MPCCI_ESHP_TET                 0x30400000 /* 304: volume cell with at least 4 vertices */
#define MPCCI_ESHP_PYRAM               0x30500000 /* 305: volume cell with at least 5 vertices */
#define MPCCI_ESHP_WEDGE               0x30600000 /* 306: volume cell with at least 6 vertices */
#define MPCCI_ESHP_HEX                 0x30800000 /* 308: volume cell with at least 8 vertices */
#define MPCCI_ESHP_PCELL               0x3f900000 /* sig=0xf: CFD polyhedron > 8(?) nodes (in fact face definitions) */


/* single integration point: network code only */
#define MPCCI_ETYP_IPNT                ( MPCCI_ESHP_IPNT  |  0 )

/* point or point clouds on a line, surface or inside a volume */
#define MPCCI_ETYP_CPNT                ( MPCCI_ESHP_CPNT  |  1 )

/* some beam elements */
#define MPCCI_ETYP_LINE2               ( MPCCI_ESHP_LINE  |  2 )
#define MPCCI_ETYP_LINE3               ( MPCCI_ESHP_LINE  |  3 )
#define MPCCI_ETYP_LINE4               ( MPCCI_ESHP_LINE  |  4 )

/* supported TRIA elements */
#define MPCCI_ETYP_TRIA3               ( MPCCI_ESHP_TRIA  |  3 )
#define MPCCI_ETYP_TRIA4               ( MPCCI_ESHP_TRIA  |  4 )
#define MPCCI_ETYP_TRIA6               ( MPCCI_ESHP_TRIA  |  6 )
#define MPCCI_ETYP_TRIA7               ( MPCCI_ESHP_TRIA  |  7 )
#define MPCCI_ETYP_TRIA10              ( MPCCI_ESHP_TRIA  | 10 )

/* supported QUAD elements */
#define MPCCI_ETYP_QUAD4               ( MPCCI_ESHP_QUAD  |  4 )
#define MPCCI_ETYP_QUAD5               ( MPCCI_ESHP_QUAD  |  5 )
#define MPCCI_ETYP_QUAD8               ( MPCCI_ESHP_QUAD  |  8 )
#define MPCCI_ETYP_QUAD9               ( MPCCI_ESHP_QUAD  |  9 )
#define MPCCI_ETYP_QUAD12              ( MPCCI_ESHP_QUAD  | 12 )
#define MPCCI_ETYP_QUAD16              ( MPCCI_ESHP_QUAD  | 16 )

#define MPCCI_ETYP_MEMB4               ( MPCCI_ESHP_MEMB  |  4 )
#define MPCCI_ETYP_MEMB5               ( MPCCI_ESHP_MEMB  |  5 )
#define MPCCI_ETYP_MEMB8               ( MPCCI_ESHP_MEMB  |  8 )
#define MPCCI_ETYP_MEMB9               ( MPCCI_ESHP_MEMB  |  9 )

#define MPCCI_ETYP_SHELL4              ( MPCCI_ESHP_SHELL |  4 )
#define MPCCI_ETYP_SHELL5              ( MPCCI_ESHP_SHELL |  5 )
#define MPCCI_ETYP_SHELL8              ( MPCCI_ESHP_SHELL |  8 )
#define MPCCI_ETYP_SHELL9              ( MPCCI_ESHP_SHELL |  9 )

/* supported TET elements */
#define MPCCI_ETYP_TET4                ( MPCCI_ESHP_TET   |  4 )
#define MPCCI_ETYP_TET5                ( MPCCI_ESHP_TET   |  5 )
#define MPCCI_ETYP_TET10               ( MPCCI_ESHP_TET   | 10 )
#define MPCCI_ETYP_TET11               ( MPCCI_ESHP_TET   | 11 )

/* supported PYRAM elements */
#define MPCCI_ETYP_PYRAM5              ( MPCCI_ESHP_PYRAM |  5 )
#define MPCCI_ETYP_PYRAM6              ( MPCCI_ESHP_PYRAM |  6 )
#define MPCCI_ETYP_PYRAM13             ( MPCCI_ESHP_PYRAM | 13 )
#define MPCCI_ETYP_PYRAM14             ( MPCCI_ESHP_PYRAM | 14 )

/* supported WEDGE elements */
#define MPCCI_ETYP_WEDGE6              ( MPCCI_ESHP_WEDGE |  6 )
#define MPCCI_ETYP_WEDGE7              ( MPCCI_ESHP_WEDGE |  7 )
#define MPCCI_ETYP_WEDGE15             ( MPCCI_ESHP_WEDGE | 15 )
#define MPCCI_ETYP_WEDGE16             ( MPCCI_ESHP_WEDGE | 16 )

/* supported HEX elements */
#define MPCCI_ETYP_HEX8                ( MPCCI_ESHP_HEX   |  8 )
#define MPCCI_ETYP_HEX9                ( MPCCI_ESHP_HEX   |  9 )
#define MPCCI_ETYP_HEX20               ( MPCCI_ESHP_HEX   | 20 )
#define MPCCI_ETYP_HEX21               ( MPCCI_ESHP_HEX   | 21 )
#define MPCCI_ETYP_HEX27               ( MPCCI_ESHP_HEX   | 27 )

/* polyhedrons can never be predefined. Some examples ... */
#define MPCCI_ETYP_PCELL21             ( MPCCI_ESHP_PCELL | 21 )
#define MPCCI_ETYP_PCELL99             ( MPCCI_ESHP_PCELL | 99 )

/* polygonal faces NNE>=5 and NVX==5 */
#define MPCCI_ETYP_POLY3               ( MPCCI_ESHP_POLY  |   3 ) /* in fact TRIA3 */
#define MPCCI_ETYP_POLY4               ( MPCCI_ESHP_POLY  |   4 ) /* in fact QUAD4 */
#define MPCCI_ETYP_POLY5               ( MPCCI_ESHP_POLY  |   5 )
#define MPCCI_ETYP_POLY6               ( MPCCI_ESHP_POLY  |   6 )
#define MPCCI_ETYP_POLY7               ( MPCCI_ESHP_POLY  |   7 )
#define MPCCI_ETYP_POLY8               ( MPCCI_ESHP_POLY  |   8 )
#define MPCCI_ETYP_POLY9               ( MPCCI_ESHP_POLY  |   9 )
#define MPCCI_ETYP_POLY10              ( MPCCI_ESHP_POLY  |  10 )
#define MPCCI_ETYP_POLY11              ( MPCCI_ESHP_POLY  |  11 )
#define MPCCI_ETYP_POLY12              ( MPCCI_ESHP_POLY  |  12 )
#define MPCCI_ETYP_POLY13              ( MPCCI_ESHP_POLY  |  13 )
#define MPCCI_ETYP_POLY14              ( MPCCI_ESHP_POLY  |  14 )
#define MPCCI_ETYP_POLY15              ( MPCCI_ESHP_POLY  |  15 )
#define MPCCI_ETYP_POLY16              ( MPCCI_ESHP_POLY  |  16 )
#define MPCCI_ETYP_POLY17              ( MPCCI_ESHP_POLY  |  17 )
#define MPCCI_ETYP_POLY18              ( MPCCI_ESHP_POLY  |  18 )
#define MPCCI_ETYP_POLY19              ( MPCCI_ESHP_POLY  |  19 )
#define MPCCI_ETYP_POLY20              ( MPCCI_ESHP_POLY  |  20 )
#define MPCCI_ETYP_POLY21              ( MPCCI_ESHP_POLY  |  21 )
#define MPCCI_ETYP_POLY22              ( MPCCI_ESHP_POLY  |  22 )
#define MPCCI_ETYP_POLY23              ( MPCCI_ESHP_POLY  |  23 )
#define MPCCI_ETYP_POLY24              ( MPCCI_ESHP_POLY  |  24 )
#define MPCCI_ETYP_POLY25              ( MPCCI_ESHP_POLY  |  25 )
#define MPCCI_ETYP_POLY26              ( MPCCI_ESHP_POLY  |  26 )
#define MPCCI_ETYP_POLY27              ( MPCCI_ESHP_POLY  |  27 )
#define MPCCI_ETYP_POLY28              ( MPCCI_ESHP_POLY  |  28 )
#define MPCCI_ETYP_POLY29              ( MPCCI_ESHP_POLY  |  29 )
#define MPCCI_ETYP_POLY30              ( MPCCI_ESHP_POLY  |  30 )
#define MPCCI_ETYP_POLY31              ( MPCCI_ESHP_POLY  |  31 )
#define MPCCI_ETYP_POLY32              ( MPCCI_ESHP_POLY  |  32 )
#define MPCCI_ETYP_POLY33              ( MPCCI_ESHP_POLY  |  33 )
#define MPCCI_ETYP_POLY34              ( MPCCI_ESHP_POLY  |  34 )
#define MPCCI_ETYP_POLY35              ( MPCCI_ESHP_POLY  |  35 )
#define MPCCI_ETYP_POLY36              ( MPCCI_ESHP_POLY  |  36 )
#define MPCCI_ETYP_POLY37              ( MPCCI_ESHP_POLY  |  37 )
#define MPCCI_ETYP_POLY38              ( MPCCI_ESHP_POLY  |  38 )
#define MPCCI_ETYP_POLY39              ( MPCCI_ESHP_POLY  |  39 )
#define MPCCI_ETYP_POLY40              ( MPCCI_ESHP_POLY  |  40 )
#define MPCCI_ETYP_POLY41              ( MPCCI_ESHP_POLY  |  41 )
#define MPCCI_ETYP_POLY42              ( MPCCI_ESHP_POLY  |  42 )
#define MPCCI_ETYP_POLY43              ( MPCCI_ESHP_POLY  |  43 )
#define MPCCI_ETYP_POLY44              ( MPCCI_ESHP_POLY  |  44 )
#define MPCCI_ETYP_POLY45              ( MPCCI_ESHP_POLY  |  45 )
#define MPCCI_ETYP_POLY46              ( MPCCI_ESHP_POLY  |  46 )
#define MPCCI_ETYP_POLY47              ( MPCCI_ESHP_POLY  |  47 )
#define MPCCI_ETYP_POLY48              ( MPCCI_ESHP_POLY  |  48 )
#define MPCCI_ETYP_POLY49              ( MPCCI_ESHP_POLY  |  49 )
#define MPCCI_ETYP_POLY50              ( MPCCI_ESHP_POLY  |  50 )
#define MPCCI_ETYP_POLY51              ( MPCCI_ESHP_POLY  |  51 )
#define MPCCI_ETYP_POLY52              ( MPCCI_ESHP_POLY  |  52 )
#define MPCCI_ETYP_POLY53              ( MPCCI_ESHP_POLY  |  53 )
#define MPCCI_ETYP_POLY54              ( MPCCI_ESHP_POLY  |  54 )
#define MPCCI_ETYP_POLY55              ( MPCCI_ESHP_POLY  |  55 )
#define MPCCI_ETYP_POLY56              ( MPCCI_ESHP_POLY  |  56 )
#define MPCCI_ETYP_POLY57              ( MPCCI_ESHP_POLY  |  57 )
#define MPCCI_ETYP_POLY58              ( MPCCI_ESHP_POLY  |  58 )
#define MPCCI_ETYP_POLY59              ( MPCCI_ESHP_POLY  |  59 )
#define MPCCI_ETYP_POLY60              ( MPCCI_ESHP_POLY  |  60 )
#define MPCCI_ETYP_POLY61              ( MPCCI_ESHP_POLY  |  61 )
#define MPCCI_ETYP_POLY62              ( MPCCI_ESHP_POLY  |  62 )
#define MPCCI_ETYP_POLY63              ( MPCCI_ESHP_POLY  |  63 )
#define MPCCI_ETYP_POLY64              ( MPCCI_ESHP_POLY  |  64 )
#define MPCCI_ETYP_POLY65              ( MPCCI_ESHP_POLY  |  65 )
#define MPCCI_ETYP_POLY66              ( MPCCI_ESHP_POLY  |  66 )
#define MPCCI_ETYP_POLY67              ( MPCCI_ESHP_POLY  |  67 )
#define MPCCI_ETYP_POLY68              ( MPCCI_ESHP_POLY  |  68 )
#define MPCCI_ETYP_POLY69              ( MPCCI_ESHP_POLY  |  69 )
#define MPCCI_ETYP_POLY70              ( MPCCI_ESHP_POLY  |  70 )
#define MPCCI_ETYP_POLY71              ( MPCCI_ESHP_POLY  |  71 )
#define MPCCI_ETYP_POLY72              ( MPCCI_ESHP_POLY  |  72 )
#define MPCCI_ETYP_POLY73              ( MPCCI_ESHP_POLY  |  73 )
#define MPCCI_ETYP_POLY74              ( MPCCI_ESHP_POLY  |  74 )
#define MPCCI_ETYP_POLY75              ( MPCCI_ESHP_POLY  |  75 )
#define MPCCI_ETYP_POLY76              ( MPCCI_ESHP_POLY  |  76 )
#define MPCCI_ETYP_POLY77              ( MPCCI_ESHP_POLY  |  77 )
#define MPCCI_ETYP_POLY78              ( MPCCI_ESHP_POLY  |  78 )
#define MPCCI_ETYP_POLY79              ( MPCCI_ESHP_POLY  |  79 )
#define MPCCI_ETYP_POLY80              ( MPCCI_ESHP_POLY  |  80 )
#define MPCCI_ETYP_POLY81              ( MPCCI_ESHP_POLY  |  81 )
#define MPCCI_ETYP_POLY82              ( MPCCI_ESHP_POLY  |  82 )
#define MPCCI_ETYP_POLY83              ( MPCCI_ESHP_POLY  |  83 )
#define MPCCI_ETYP_POLY84              ( MPCCI_ESHP_POLY  |  84 )
#define MPCCI_ETYP_POLY85              ( MPCCI_ESHP_POLY  |  85 )
#define MPCCI_ETYP_POLY86              ( MPCCI_ESHP_POLY  |  86 )
#define MPCCI_ETYP_POLY87              ( MPCCI_ESHP_POLY  |  87 )
#define MPCCI_ETYP_POLY88              ( MPCCI_ESHP_POLY  |  88 )
#define MPCCI_ETYP_POLY89              ( MPCCI_ESHP_POLY  |  89 )
#define MPCCI_ETYP_POLY90              ( MPCCI_ESHP_POLY  |  90 )
#define MPCCI_ETYP_POLY91              ( MPCCI_ESHP_POLY  |  91 )
#define MPCCI_ETYP_POLY92              ( MPCCI_ESHP_POLY  |  92 )
#define MPCCI_ETYP_POLY93              ( MPCCI_ESHP_POLY  |  93 )
#define MPCCI_ETYP_POLY94              ( MPCCI_ESHP_POLY  |  94 )
#define MPCCI_ETYP_POLY95              ( MPCCI_ESHP_POLY  |  95 )
#define MPCCI_ETYP_POLY96              ( MPCCI_ESHP_POLY  |  96 )
#define MPCCI_ETYP_POLY97              ( MPCCI_ESHP_POLY  |  97 )
#define MPCCI_ETYP_POLY98              ( MPCCI_ESHP_POLY  |  98 )
#define MPCCI_ETYP_POLY99              ( MPCCI_ESHP_POLY  |  99 )
#define MPCCI_ETYP_POLY100             ( MPCCI_ESHP_POLY  | 100 )
#define MPCCI_ETYP_POLY101             ( MPCCI_ESHP_POLY  | 101 )
#define MPCCI_ETYP_POLY102             ( MPCCI_ESHP_POLY  | 102 )
#define MPCCI_ETYP_POLY103             ( MPCCI_ESHP_POLY  | 103 )
#define MPCCI_ETYP_POLY104             ( MPCCI_ESHP_POLY  | 104 )
#define MPCCI_ETYP_POLY105             ( MPCCI_ESHP_POLY  | 105 )
#define MPCCI_ETYP_POLY106             ( MPCCI_ESHP_POLY  | 106 )
#define MPCCI_ETYP_POLY107             ( MPCCI_ESHP_POLY  | 107 )
#define MPCCI_ETYP_POLY108             ( MPCCI_ESHP_POLY  | 108 )
#define MPCCI_ETYP_POLY109             ( MPCCI_ESHP_POLY  | 109 )
#define MPCCI_ETYP_POLY110             ( MPCCI_ESHP_POLY  | 110 )
#define MPCCI_ETYP_POLY111             ( MPCCI_ESHP_POLY  | 111 )
#define MPCCI_ETYP_POLY112             ( MPCCI_ESHP_POLY  | 112 )
#define MPCCI_ETYP_POLY113             ( MPCCI_ESHP_POLY  | 113 )
#define MPCCI_ETYP_POLY114             ( MPCCI_ESHP_POLY  | 114 )
#define MPCCI_ETYP_POLY115             ( MPCCI_ESHP_POLY  | 115 )
#define MPCCI_ETYP_POLY116             ( MPCCI_ESHP_POLY  | 116 )
#define MPCCI_ETYP_POLY117             ( MPCCI_ESHP_POLY  | 117 )
#define MPCCI_ETYP_POLY118             ( MPCCI_ESHP_POLY  | 118 )
#define MPCCI_ETYP_POLY119             ( MPCCI_ESHP_POLY  | 119 )
#define MPCCI_ETYP_POLY120             ( MPCCI_ESHP_POLY  | 120 )
#define MPCCI_ETYP_POLY121             ( MPCCI_ESHP_POLY  | 121 )
#define MPCCI_ETYP_POLY122             ( MPCCI_ESHP_POLY  | 122 )
#define MPCCI_ETYP_POLY123             ( MPCCI_ESHP_POLY  | 123 )
#define MPCCI_ETYP_POLY124             ( MPCCI_ESHP_POLY  | 124 )
#define MPCCI_ETYP_POLY125             ( MPCCI_ESHP_POLY  | 125 )
#define MPCCI_ETYP_POLY126             ( MPCCI_ESHP_POLY  | 126 )
#define MPCCI_ETYP_POLY127             ( MPCCI_ESHP_POLY  | 127 )
#define MPCCI_ETYP_POLY128             ( MPCCI_ESHP_POLY  | 128 )
#define MPCCI_ETYP_POLY129             ( MPCCI_ESHP_POLY  | 129 )
#define MPCCI_ETYP_POLY130             ( MPCCI_ESHP_POLY  | 130 )
#define MPCCI_ETYP_POLY131             ( MPCCI_ESHP_POLY  | 131 )
#define MPCCI_ETYP_POLY132             ( MPCCI_ESHP_POLY  | 132 )
#define MPCCI_ETYP_POLY133             ( MPCCI_ESHP_POLY  | 133 )
#define MPCCI_ETYP_POLY134             ( MPCCI_ESHP_POLY  | 134 )
#define MPCCI_ETYP_POLY135             ( MPCCI_ESHP_POLY  | 135 )
#define MPCCI_ETYP_POLY136             ( MPCCI_ESHP_POLY  | 136 )
#define MPCCI_ETYP_POLY137             ( MPCCI_ESHP_POLY  | 137 )
#define MPCCI_ETYP_POLY138             ( MPCCI_ESHP_POLY  | 138 )
#define MPCCI_ETYP_POLY139             ( MPCCI_ESHP_POLY  | 139 )
#define MPCCI_ETYP_POLY140             ( MPCCI_ESHP_POLY  | 140 )
#define MPCCI_ETYP_POLY141             ( MPCCI_ESHP_POLY  | 141 )
#define MPCCI_ETYP_POLY142             ( MPCCI_ESHP_POLY  | 142 )
#define MPCCI_ETYP_POLY143             ( MPCCI_ESHP_POLY  | 143 )
#define MPCCI_ETYP_POLY144             ( MPCCI_ESHP_POLY  | 144 )
#define MPCCI_ETYP_POLY145             ( MPCCI_ESHP_POLY  | 145 )
#define MPCCI_ETYP_POLY146             ( MPCCI_ESHP_POLY  | 146 )
#define MPCCI_ETYP_POLY147             ( MPCCI_ESHP_POLY  | 147 )
#define MPCCI_ETYP_POLY148             ( MPCCI_ESHP_POLY  | 148 )
#define MPCCI_ETYP_POLY149             ( MPCCI_ESHP_POLY  | 149 )
#define MPCCI_ETYP_POLY150             ( MPCCI_ESHP_POLY  | 150 )
#define MPCCI_ETYP_POLY151             ( MPCCI_ESHP_POLY  | 151 )
#define MPCCI_ETYP_POLY152             ( MPCCI_ESHP_POLY  | 152 )
#define MPCCI_ETYP_POLY153             ( MPCCI_ESHP_POLY  | 153 )
#define MPCCI_ETYP_POLY154             ( MPCCI_ESHP_POLY  | 154 )
#define MPCCI_ETYP_POLY155             ( MPCCI_ESHP_POLY  | 155 )
#define MPCCI_ETYP_POLY156             ( MPCCI_ESHP_POLY  | 156 )
#define MPCCI_ETYP_POLY157             ( MPCCI_ESHP_POLY  | 157 )
#define MPCCI_ETYP_POLY158             ( MPCCI_ESHP_POLY  | 158 )
#define MPCCI_ETYP_POLY159             ( MPCCI_ESHP_POLY  | 159 )
#define MPCCI_ETYP_POLY160             ( MPCCI_ESHP_POLY  | 160 )
#define MPCCI_ETYP_POLY161             ( MPCCI_ESHP_POLY  | 161 )
#define MPCCI_ETYP_POLY162             ( MPCCI_ESHP_POLY  | 162 )
#define MPCCI_ETYP_POLY163             ( MPCCI_ESHP_POLY  | 163 )
#define MPCCI_ETYP_POLY164             ( MPCCI_ESHP_POLY  | 164 )
#define MPCCI_ETYP_POLY165             ( MPCCI_ESHP_POLY  | 165 )
#define MPCCI_ETYP_POLY166             ( MPCCI_ESHP_POLY  | 166 )
#define MPCCI_ETYP_POLY167             ( MPCCI_ESHP_POLY  | 167 )
#define MPCCI_ETYP_POLY168             ( MPCCI_ESHP_POLY  | 168 )
#define MPCCI_ETYP_POLY169             ( MPCCI_ESHP_POLY  | 169 )
#define MPCCI_ETYP_POLY170             ( MPCCI_ESHP_POLY  | 170 )
#define MPCCI_ETYP_POLY171             ( MPCCI_ESHP_POLY  | 171 )
#define MPCCI_ETYP_POLY172             ( MPCCI_ESHP_POLY  | 172 )
#define MPCCI_ETYP_POLY173             ( MPCCI_ESHP_POLY  | 173 )
#define MPCCI_ETYP_POLY174             ( MPCCI_ESHP_POLY  | 174 )
#define MPCCI_ETYP_POLY175             ( MPCCI_ESHP_POLY  | 175 )
#define MPCCI_ETYP_POLY176             ( MPCCI_ESHP_POLY  | 176 )
#define MPCCI_ETYP_POLY177             ( MPCCI_ESHP_POLY  | 177 )
#define MPCCI_ETYP_POLY178             ( MPCCI_ESHP_POLY  | 178 )
#define MPCCI_ETYP_POLY179             ( MPCCI_ESHP_POLY  | 179 )
#define MPCCI_ETYP_POLY180             ( MPCCI_ESHP_POLY  | 180 )
#define MPCCI_ETYP_POLY181             ( MPCCI_ESHP_POLY  | 181 )
#define MPCCI_ETYP_POLY182             ( MPCCI_ESHP_POLY  | 182 )
#define MPCCI_ETYP_POLY183             ( MPCCI_ESHP_POLY  | 183 )
#define MPCCI_ETYP_POLY184             ( MPCCI_ESHP_POLY  | 184 )
#define MPCCI_ETYP_POLY185             ( MPCCI_ESHP_POLY  | 185 )
#define MPCCI_ETYP_POLY186             ( MPCCI_ESHP_POLY  | 186 )
#define MPCCI_ETYP_POLY187             ( MPCCI_ESHP_POLY  | 187 )
#define MPCCI_ETYP_POLY188             ( MPCCI_ESHP_POLY  | 188 )
#define MPCCI_ETYP_POLY189             ( MPCCI_ESHP_POLY  | 189 )
#define MPCCI_ETYP_POLY190             ( MPCCI_ESHP_POLY  | 190 )
#define MPCCI_ETYP_POLY191             ( MPCCI_ESHP_POLY  | 191 )
#define MPCCI_ETYP_POLY192             ( MPCCI_ESHP_POLY  | 192 )
#define MPCCI_ETYP_POLY193             ( MPCCI_ESHP_POLY  | 193 )
#define MPCCI_ETYP_POLY194             ( MPCCI_ESHP_POLY  | 194 )
#define MPCCI_ETYP_POLY195             ( MPCCI_ESHP_POLY  | 195 )
#define MPCCI_ETYP_POLY196             ( MPCCI_ESHP_POLY  | 196 )
#define MPCCI_ETYP_POLY197             ( MPCCI_ESHP_POLY  | 197 )
#define MPCCI_ETYP_POLY198             ( MPCCI_ESHP_POLY  | 198 )
#define MPCCI_ETYP_POLY199             ( MPCCI_ESHP_POLY  | 199 )
#define MPCCI_ETYP_POLY200             ( MPCCI_ESHP_POLY  | 200 )
#define MPCCI_ETYP_POLY201             ( MPCCI_ESHP_POLY  | 201 )
#define MPCCI_ETYP_POLY202             ( MPCCI_ESHP_POLY  | 202 )
#define MPCCI_ETYP_POLY203             ( MPCCI_ESHP_POLY  | 203 )
#define MPCCI_ETYP_POLY204             ( MPCCI_ESHP_POLY  | 204 )
#define MPCCI_ETYP_POLY205             ( MPCCI_ESHP_POLY  | 205 )
#define MPCCI_ETYP_POLY206             ( MPCCI_ESHP_POLY  | 206 )
#define MPCCI_ETYP_POLY207             ( MPCCI_ESHP_POLY  | 207 )
#define MPCCI_ETYP_POLY208             ( MPCCI_ESHP_POLY  | 208 )
#define MPCCI_ETYP_POLY209             ( MPCCI_ESHP_POLY  | 209 )
#define MPCCI_ETYP_POLY210             ( MPCCI_ESHP_POLY  | 210 )
#define MPCCI_ETYP_POLY211             ( MPCCI_ESHP_POLY  | 211 )
#define MPCCI_ETYP_POLY212             ( MPCCI_ESHP_POLY  | 212 )
#define MPCCI_ETYP_POLY213             ( MPCCI_ESHP_POLY  | 213 )
#define MPCCI_ETYP_POLY214             ( MPCCI_ESHP_POLY  | 214 )
#define MPCCI_ETYP_POLY215             ( MPCCI_ESHP_POLY  | 215 )
#define MPCCI_ETYP_POLY216             ( MPCCI_ESHP_POLY  | 216 )
#define MPCCI_ETYP_POLY217             ( MPCCI_ESHP_POLY  | 217 )
#define MPCCI_ETYP_POLY218             ( MPCCI_ESHP_POLY  | 218 )
#define MPCCI_ETYP_POLY219             ( MPCCI_ESHP_POLY  | 219 )
#define MPCCI_ETYP_POLY220             ( MPCCI_ESHP_POLY  | 220 )
#define MPCCI_ETYP_POLY221             ( MPCCI_ESHP_POLY  | 221 )
#define MPCCI_ETYP_POLY222             ( MPCCI_ESHP_POLY  | 222 )
#define MPCCI_ETYP_POLY223             ( MPCCI_ESHP_POLY  | 223 )
#define MPCCI_ETYP_POLY224             ( MPCCI_ESHP_POLY  | 224 )
#define MPCCI_ETYP_POLY225             ( MPCCI_ESHP_POLY  | 225 )
#define MPCCI_ETYP_POLY226             ( MPCCI_ESHP_POLY  | 226 )
#define MPCCI_ETYP_POLY227             ( MPCCI_ESHP_POLY  | 227 )
#define MPCCI_ETYP_POLY228             ( MPCCI_ESHP_POLY  | 228 )
#define MPCCI_ETYP_POLY229             ( MPCCI_ESHP_POLY  | 229 )
#define MPCCI_ETYP_POLY230             ( MPCCI_ESHP_POLY  | 230 )
#define MPCCI_ETYP_POLY231             ( MPCCI_ESHP_POLY  | 231 )
#define MPCCI_ETYP_POLY232             ( MPCCI_ESHP_POLY  | 232 )
#define MPCCI_ETYP_POLY233             ( MPCCI_ESHP_POLY  | 233 )
#define MPCCI_ETYP_POLY234             ( MPCCI_ESHP_POLY  | 234 )
#define MPCCI_ETYP_POLY235             ( MPCCI_ESHP_POLY  | 235 )
#define MPCCI_ETYP_POLY236             ( MPCCI_ESHP_POLY  | 236 )
#define MPCCI_ETYP_POLY237             ( MPCCI_ESHP_POLY  | 237 )
#define MPCCI_ETYP_POLY238             ( MPCCI_ESHP_POLY  | 238 )
#define MPCCI_ETYP_POLY239             ( MPCCI_ESHP_POLY  | 239 )
#define MPCCI_ETYP_POLY240             ( MPCCI_ESHP_POLY  | 240 )
#define MPCCI_ETYP_POLY241             ( MPCCI_ESHP_POLY  | 241 )
#define MPCCI_ETYP_POLY242             ( MPCCI_ESHP_POLY  | 242 )
#define MPCCI_ETYP_POLY243             ( MPCCI_ESHP_POLY  | 243 )
#define MPCCI_ETYP_POLY244             ( MPCCI_ESHP_POLY  | 244 )
#define MPCCI_ETYP_POLY245             ( MPCCI_ESHP_POLY  | 245 )
#define MPCCI_ETYP_POLY246             ( MPCCI_ESHP_POLY  | 246 )
#define MPCCI_ETYP_POLY247             ( MPCCI_ESHP_POLY  | 247 )
#define MPCCI_ETYP_POLY248             ( MPCCI_ESHP_POLY  | 248 )
#define MPCCI_ETYP_POLY249             ( MPCCI_ESHP_POLY  | 249 )
#define MPCCI_ETYP_POLY250             ( MPCCI_ESHP_POLY  | 250 )
#define MPCCI_ETYP_POLY251             ( MPCCI_ESHP_POLY  | 251 )
#define MPCCI_ETYP_POLY252             ( MPCCI_ESHP_POLY  | 252 )
#define MPCCI_ETYP_POLY253             ( MPCCI_ESHP_POLY  | 253 )
#define MPCCI_ETYP_POLY254             ( MPCCI_ESHP_POLY  | 254 )
#define MPCCI_ETYP_POLY255             ( MPCCI_ESHP_POLY  | 255 )
#define MPCCI_ETYP_POLY256             ( MPCCI_ESHP_POLY  | 256 )

/*
 * Helper macros for the special treatment of polyhedrons:
 */
#define MPCCI_ETYP_HAS_PSIG(_etyp)     ( MPCCI_ETYP_SIG(_etyp) == 0x0f000000 )

#define MPCCI_ESHP_IS_POLY(_eshp)      (       (_eshp)         == MPCCI_ESHP_POLY )
#define MPCCI_ETYP_IS_POLY(_etyp)      ( MPCCI_ETYP_SHP(_etyp) == MPCCI_ESHP_POLY )

#define MPCCI_ESHP_IS_PCELL(_eshp)     (       (_eshp)         == MPCCI_ESHP_PCELL )
#define MPCCI_ETYP_IS_PCELL(_etyp)     ( MPCCI_ETYP_SHP(_etyp) == MPCCI_ESHP_PCELL )

#define MPCCI_ETYP_NVX(_etyp)     /* like NVX_BASE, but with poly check */ \
(\
   MPCCI_ETYP_HAS_PSIG(_etyp) /* this is a polygon or pcell */\
      ? MPCCI_ETYP_NNE(_etyp) /* ... and we have NVX := NNE */\
      : MPCCI_ETYP_NVX_BASE(_etyp)\
)

/*
 * The nodeids array of a PCELL is in fact organised into two parts:
 *
 *    a) nodeids[0 ... nodeids[0]-1]   lookup table (header part of nodeids)
 *    b) nodeids[nodeids[0] .... nne]  vertex index list
 *
 *    a)
 *       nodeids[0]              => start index of 1st vertex of face[0] in nodeids
 *       nodeids[1]              => start index of 1st vertex of face[1] in nodeids
 *       nodeids[2]              => start index of 1st vertex of face[2] in nodeids
 *       .........
 *       nodeids[nface-1]        => start index of 1st vertex of last face in nodeids
 *       nodeids[nface]          => index of last vertex+1 (just points behind the last face)
 *
 *    b)
 *       nodeids[nodeids[0]  ]   => vertex id of 1st vertex of face 0
 *       nodeids[nodeids[0]+1]   => vertex id of 2nd vertex of face 0
 *       nodeids[nodeids[0]+2]   => vertex id of 3rd vertex of face 0
 *       .........
 *       nodeids[nodeids[1]  ]   => vertex id of 1st vertex of face 1
 *       nodeids[nodeids[1]+1]   => vertex id of 2nd vertex of face 1
 *       .........
 *
 *       nodeids[0] - 1          => no. of faces in PCELL
 */


#define MPCCI_PCELL_NFACES(_n)      ((_n)[0] - 1 )        /* no. of faces in polyhedron */
#define MPCCI_PCELL_FNVERT(_n,_i)   ((_n)[_i+1]-(_n)[_i]) /* no. of vertices of face i */
#define MPCCI_PCELL_ISTART(_n,_i)   (         (_n)[_i] )  /* index   of 1st vertex of face i */
#define MPCCI_PCELL_VSTART(_n,_i)   ((_n)[(_n)[_i]])      /* value   of 1st vertex of face i */
#define MPCCI_PCELL_PSTART(_n,_i)   ((_n)+(_n)[_i] )      /* pointer to 1st vertex of face i */

#endif
