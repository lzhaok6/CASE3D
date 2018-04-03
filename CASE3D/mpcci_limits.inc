C****************************************************************************************
C mpcci_limits.inc
C
C Purpose:
C    Defines some basic limits used to dimension some arrays.
C
C    In principle (by design) there are no limits (except memory), neither
C    for the MpCCI server nor the client library since all objects are managed as chains.
C    However for performance reasons sometimes the chained objects (or pointers)
C    temporarily need to be copied/compressd/collected into a simple vector.
C    In this case the arrays needs to have a dimension which is defined here.
C
C    In case of trouble or error messages simply increase the limits and rebuild
C    the server and client and all adapters.
C
C Reviews/changes:
C    2008/Mar: Carsten Dehning, Initial release
C    $Id: mpcci_limits.inc 104 2013-08-16 11:56:59Z pbayrasy $
C
C****************************************************************************************
C

C
C Define the max. no. of coupled codes handles be the server.
C This parameter is only used for dimensioning some temporary arrays used within
C collect operations.
C
      INTEGER   MPCCI_NCODES_MAX
      PARAMETER(MPCCI_NCODES_MAX=64)


C
C Define the max. no. of quantities send or received on a single mesh.
C This parameter is only used for dimensioning arrays within some structures
C with fixed length to avoid malloc() and free() and temporary arrays used within
C the "collect special quantities" operations
C
      INTEGER   MPCCI_MESHQ_MAX
      PARAMETER(MPCCI_MESHQ_MAX=32)


C
C Define the max. no. of chars used for a mesh/part name.
C This parameter is only used for dimensioning some char arrays within some structures
C with fixed length to avoid strdup() and free()
C
      INTEGER   MPCCI_MESHNAME_MAX
      PARAMETER(MPCCI_MESHNAME_MAX=128)


C
C Define the max. no. of chars used for quantity names.
C This parameter is only used for dimensioning some char arrays within some structures
C with fixed length to avoid strdup() and free()
C
      INTEGER   MPCCI_QUANTNAME_MAX
      PARAMETER(MPCCI_QUANTNAME_MAX=64)


C
C Define the max. no. of chars used for the auxialiary string in the MPCCI_PINFO.
C This parameter is only used for dimensioning some char arrays within some structures
C with fixed length to avoid strdup() and free()
C
      INTEGER   MPCCI_AUXSTRING_MAX
      PARAMETER(MPCCI_AUXSTRING_MAX=256)


C
C Define the max. no. of faces on a polyhedron cell (PCELL)
C
      INTEGER   MPCCI_PCELLFACES_MAX
      PARAMETER(MPCCI_PCELLFACES_MAX=256)


C
C Define the max. no. of vertices of a polygonal face (POLY)
C
      INTEGER   MPCCI_POLYVERTICES_MAX
      PARAMETER(MPCCI_POLYVERTICES_MAX=256)


C
C Define the max. no. of user pointers in the MPCCI_QUANT & MPCCI_PART
C
      INTEGER   MPCCI_USERPOINTERS_MAX
      PARAMETER(MPCCI_USERPOINTERS_MAX=8)
