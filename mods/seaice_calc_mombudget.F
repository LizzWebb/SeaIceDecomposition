#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_CALC_MOMBUDGET
C     !INTERFACE:
      SUBROUTINE SEAICE_CALC_MOMBUDGET(
     I     uIceLoc, vIceLoc,
     &     TILTx, TILTy, TAU_AI_X, TAU_AI_Y,uIce_prev,vIce_prev,
     I     myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_CALC_MOMBUDGET
C     *==========================================================*
C     | written by Elizabeth Webb, March 2025
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"      
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"
#include "EXF_FIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

#if ( defined ALLOW_DIAGNOSTICS && defined SEAICE_CGRID ) 
C     !LOCAL VARIABLES:                                           
C     === Local variables ===                                     
C     i,j    :: Loop counters                                     
C     bi,bj  :: tile counters                                     
      INTEGER i, j, bi, bj
      INTEGER RHSiMax, RHSiMin, RHSjMax, RHSjMin
C     Vertical integers: k from LHS, ksurf from RHS 
      INTEGER kSrf
      
C     fractional area at velocity points
C     Note, we follow the below from LHS      
C      _RL areaW(1:sNx,1:sNy)
C      _RL areaS(1:sNx,1:sNy)
      _RL RHSareaW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RHSareaS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS SINWAT
      _RL COSWAT, recip_deltaT
      _RL RHSareaMinLoc
C     backward difference extrapolation factor (includes 1/deltaT)  
      _RL bdfAlphaOverDt
      
C     symmetric drag coefficient: from LHS 
      _RL RHSdragSym(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      
C     u/vIceLoc :: local copies of the current ice velocity
      _RL uIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uIce_prev(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIce_prev(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     u/vIceLHS :: LHS of momentum equations
      _RL uIceLHSLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLHSLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uIceRHSLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceRHSLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     Residual 
      _RL uIceRESLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceRESLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C      Each term 
      _RL MunewDt (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL MvnewDt (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL MuoldDt (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL MvoldDt (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL coriX (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL coriY (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uDWUi (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vDWUi (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uDWUo (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vDWUo (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL TAU_AI_X   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL TAU_AI_Y   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL TILTx (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL TILTy (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dsigX (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dsigY (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      
C     EAW[DO: Re-calculate RHS and LHS of momentum equation.                           
      CALL SEAICE_CALC_MOMBUDGET_RHS(
     I     uIce_prev, vIce_prev, TILTx, TILTy, TAU_AI_X, TAU_AI_Y,
     O     uIceRHSLoc, vIceRHSLoc, MuoldDt, MvoldDt, uDWUo, vDWUo,
     I     0, myTime, myIter, myThid )

       CALL SEAICE_CALC_MOMBUDGET_LHS(
     I     uIce, vIce,
     O     uIceLHSLoc, vIceLHSLoc, MunewDt, MvnewDt, dsigX, dsigY,                                                                           
     &     coriX, coriY, uDWUi, vDWUi,
     I     myTime, myIter, myThid )

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO J=1,sNy
          DO I=1,sNx
             uIceRESLoc(I,J,bi,bj) =
     &               uIceLHSLoc(I,J,bi,bj) - uIceRHSLoc(I,J,bi,bj)
             vIceRESLoc(I,J,bi,bj) =
     &               vIceLHSLoc(I,J,bi,bj) - vIceRHSLoc(I,J,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C     Writing everything to diagnostics
       IF ( useDiagnostics ) THEN
C         CALL DIAGNOSTICS_FILL(dUdtX,'mdUdtX  ',0,1,0,1,1,myThid)
C         CALL DIAGNOSTICS_FILL(dUdtY,'mdUdtY  ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(MuoldDt,'MuoldDt ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(MunewDt,'MunewDt ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(MvoldDt,'MvoldDt ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(MvnewDt,'MvnewDt ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(coriX,'coriX   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(coriY,'coriY   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(TAU_AI_X,'TAU_AI_X',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(TAU_AI_Y,'TAU_AI_Y',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(uDWUi,'uDWUi   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(vDWUi,'vDWUi   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(uDWUo,'uDWUo   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(vDWUo,'vDWUo   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(TILTx,'TILTx   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(TILTy,'TILTy   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(dsigX,'dsigX   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(dsigY,'dsigY   ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(uIceLHSLoc,'uIceLHS ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(vIceLHSLoc,'vIceLHS ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(uIceRHSLoc,'uIceRHS ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(vIceRHSLoc,'vIceRHS ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(uIceRESLoc,'uIceRES ',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(vIceRESLoc,'vIceRES ',0,1,0,1,1,myThid)
       ENDIF 

#endif /* (  ALLOW_DIAGNOSTICS && defined SEAICE_CGRID )*/
      RETURN
      END



C --------------------------------------------------------------------------------------------
C --------------------------------------------------------------------------------------------                                                 
#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_CALC_MOMBUDGET_RHS
C     !INTERFACE:
      SUBROUTINE SEAICE_CALC_MOMBUDGET_RHS(
     I     uIce_prev, vIce_prev, TILTx, TILTy, TAU_AI_X, TAU_AI_Y,
     O     uIceRHS, vIceRHS, MuoldDt, MvoldDt, uDWUo, vDWUo,
     I     krylovIter, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_CALC_MOMBUDGET_RHS
C     | o Right-hand side of momentum equations after
C     |   solver has updated velocities.
C     *==========================================================*
C     | written by Martin Losch, Oct 2012
C     | modified by Elizabeth Webb, Sept 2025                                                                                                
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C     newtonIter :: current iterate of Newton iteration
C     krylovIter :: current iterate of Krylov iteration
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER newtonIter
      INTEGER krylovIter
C     u/vIceRHS :: RHS of momentum equations
      _RL uIceRHS   (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL vIceRHS   (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL uIce_prev (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL vIce_prev (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL MuoldDt   (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL MvoldDt   (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL uDWUo     (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL vDWUo     (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL TAU_AI_X  (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL TAU_AI_Y  (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL TILTx     (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL TILTy     (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)

      
#if ( defined SEAICE_CGRID \
      && ( defined SEAICE_ALLOW_JFNK || defined SEAICE_ALLOW_KRYLOV ) )
C     i,j,bi,bj,kSrf :: loop indices
      INTEGER i,j,bi,bj
      INTEGER kSrf
      _RS     SINWAT
      _RL     COSWAT
      _RL recip_deltaT
C     fractional area at velocity points
      _RL areaW(1:sNx,1:sNy)
      _RL areaS(1:sNx,1:sNy)
CEOP

      IF ( usingPCoords ) THEN
       kSrf = Nr
      ELSE
       kSrf = 1
      ENDIF
C--   introduce turning angles
      SINWAT=SIN(SEAICE_waterTurnAngle*deg2rad)
      COSWAT=COS(SEAICE_waterTurnAngle*deg2rad)

C     initialise fractional areas at velocity points
      DO J=1,sNy
       DO I=1,sNx
        areaW(I,J) = 1. _d 0
        areaS(I,J) = 1. _d 0
       ENDDO
      ENDDO

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        IF ( SEAICEscaleSurfStress ) THEN
         DO J=1,sNy
          DO I=1,sNx
           areaW(I,J) = 0.5 _d 0*(AREA(I,J,bi,bj)+AREA(I-1,J,bi,bj))
           areaS(I,J) = 0.5 _d 0*(AREA(I,J,bi,bj)+AREA(I,J-1,bi,bj))
          ENDDO
         ENDDO
        ENDIF
        DO J=1,sNy
         DO I=1,sNx
      recip_deltaT = 1. _d 0 / SEAICE_deltaTdyn

C     ice-velocity independent contribution to drag terms
C     - Cdrag*(uVel*cosWat - vVel*sinWat)/(vVel*cosWat + uVel*sinWat)
C     ( remember to average to correct velocity points )
C
      MuoldDt(I,J,bi,bj)   =
     &           seaiceMassU(I,J,bi,bj)*recip_deltaT
     &         *uIce_prev(I,J,bi,bj) *maskinW(I,J,bi,bj)
      uDWUo(I,J,bi,bj)   = ( 
     &         0.5 _d 0 * ( DWATN(I,J,bi,bj)+DWATN(I-1,J,bi,bj) )*
     &         COSWAT * uVel(I,J,kSrf,bi,bj)
     &         - SIGN(SINWAT, _fCori(I,J,bi,bj))* 0.5 _d 0 *
     &         ( DWATN(I  ,J,bi,bj) * 0.5 _d 0 *
     &          ( vVel(I  ,J,kSrf,bi,bj)+vVel(I  ,J+1,kSrf,bi,bj) )
     &         + DWATN(I-1,J,bi,bj) * 0.5 _d 0 *
     &          ( vVel(I-1,J,kSrf,bi,bj)+vVel(I-1,J+1,kSrf,bi,bj) )
     &           ) )*areaW(I,J)*maskinW(I,J,bi,bj)
C     V component of RHS of the model momentum equation 
      MvoldDt(I,J,bi,bj) =
     &         seaiceMassV(I,J,bi,bj)*recip_deltaT
     &         *vIce_prev(I,J,bi,bj) *maskinS(I,J,bi,bj)
      vDWUo(I,J,bi,bj)    = (
     &         0.5 _d 0 * (DWATN(I,J,bi,bj)+DWATN(I,J-1,bi,bj) ) *
     &         COSWAT * vVel(I,J,kSrf,bi,bj)
     &         + SIGN(SINWAT, _fCori(I,J,bi,bj)) * 0.5 _d 0 *
     &         ( DWATN(I,J  ,bi,bj) * 0.5 _d 0 *
     &          ( uVel(I,J  ,kSrf,bi,bj)+uVel(I+1,J  ,kSrf,bi,bj))
     &         + DWATN(I,J-1,bi,bj) * 0.5 _d 0 *
     &          ( uVel(I,J-1,kSrf,bi,bj)+uVel(I+1,J-1,kSrf,bi,bj))
     &     ) )*maskinS(I,J,bi,bj)
C       
         uIceRHS(I,J,bi,bj) = (MuoldDt(I,J,bi,bj) + uDWUo(I,J,bi,bj)
     &     + TILTx(I,J,bi,bj) + TAU_AI_X(I,J,bi,bj) )*maskinW(I,J,bi,bj)
         vIceRHS(I,J,bi,bj) = (MvoldDt(I,J,bi,bj) + vDWUo(I,J,bi,bj)
     &     + TILTy(I,J,bi,bj) + TAU_AI_Y(I,J,bi,bj) )*maskinS(I,J,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* SEAICE_CGRID, SEAICE_ALLOW_JFNK and KRYLOV */
      RETURN
      END


C --------------------------------------------------------------------------------------------                                                
C --------------------------------------------------------------------------------------------

#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_CALC_MOMBUDGET_LHS
C     !INTERFACE:
      SUBROUTINE SEAICE_CALC_MOMBUDGET_LHS(
     I     uIceLoc, vIceLoc,
     O     uIceLHS, vIceLHS,MunewDt,MvnewDt,dsigX,dsigY,                                                                                     
     &     coriX, coriY,uDWUi,vDWUi,
     I     myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_CALC_LHS
C     | o Left-hand side of momentum equations, i.e. all terms
C     |   that depend on the ice velocities of the current
C     |   iterate of the Newton-Krylov iteration
C     |
C     | o The scheme is backward Euler in time, i.e. the
C     |   rhs-vector contains only terms that are independent
C     |   of u/vIce, except for the time derivative part
C     |   mass*(u/vIce-u/vIceNm1)/deltaT
C     | o Left-hand side contributions
C     |   + mass*(u/vIce)/deltaT
C     |   + Cdrag*(uIce*cosWat - vIce*sinWat)
C     |          /(vIce*cosWat + uIce*sinWat)
C     |   - mass*f*vIce/+mass*f*uIce
C     |   - dsigma/dx / -dsigma/dy, eta and zeta are
C     |                   computed only once per Newton iterate
C     *==========================================================*
C     | written by Martin Losch, Oct 2012
C     | modified by Elizabeth Webb, Sept 2025      
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C     newtonIter :: current iterate of Newton iteration
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER newtonIter
C     u/vIceLoc :: local copies of the current ice velocity
      _RL uIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     u/vIceLHS :: LHS of momentum equations
      _RL uIceLHS  (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL vIceLHS  (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL MunewDt  (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL MvnewDt  (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL dsigX    (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL dsigY    (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL coriX    (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL coriY    (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL uDWUi    (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)
      _RL vDWUi    (1-OLx:sNx+OLx, 1-OLy:sNy+OLy, nSx, nSy)

      

#if ( defined SEAICE_CGRID \
      && ( defined SEAICE_ALLOW_JFNK || defined SEAICE_ALLOW_KRYLOV ) )
C     i,j,bi,bj,k :: loop indices
      INTEGER i,j,bi,bj
      INTEGER k
      _RS     SINWAT
      _RL     COSWAT
      _RL recip_deltaT
C     backward difference extrapolation factor
      _RL bdfAlpha
C     symmetric drag coefficient
      _RL dragSym(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C     fractional area at velocity points
      _RL areaW(1:sNx,1:sNy)
      _RL areaS(1:sNx,1:sNy)

#ifdef SEAICE_ALLOW_MOM_ADVECTION
C     tendency due to advection of momentum
      _RL gUmom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gVmom(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif /*  SEAICE_ALLOW_MOM_ADVECTION */
CEOP

      k=1
      recip_deltaT = 1. _d 0 / SEAICE_deltaTdyn
C--   introduce turning angles
      SINWAT=SIN(SEAICE_waterTurnAngle*deg2rad)
      COSWAT=COS(SEAICE_waterTurnAngle*deg2rad)
C     backward difference extrapolation factor
      bdfAlpha = 1. _d 0
      IF ( SEAICEuseBDF2 ) THEN
       IF ( myIter.EQ.nIter0 .AND. SEAICEmomStartBDF.EQ.0 ) THEN
        bdfAlpha = 1. _d 0
       ELSE
        bdfAlpha = 1.5 _d 0
       ENDIF
      ENDIF

C     initialise fractional areas at velocity points
      DO J=1,sNy
       DO I=1,sNx
        areaW(I,J) = 1. _d 0
        areaS(I,J) = 1. _d 0
       ENDDO
      ENDDO

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C     symmetric drag coefficient may include bottomdrag for grounded ice
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          dragSym(I,J) = DWATN(I,J,bi,bj)*COSWAT
#ifdef SEAICE_ALLOW_BOTTOMDRAG
     &         +CbotC(I,J,bi,bj)
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
         ENDDO
        ENDDO
C
C     compute components of stress tensor from current velocity field
C     and compute divergence of stress tensor
C
        CALL SEAICE_CALC_STRESSDIV(
     I       e11, e22, e12, press, zeta, eta, etaZ,
     O       dsigX, dsigY,
     I       bi, bj, myTime, myIter, myThid )
C     compute lhs side of momentum equations
        IF ( SEAICEscaleSurfStress ) THEN
         DO J=1,sNy
          DO I=1,sNx
           areaW(I,J) = 0.5 _d 0*(AREA(I,J,bi,bj)+AREA(I-1,J,bi,bj))
           areaS(I,J) = 0.5 _d 0*(AREA(I,J,bi,bj)+AREA(I,J-1,bi,bj))
          ENDDO
         ENDDO
        ENDIF
        DO J=1,sNy
         DO I=1,sNx
C     mass*(uIce)/deltaT - dsigma/dx
          dsigX(I,J,bi,bj) = -dsigX(I,J,bi,bj)
     &                 * maskinW(I,J,bi,bj)
          MunewDt(I,J,bi,bj) =
     &         seaiceMassU(I,J,bi,bj)*recip_deltaT
     &         *uIceLoc(I,J,bi,bj) *maskinW(I,J,bi,bj)
C     coriols terms: - mass*f*vIce
          coriX(I,J,bi,bj) = (- 0.5 _d 0*(
     &         seaiceMassC(I  ,J,bi,bj) * _fCori(I  ,J,bi,bj)
     &       * 0.5 _d 0*( vIceLoc(I  ,J,bi,bj)+vIceLoc(I  ,J+1,bi,bj) )
     &       + seaiceMassC(I-1,J,bi,bj) * _fCori(I-1,J,bi,bj)
     &       * 0.5 _d 0*( vIceLoc(I-1,J,bi,bj)+vIceLoc(I-1,J+1,bi,bj) )
     &           ))*maskinW(I,J,bi,bj)
C     ocean-ice and bottom drag terms: + (Cdrag*cosWat+Cb)*uIce - vIce*sinWat)
          uDWUi(I,J,bi,bj) =  (
     &         0.5 _d 0 * ( dragSym(I,J)+dragSym(I-1,J) )
     &         * uIceLoc(I,J,bi,bj)
     &         - SIGN(SINWAT, _fCori(I,J,bi,bj))* 0.5 _d 0 *
     &         (  DWATN(I  ,J,bi,bj) * 0.5 _d 0 *
     &         (vIceLoc(I  ,J,bi,bj)+vIceLoc(I  ,J+1,bi,bj))
     &         +  DWATN(I-1,J,bi,bj) * 0.5 _d 0 *
     &         (vIceLoc(I-1,J,bi,bj)+vIceLoc(I-1,J+1,bi,bj))
     &         ) ) *maskinW(I,J,bi,bj)
          
C     mass*(vIce)/deltaT - dsigma/dy
          dsigY(I,J,bi,bj) = - dsigY(I,J,bi,bj)
     &                 *maskinS(I,J,bi,bj)
          MvnewDt(I,J,bi,bj) =
     &         seaiceMassV(I,J,bi,bj)*recip_deltaT
     &         *vIceLoc(I,J,bi,bj)*maskinS(I,J,bi,bj)
C     + mass*f*uIce
          coriY(I,J,bi,bj) = (0.5 _d 0*(
     &         seaiceMassC(I,J  ,bi,bj) * _fCori(I,J  ,bi,bj)
     &       * 0.5 _d 0*( uIceLoc(I,J  ,bi,bj)+uIceLoc(I+1,  J,bi,bj) )
     &       + seaiceMassC(I,J-1,bi,bj) * _fCori(I,J-1,bi,bj)
     &       * 0.5 _d 0*( uIceLoc(I,J-1,bi,bj)+uIceLoc(I+1,J-1,bi,bj) )
     &           ))*maskinS(I,J,bi,bj)
C    + (Cdrag*cosWat+Cb)*uIce + uIce*sinWat)
          vDWUi(I,J,bi,bj) =  (
     &         0.5 _d 0 * ( dragSym(I,J)+dragSym(I,J-1) )
     &         * vIceLoc(I,J,bi,bj)
     &         + SIGN(SINWAT, _fCori(I,J,bi,bj)) * 0.5 _d 0 *
     &         (  DWATN(I,J  ,bi,bj) * 0.5 _d 0 *
     &         (uIceLoc(I,J  ,bi,bj)+uIceLoc(I+1,J  ,bi,bj))
     &         +  DWATN(I,J-1,bi,bj) * 0.5 _d 0 *
     &         (uIceLoc(I,J-1,bi,bj)+uIceLoc(I+1,J-1,bi,bj))
     &         ) ) * areaS(I,J)*maskinS(I,J,bi,bj)
C     apply masks for interior (important when we have open boundaries)
          uIceLHS(I,J,bi,bj) = dsigX(I,J,bi,bj) + MunewDt(I,J,bi,bj) +
     &         coriX(I,J,bi,bj) + uDWUi(I,J,bi,bj)
          vIceLHS(I,J,bi,bj) = dsigY(I,J,bi,bj) + MvnewDt(I,J,bi,bj) +
     &         coriY(I,J,bi,bj) + vDWUi(I,J,bi,bj)
         ENDDO
        ENDDO
#ifdef SEAICE_ALLOW_MOM_ADVECTION
        IF ( SEAICEmomAdvection ) THEN
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
           gUmom(I,J) = 0. _d 0
           gVmom(I,J) = 0. _d 0
          ENDDO
         ENDDO
         CALL SEAICE_MOM_ADVECTION(
     I        bi,bj,1,sNx,1,sNy,
     I        uIceLoc, vIceLoc,
     O        gUmom, gVmom,
     I        myTime, myIter, myThid )
C     Beware of sign! gU/Vmom is computed for the rhs of the equation;
C     therefore, we need to substract gU/Vmom from the left hand side
         DO J=1,sNy
          DO I=1,sNx
           uIceLHS(I,J,bi,bj) = uIceLHS(I,J,bi,bj) - gUmom(I,J)
           vIceLHS(I,J,bi,bj) = vIceLHS(I,J,bi,bj) - gVmom(I,J)
          ENDDO
         ENDDO
        ENDIF
#endif /* SEAICE_ALLOW_MOM_ADVECTION */
       ENDDO
      ENDDO

#endif /* SEAICE_CGRID, SEAICE_ALLOW_JFNK and KRYLOV */

      RETURN
      END


