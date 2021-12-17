! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_AOT
!
!> \brief MPAS ocean explicit bottom drag
!> \author Mark Petersen
!> \date   August 2017
!> \details
!>  This module contains the routine for computing
!>  tendencies from explicit bottom drag.
!
!-----------------------------------------------------------------------

module ocn_vel_aot

   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_forcing
!   use ocn_time_varying_forcing
   use ocn_diagnostics_variables
   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vel_aot_tend, &
             ocn_vel_aot_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      aotOff  ! on/off switch for explicit bottom drag

   real (kind=RKIND) :: &
      aotCoeff          ! drag coefficient

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_AOT_tend
!
!> \brief   Computes tendency term from explicit bottom drag
!> \author  Mark Petersen
!> \date    15 August 2017
!> \details
!>  This routine computes the explicit bottom drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_aot_tend(timeVaryingForcingPool,forcingPool, tend, err)
!normVelocity, KECell, layerThickEdge, tend, err) !{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
!      type (domain_type) intent(inout) :: domain
      type (mpas_pool_type), intent(in) :: timeVaryingForcingPool
      type (mpas_pool_type), intent(in) :: forcingPool
!!      real (kind=RKIND), dimension(:,:), intent(in) :: &
!!         normVelocity,  &!< [in] normal velocity
!!         KECell,        &!< [in] kinetic energy at cell
!!         layerThickEdge  !< [in] layer thickness at edge
!      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: &
!          velocityMeridional, &!< [in] meridional velocity at cell centers
!          velocityZonal,       &!< [in] zonal velocity at cell centers 
!          velocityMeridionalObs, &!< [in] meridional velocity at cell centers
!          velocityZonalObs        !< [in] zonal velocity at cell centers 
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iEdge,      &! loop index for edge loop
         k,          &! vertical index of lowest active layer at edge
         j,          &! 
         iSmooth,    &
         iCell,      &
         coc,        &
         eoc,        &
         nCells,     &
         cell1, cell2 ! neighbor cell addresses across edge
      real (kind=RKIND) :: aotCounter, aotCounterOffset
      integer, dimension(:), allocatable :: aotLayer ! tracks cell used for aot smoothing term
!      real (kind=RKIND), dimension(:,:), allocatable :: velocityZonalAOTNew, velocityMeridionalAOTNew
      real (kind=RKIND), dimension(:,:), pointer :: &
         velocityMeridionalObs, &!< [in] meridional velocity at cell centers
         velocityZonalObs        !< [in] zonal velocity at cell centers 
!         velocityMeridionalAOT, &!< [in] interpolated meridional velocity for AOT
!         velocityZonalAOT        !< [in] interpolated zonal velocity for AOT 
      integer, dimension(:), pointer :: aotForcingMask
      real (kind=RKIND), dimension(:,:), allocatable :: &
         velocityMeridionalAOT, &!< [in] interpolated meridional velocity for AOT
         velocityZonalAOT,      &!< [in] interpolated zonal velocity for AOT 
         velocityMeridionalAOTNew, &!< [in] interpolated meridional velocity for AOT
         velocityZonalAOTNew        !< [in] interpolated zonal velocity for AOT 
      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error code and return if not turned on
      !*** Otherwise start timer

      err = 0
      if (aotOff) return
      call mpas_timer_start('AOT u-v term')
      
      ! Explicit bottom drag term:
      ! du/dt = ... - c |u| u / h
      ! appied to bottom layer only.
      ! This term comes from the bottom boundary condition in the vertical
      ! momentum mixing, and is explicit if both |u| and u are chosen to be at
      ! time level n.

      call mpas_pool_get_array(timeVaryingForcingPool, 'velocityZonalObs', velocityZonalObs)
      call mpas_pool_get_array(timeVaryingForcingPool, 'velocityMeridionalObs', velocityMeridionalObs)
      
      call mpas_pool_get_array(forcingPool, 'aotForcingMask', aotForcingMask)


!      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      nCells = nCellsAll
!      allocate(velocityZonalAOTNew(nVertLevels,nCells))
!      allocate(velocityMeridionalAOTNew(nVertLevels,nCells))
      allocate(aotLayer(nCells))
!      allocate(aotCounter(nCells))
      allocate(velocityMeridionalAOT(nVertLevels,nCells))
      allocate(velocityZonalAOT(nVertLevels,nCells))
      allocate(velocityMeridionalAOTNew(nVertLevels,nCells))
      allocate(velocityZonalAOTNew(nVertLevels,nCells))


!      write(*,*) "velocicyZonalObs"
!      write(*,*) velocityZonalObs(1,1:10)
!      write(*,*) "velocicyMeridionalObs"
!      write(*,*) velocityMeridionalObs(1,1:10)
!      write(*,*) "velocicyZonal"
!      write(*,*) velocityZonal(1,1:10)
!      write(*,*) "velocicyMeridional"
!      write(*,*) velocityMeridional(1,1:10)
!      write(*,*) "mu"
!      write(*,*) aotCoeff
!      write(*,*) aotForcingMask(1:10)

!initialize aotLayer for non-masked cells
      do iCell = 1, nCells
        if(aotForcingMask(iCell) == 1) then
           aotLayer(iCell) = 1 !non-masked cells are in the first layer, which increments by 1 going outward to neighbors
           do k = 1, maxLevelCell(iCell)
             velocityMeridionalAOT(k,iCell) = velocityMeridionalObs(k,iCell) - velocityMeridional(k,iCell)
             velocityZonalAOT(k,iCell) = velocityZonalObs(k,iCell) - velocityZonal(k,iCell)
           enddo !k
        else
          aotLayer(iCell) = 0
            do k = 1, maxLevelCell(iCell)
             velocityMeridionalAOT(k,iCell) = 0
             velocityZonalAOT(k,iCell) = 0
           enddo
        endif !aotForcingMask(iCell)
      enddo !iCell

      ! Nearest Neighbor initializeation
      if(config_aot_interp_weight == "constant") then
        aotCounterOffset = 1
        do iSmooth = 1, config_aot_neighbor_iterations
          do iCell = 1, nCells
            if (aotLayer(iCell) == iSmooth) then
              do j = 1, nEdgesOnCell(iCell)
                coc = cellsOnCell(j,iCell)
                if (coc < nCells+1) then
                  if (aotLayer(coc) == 0 .OR. aotLayer(coc) == iSmooth + 1) then
                    aotLayer(coc) = iSmooth + 1
                    !aotCounter = aotCounter + 1
                    do k = 1, maxLevelCell(iCell)
                      velocityMeridionalAOT(k,coc) = (velocityMeridionalAOT(k,coc)*aotCounter &
                                                   + velocityMeridionalAOT(k,iCell))/(aotCounter+aotCounterOffset)
                      velocityZonalAOT(k,coc) = (velocityZonalAOT(k,coc)*aotCounter &
                                              + velocityZonalAOT(k,iCell))/(aotCounter+aotCounterOffset)
                    enddo !k
                    aotCounter = aotCounter + aotCounterOffset
                  endif !aotLayer(coc)
                endif !coc
              enddo !j
            endif !aotLayer(iCell)
          enddo !iCell
        enddo !iSmooth
      else if(config_aot_interp_weight == "geometric") then
        do iSmooth = 1, config_aot_neighbor_iterations
          do iCell = 1, nCells
            if (aotLayer(iCell) == iSmooth) then
              do j = 1, nEdgesOnCell(iCell)
                coc = cellsOnCell(j,iCell)
                eoc = edgesOnCell(j,iCell)
                if (coc < nCells+1) then
                  if (aotLayer(coc) == 0 .OR. aotLayer(coc) == iSmooth + 1) then
                    aotLayer(coc) = iSmooth + 1
                    aotCounterOffset = dcEdge(eoc)
                    !aotCounter = aotCounter + 1
                    do k = 1, maxLevelCell(iCell)
                      velocityMeridionalAOT(k,coc) = (velocityMeridionalAOT(k,coc)*aotCounter &
                                                   + velocityMeridionalAOT(k,iCell))/(aotCounter+aotCounterOffset)
                      velocityZonalAOT(k,coc) = (velocityZonalAOT(k,coc)*aotCounter &
                                              + velocityZonalAOT(k,iCell))/(aotCounter+aotCounterOffset)
                    enddo !k
                    aotCounter = aotCounter + aotCounterOffset
                  endif !aotLayer(coc)
                endif !coc
              enddo !j
            endif !aotLayer(iCell)
          enddo !iCell
        enddo !iSmooth
      endif

      ! Once data has been initialized apply iterated averaging to smooth it
      if(config_aot_interp_weight == "constant") then
        aotCounterOffset = 1
        do iSmooth = 1, config_aot_smoothing_iterations
          do iCell = 1, nCells
            if (aotForcingMask(iCell) == 0) then
              velocityZonalAOTNew(:,iCell) = velocityZonalAOT(:,iCell)
              velocityMeridionalAOTNew(:,iCell) = velocityMeridionalAOT(:,iCell)
              aotCounter = 1
              do j = 1, nEdgesOnCell(iCell)
                coc = cellsOnCell(j,iCell)
                if (coc < nCells+1) then
                  velocityMeridionalAOTNew(:,iCell) = velocityMeridionalAOTNew(:,iCell) + velocityMeridionalAOT(:,coc)
                  velocityZonalAOTNew(:,iCell) = velocityZonalAOTNew(:,iCell) + velocityZonalAOT(:,coc)
                  aotCounter = aotCounter + aotCounterOffset
                endif !coc
              enddo !j
            endif !aotLayer(iCell)
            velocityMeridionalAOTNew(:,iCell) = velocityMeridionalAOTNew(:,iCell)/aotCounter
            velocityZonalAOTNew(:,iCell) = velocityZonalAOTNew(:,iCell)/aotCounter
          enddo !iCell
          velocityMeridionalAOT = velocityMeridionalAOTNew
          velocityZonalAOT = velocityZonalAOTNew
        enddo !iSmooth
      else if(config_aot_interp_weight == "geometric") then
        do iSmooth = 1, config_aot_smoothing_iterations
          do iCell = 1, nCells
            if (aotForcingMask(iCell) == 0) then
              velocityZonalAOTNew(:,iCell) = velocityZonalAOT(:,iCell)
              velocityMeridionalAOTNew(:,iCell) = velocityMeridionalAOT(:,iCell)
              aotCounter = 1
              do j = 1, nEdgesOnCell(iCell)
                coc = cellsOnCell(j,iCell)
                if (coc < nCells+1) then
                  velocityMeridionalAOTNew(:,iCell) = velocityMeridionalAOTNew(:,iCell) + velocityMeridionalAOT(:,coc)
                  velocityZonalAOTNew(:,iCell) = velocityZonalAOTNew(:,iCell) + velocityZonalAOT(:,coc)
                  aotCounterOffset = dcEdge(eoc)
                  aotCounter = aotCounter + aotCounterOffset
                endif !coc
              enddo !j
            endif !aotLayer(iCell)
            velocityMeridionalAOTNew(:,iCell) = velocityMeridionalAOTNew(:,iCell)/aotCounter
            velocityZonalAOTNew(:,iCell) = velocityZonalAOTNew(:,iCell)/aotCounter
          enddo !iCell
          velocityMeridionalAOT = velocityMeridionalAOTNew
          velocityZonalAOT = velocityZonalAOTNew
        enddo !iSmooth
      endif









!      call mpas_log_write("Testing observational data loading: $r", MPAS_LOG_WARN, realArgs=(/ velocityZonalObs(1,1,1)))
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, tend, latCell,&
      !$acc            velocityMeridionalAOT, velocityZonalAOT) &
      !$acc    private(k, cell1, cell2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(, cell1, cell2)
#endif
      do iEdge = 1, nEdgesOwned
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
!        k =  maxLevelEdgeTop(iEdge)
        if (latCell(cell1)<=1.48353 .AND. latCell(cell2) <=1.48353) then
          do k = 1, maxLevelEdgeTop(iEdge)
!        if (k > 0) then
!           tend(k,iEdge) = tend(k,iEdge) - aotCoeff* &
!                           sqrt(KECell(k,cell1) + KECell(k,cell2))* &
!                           normVelocity(k,iEdge)/layerThickEdge(k,iEdge)
!        end if
            tend(k,iEdge) = tend(k,iEdge) + aotCoeff*&
              (  sin(angleEdge(iEdge)) * (0.5*(velocityMeridionalAOT(k,cell1) + velocityMeridionalAOT(k,cell2))) + &
                 cos(angleEdge(iEdge)) * (0.5*(velocityZonalAOT(k,cell1) + velocityZonalAOT(k,cell2)))     )
          enddo
        end if
      enddo
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      deallocate( aotLayer, velocityZonalAOT, velocityMeridionalAOT, velocityZonalAOTNew, velocityMeridionalAOTNew)
      call mpas_timer_stop('AOT u-v term')

   !--------------------------------------------------------------------

   end subroutine ocn_vel_aot_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_AOT_init
!
!> \brief   Initializes ocean explicit bottom drag forcing
!> \author  Mark Petersen
!> \date    August 2017
!> \details
!>  This routine initializes quantities related to explicit bottom drag
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_aot_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize return error code and set module defaults
      err = 0

      aotOff  = .true.
      aotCoeff = 0.0_RKIND


      !*** Reset values based on input model configuration

      if (config_use_aot) then
         aotOff = .false.
         aotCoeff = config_aot_coeff
!         write(*,*) "here"
      endif

      if (config_disable_vel_aot) then
         aotOff = .true.
         aotCoeff = 0.0_RKIND
!         write(*,*) "no, here"
      endif
      
!      write(*,*) "Mu:"
!      write(*,*) aotCoeff
   !--------------------------------------------------------------------

   end subroutine ocn_vel_aot_init!}}}

!***********************************************************************

end module ocn_vel_aot

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
