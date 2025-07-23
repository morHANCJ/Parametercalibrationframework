!     #######################################################
!     ##                                                                                                      ##
!     ##    Urban Flood Model                                                                  ##
!     ##                                                                                                      ##
!     #######################################################
!     ======= HAN Congji                     ======================= 
!---- HLLC Solver & Hancock: Toro E F. Shock-capturing methods for ---
!----     free-surface shallow flows[M]. Wiley and Sons Ltd., 2001. ---
!---- Variable Reconstruction & Slope Limiter: Jawahar P, Kamath H. --- 
!----     A high-resolution procedure for Euler and Navier–Stokes ---
!----     computations on unstructured grids[J]. Journal of Computational --- 
!----     Physics, 2000, 164(1): 165-203. --- 
!---- Xia, X., Liang, Q., Ming, X., & Hou, J. (2017). An efficient and ---
!----     stable hydrodynamic model with novel source term discretization ---
!----     schemes for overland flow and flood simulations. ---
!----     Water resources research, 53(5), 3730-3759. 

module globals_condition
  integer type_r, type_fl
  integer type_tech, type_techm, type_tecbs, type_tecinf 
  real(8) timmax 
  real(8) dt, dkout, dpout 
  integer num_inf
  integer , allocatable :: n_inf(:) 
  real(8) , allocatable :: inf_mn(:), inf_cd(:) 
  real(8) gg, dt2, fita, pi
  real(8) time 
  integer mstep 
end module globals_condition

module globals_meshdata
  integer node, mesh, link 
  integer , allocatable :: ko(:) 
  integer , allocatable :: inf(:) 
  integer , allocatable :: buildmark(:)
  integer , allocatable :: menode(:,:), melink(:,:), limesh(:,:), linode(:,:)
  real(8) , allocatable :: baseo(:)
  real(8) , allocatable :: dnox(:), dnoy(:) 
  real(8) , allocatable :: smesh(:), scv(:), rthl(:, :), ux(:), uy(:) 
  real(8) , allocatable :: mn(:), cd(:), lambda(:), rbeta(:), aj(:)
  real(8) , allocatable :: xmesh(:), ymesh(:), rtuv(:, :) 
end module globals_meshdata

module globals_ground
  real(8) :: th = 1.0d-4
  real(8) , allocatable :: h(:), ho(:), hmax(:)
  real(8) , allocatable :: hb(:), hbo(:)
  real(8) , allocatable :: um(:), umo(:)
  real(8) , allocatable :: vn(:), vno(:)
  real(8) , allocatable :: uum(:), vvm(:)
  real(8) , allocatable :: hpre(:), upre(:), vpre(:), ufpre(:), vfpre(:)
  real(8) , allocatable :: deltabx(:), deltaby(:), deltaux(:), deltauy(:), deltavx(:), deltavy(:)
  real(8) , allocatable :: wedeltabx(:), wedeltaby(:), wedeltaux(:), wedeltauy(:), wedeltavx(:), wedeltavy(:)
  real(8) , allocatable :: reh(:,:), reu(:,:), rev(:,:), reb(:,:), rehh(:,:), bf(:), bfk(:,:), db(:)
  real(8) , allocatable :: h11(:), u11(:), v11(:), u11g(:), v11g(:), u13h(:), v13h(:)
  real(8) , allocatable :: sumh(:), sumu(:), sumv(:)
  real(8) , allocatable :: corrh(:), corru(:), corrv(:)
  real(8) , allocatable :: opx(:), opy(:), opz(:), fm(:), recz(:)
  integer :: op
  integer , allocatable :: opmesh(:)
end module globals_ground
    
module globals_rain
  real(8) dtrain, vrain
  real(8) , allocatable :: rr(:), rain(:)
end module globals_rain

module globals_kyokai
  integer qnb, inl1, inl2
  real(8) dtq, vkyokai
  real(8) , allocatable :: qin1(:), qin2(:) 
end module globals_kyokai
    
module globals_okyokai
    integer oid
    integer, allocatable :: outl(:)
end module globals_okyokai
 !===========================   
  
module subprogs

  implicit none
  contains
!============================================================
!                        READ DATA
!============================================================
!====================================
!         Ground Mesh Data
!====================================
  subroutine rdat 
  use globals_condition
  use globals_meshdata
  implicit none
  integer no, me, li, i, tmp, k, k2
  real(8) dx1, dy1, llink, stmp
!=== NODE ===========================
  read(50, *) node 
  allocate(dnox(node), dnoy(node)) 
  read(50, *)
  do no = 1, node
    read(50, *) tmp, dnox(no), dnoy(no) 
  enddo
!=== MESH ===========================
  read(53, *) mesh 
  allocate(ko(mesh), menode(mesh, 3)) 
  allocate(melink(mesh, 3), smesh(mesh), xmesh(mesh), ymesh(mesh))
  allocate(rtuv(mesh, 3))
    do me = 1, mesh
      read(53, *) tmp, ko(me), (menode(me, k), k = 1, ko(me)) 
      read(53, *) (melink(me, k), k = 1, ko(me))
      read(53, *) smesh(me), xmesh(me), ymesh(me) 
      read(53, *) (rtuv(me, k), k = 1, ko(me)) 
    enddo
!=== LINK ===========================
  read(54, *) link 
  allocate(limesh(link, 2), linode(link, 2))
  allocate(scv(link), rthl(link, 2))
  allocate(ux(link), uy(link))
  do li = 1, link
    read(54, *) tmp, limesh(li, 1), limesh(li, 2), linode(li, 1), linode(li, 2) 
    read(54, *) scv(li), rthl(li, 1), rthl(li, 2) 
    read(54, *) ux(li), uy(li) 
  enddo
!=== INF ============================
  allocate(inf(mesh))
  read(51, *)
  read(51, *)
  do me = 1, mesh
    read(51, *) tmp, inf(me) 
  enddo
  allocate(buildmark(link))
  buildmark = 0
  ! If the BH method is used 
  !do li = 1, link
  !    if (limesh(li, 2) /= 0 .and. inf(limesh(li, 1)) == 2 .and. inf(limesh(li, 1)) /= inf(limesh(li, 2))) buildmark(li) = 1
  !    if (limesh(li, 2) /= 0 .and. inf(limesh(li, 2)) == 2 .and. inf(limesh(li, 1)) /= inf(limesh(li, 2))) buildmark(li) = 1
  !enddo
!=== BASEO ==========================
  allocate(baseo(mesh))
  do me = 1, mesh
    read(52, *) tmp, baseo(me) 
  enddo
!=== INPUT DATA ======================
  read(11, *) inf_mn(1), inf_mn(3), inf_mn(4), inf_cd(4)
  write(*, *) 'mn',inf_mn 
  write(*, *) 'cd', inf_cd
  allocate(cd(mesh), mn(mesh))
  do me = 1, mesh
      do i = 1, num_inf
          if (inf(me) == n_inf(i)) then
              cd(me) = inf_cd(i)
              mn(me) = inf_mn(i)
          endif
      enddo
  enddo
!=== LAMBDA =========================
  allocate(lambda(mesh))
  do me = 1, mesh
      read(201, *) lambda(me) ! Here lambda represents the building coverage ratio, and (1 - lambda) represents the porosity.
  enddo
!=== RBETA ==========================
  allocate(rbeta(link), aj(mesh))
  do li = 1, link
    rbeta(li) = 1.0d0
    if(limesh(li, 2) == 0) then 
      rbeta(li) = 1.0d0 
    elseif(limesh(li, 1) /= 0 .and. limesh(li, 2) /= 0) then 
        rbeta(li) = 0.5d0 * (1-sqrt(lambda(limesh(li, 1))) + 1-sqrt(lambda(limesh(li, 2))))
    endif
  enddo
  do me = 1, mesh
      lambda(me) = 1 - lambda(me)
      stmp = 0.0d0
      do k =1, ko(me)
          k2 = mod(k, ko(me)) + 1
          dx1 = dnox(menode(me, k2)) - dnox(menode(me, k))
          dy1 = dnoy(menode(me, k2)) - dnoy(menode(me, k))
          llink = sqrt(dx1**2 + dy1**2)
          stmp =  stmp + (1-rbeta(melink(me, k))) * llink
      enddo 
      aj(me) = stmp/(3*smesh(me))
  enddo
  end subroutine rdat
!====================================
!       External Volume Data
!====================================
!=== RAINDATA =======================
  subroutine rraindat 
  use globals_rain
  implicit none
  integer irain, i
  read(56, *) irain, dtrain 
  write(*, *) 'rain', irain, dtrain
  allocate(rain(irain))
  do i = 1, irain
    read(56, *) rain(i) 
  enddo
  endsubroutine rraindat
!=== KYOKAI VOLUME ==================
  subroutine rkyokaidat 
  use globals_kyokai
  use globals_meshdata
  implicit none
  integer id, i
  read(57, *) qnb, dtq
  read(57, *) inl1, inl2
  allocate(qin1(qnb), qin2(qnb))
  do i = 1, qnb
    read(57, *) qin1(i), qin2(i)  
  enddo
  rbeta(inl1) = 1.0d0
  rbeta(inl2) = 1.0d0
  endsubroutine rkyokaidat
!!=== OUTFLOW =======================
!  subroutine rokyokaidat
!  use globals_okyokai
!  use globals_meshdata
!  implicit none
!  integer i
!  read(58, *) oid
!  allocate(outl(oid))
!  do i = 1, oid
!      read(58, *) outl(i)
!  enddo
!  rbeta(outl) = 1.0d0
!  endsubroutine rokyokaidat
!====================================
!          Observation Points
!====================================
  subroutine observationpoints
  use globals_ground
  use globals_meshdata
  implicit none
  integer i, j
  real(8) alpha, beta, gamma, denom
  read(202, *) op
  allocate(opx(op), opy(op), opz(op), fm(op), opmesh(op), recz(op))
  opmesh = 1
  do i = 1, op
      read(202, *) opx(i), opy(i), opz(i), fm(i)
      do j = 1, mesh ! Find the meshid of op
          denom = (dnoy(menode(j, 2)) - dnoy(menode(j, 3)))*(dnox(menode(j, 1)) - dnox(menode(j, 3))) + (dnox(menode(j, 3)) - dnox(menode(j, 2)))*(dnoy(menode(j, 1)) - dnoy(menode(j, 3)))
          alpha = ((dnoy(menode(j, 2)) - dnoy(menode(j, 3)))*(opx(i) - dnox(menode(j, 3))) + (dnox(menode(j, 3)) - dnox(menode(j, 2)))*(opy(i) - dnoy(menode(j, 3)))) / denom
          beta  = ((dnoy(menode(j, 3)) - dnoy(menode(j, 1)))*(opx(i) - dnox(menode(j, 3))) + (dnox(menode(j, 1)) - dnox(menode(j, 3)))*(opy(i) - dnoy(menode(j, 3)))) / denom
          gamma = 1.0 - alpha - beta
          if (alpha >= 0 .and. beta >= 0 .and. gamma >= 0) then
              opmesh(i) = j
              if (opz(i) > 0.1d0) baseo(j) = opz(i)
              exit
          endif
      enddo
  enddo
  endsubroutine observationpoints 
!============================================================
!                   Allocate array size
!============================================================
!=== GROUND =========================
  subroutine allocate_ground
  use globals_meshdata
  use globals_ground
  !use globals_concentration
  implicit none
  allocate(h(mesh),  ho(mesh), hmax(mesh)) 
  allocate(hb(mesh), hbo(mesh)) 
  allocate(um(mesh), umo(mesh), uum(mesh))
  allocate(vn(mesh), vno(mesh), vvm(mesh))
  allocate(hpre(mesh), upre(mesh), vpre(mesh), ufpre(mesh), vfpre(mesh))
  allocate(deltabx(mesh), deltaby(mesh), deltaux(mesh), deltauy(mesh), deltavx(mesh), deltavy(mesh))
  allocate(wedeltabx(mesh), wedeltaby(mesh), wedeltaux(mesh), wedeltauy(mesh), wedeltavx(mesh), wedeltavy(mesh))
  allocate(reh(mesh,3), reu(mesh,3), rev(mesh,3), reb(mesh,3), rehh(mesh,3))
  allocate(bf(link), db(link), bfk(mesh,3))
  allocate(h11(link), u11(link), v11(link), u11g(link), v11g(link), u13h(mesh), v13h(mesh))
  allocate(sumh(mesh), sumu(mesh), sumv(mesh))
  allocate(corrh(link), corru(link), corrv(link))
  end subroutine allocate_ground
!=== RAINFALL =======================
  subroutine allocate_rain
  use globals_meshdata
  use globals_rain
  implicit none
  allocate(rr(mesh))
  endsubroutine allocate_rain
  
!============================================================
!                      Initial Value
!============================================================
!=== COMMON INITIAL CONDITIONS ======
  subroutine initial
  use globals_condition
  use globals_meshdata
  use globals_ground
  use globals_rain
  implicit none
  integer i, me, me1, me2, me3, tme, k
  real(8) x1, x2, x3, y1, y2, y3, denominator, rx, ry
  
  do me = 1, mesh   
    h(me) = 0.0d0
    ho(me) = 0.0d0
    hb(me) = baseo(me) + h(me)
    hbo(me) = baseo(me) + ho(me)
    hmax(me) = 0.0d0  
    um(me)  = 0.0d0
    umo(me) = 0.0d0
    uum(me) = 0.0d0    
    vn(me)  = 0.0d0
    vno(me) = 0.0d0
    vvm(me) = 0.0d0
  enddo
  vrain = 0.0d0
  time  = 0
  mstep = 0
!--- Gradient ---  
  do me = 1, mesh
      call get_meid(me, 1, me1, tme, x1, y1)
      call get_meid(me, 2, me2, tme, x2, y2)
      call get_meid(me, 3, me3, tme, x3, y3)
      denominator = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1)
!--- Bed elevation --- 
      deltabx(me) = ((y3 - y1) * (baseo(me2) - baseo(me1)) - (y2 - y1) * (baseo(me3) - baseo(me1))) / denominator
      deltaby(me) = ((x2 - x1) * (baseo(me3) - baseo(me1)) - (x3 - x1) * (baseo(me2) - baseo(me1))) / denominator
    call calculate_wedelta(me1, me2, me3, deltabx, deltaby, wedeltabx(me), wedeltaby(me))      
  enddo
do me = 1, mesh
    do k = 1, ko(me)
            rx = ((dnox(linode(melink(me, k), 2))) +  (dnox(linode(melink(me, k), 1)))) * 0.5d0
            ry = ((dnoy(linode(melink(me, k), 2))) +  (dnoy(linode(melink(me, k), 1)))) * 0.5d0  
            reb(me, k) = baseo(me) + wedeltabx(me) * (rx - xmesh(me)) + wedeltaby(me) * (ry - ymesh(me))
            if (reb(me, k) < 0.0d0) reb(me, k) = 0.0d0          
    enddo
enddo  
!--- Observation point ---
do i = 1, op
    if (opz(i) > 0.1d0) then
        recz(i) = baseo(opmesh(i))
    else
        recz(i) = baseo(opmesh(i)) + wedeltabx(opmesh(i)) * (opx(i) - xmesh(opmesh(i))) + wedeltaby(opmesh(i)) * (opy(i) - ymesh(opmesh(i)))
    endif
enddo
  endsubroutine initial
  
!============================================================
!                          Predictor Step 
!============================================================   
  subroutine predictor
  use globals_condition
  use globals_meshdata
  use globals_ground
  use globals_kyokai
  !use globals_okyokai
  implicit none
  integer me, i, k, k2
  real(8) sumhpre, sumupre, sumvpre, sfx, dx1, dy1, dx, dy, nx, ny
  
  do me = 1, mesh
      sumhpre = 0.0d0
      sumupre = 0.0d0
      sumvpre = 0.0d0
      sumcpre = 0.0d0
      if (h(me) >= th) then
          do k = 1, ko(me)
              k2 = mod(k, ko(me)) + 1
              dx1 = dnox(menode(me, k2)) - dnox(menode(me, k))
              dy1 = dnoy(menode(me, k2)) - dnoy(menode(me, k))
              dx = dx1 / sqrt(dx1**2 + dy1**2)
              dy = dy1 / sqrt(dx1**2 + dy1**2)
              nx = dy !--- for test ---
              ny = -dx
!--- Flux Calculation ---
              sumhpre = sumhpre + (h(me)*uum(me)*nx + h(me)*vvm(me)*ny) * sqrt(dx1**2 + dy1**2)
              sumupre = sumupre + ((h(me)*uum(me)**2 + 0.5d0*gg*h(me)**2)*nx + h(me)*uum(me)*vvm(me)*ny) * sqrt(dx1**2 + dy1**2)
              sumvpre = sumvpre + (h(me)*uum(me)*vvm(me)*nx + (h(me)*vvm(me)**2 + 0.5d0*gg*h(me)**2)*ny) * sqrt(dx1**2 + dy1**2)
          enddo
          hpre(me) = h(me) - dt*sumhpre/smesh(me)
          upre(me) = uum(me) - (dt*sumupre/smesh(me) - dt*gg*h(me)*deltabx(me)) / h(me)
          vpre(me) = vvm(me) - (dt*sumvpre/smesh(me) - dt*gg*h(me)*deltaby(me)) / h(me)
          if (hpre(me) > th) then
              sfx = gg * (mn(me)**2) * (h(me)**(-1.333333)) * sqrt(uum(me)**2+vvm(me)**2) * dt
              ufpre(me) = 1.0d0 / (1.0d0 + sfx) * upre(me)
              vfpre(me) = 1.0d0 / (1.0d0 + sfx) * vpre(me)
          elseif (hpre(me) > 0.0d0 .and. hpre(me) < th) then
              ufpre(me) = 0.0d0
              vfpre(me) = 0.0d0
          elseif (hpre(me) < 0.0d0) then
              hpre(me) = 0.0d0
              ufpre(me) = 0.0d0
              vfpre(me) = 0.0d0
          endif
      else
          hpre(me) = h(me)
          ufpre(me) = 0.0d0
          vfpre(me) = 0.0d0
      endif
!--- replace ---
      h(me) = hpre(me)
      uum(me) = ufpre(me)
      vvm(me) = vfpre(me)
  enddo
  endsubroutine predictor  
  
!============================================================
!                           Linear Reconstruction 
!============================================================ 
  subroutine reconstruction
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer i, k, me, me1, me2, me3, tme
  real(8) rx, ry, x1, y1, x2, y2, x3, y3, denominator

!--- Gradient ---  
  do me = 1, mesh
      call get_meid(me, 1, me1, tme, x1, y1)
      call get_meid(me, 2, me2, tme, x2, y2)
      call get_meid(me, 3, me3, tme, x3, y3)
      denominator = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1)
!--- u ---            
      deltaux(me) = ((y3 - y1) * (uum(me2) - uum(me1)) - (y2 - y1) * (uum(me3) - uum(me1))) / denominator
      deltauy(me) = ((x2 - x1) * (uum(me3) - uum(me1)) - (x3 - x1) * (uum(me2) - uum(me1))) / denominator
!--- v ---            
      deltavx(me) = ((y3 - y1) * (vvm(me2) - vvm(me1)) - (y2 - y1) * (vvm(me3) - vvm(me1))) / denominator
      deltavy(me) = ((x2 - x1) * (vvm(me3) - vvm(me1)) - (x3 - x1) * (vvm(me2) - vvm(me1))) / denominator
  enddo
do me = 1, mesh
      call get_meid(me, 1, me1, tme, x1, y1)
      call get_meid(me, 2, me2, tme, x2, y2)
      call get_meid(me, 3, me3, tme, x3, y3)
!--- u ---    
    call calculate_wedelta1(me1, me2, me3, deltaux, deltauy, wedeltaux(me), wedeltauy(me))
!--- v ---
    call calculate_wedelta1(me1, me2, me3, deltavx, deltavy, wedeltavx(me), wedeltavy(me))
    do k = 1, ko(me)
        if (h(me) < th) then
            reu(me, k) = 0.0d0
            rev(me, k) = 0.0d0
        else
            rx = ((dnox(linode(melink(me, k), 2))) +  (dnox(linode(melink(me, k), 1)))) * 0.5d0
            ry = ((dnoy(linode(melink(me, k), 2))) +  (dnoy(linode(melink(me, k), 1)))) * 0.5d0  
            reu(me, k) = uum(me) + wedeltaux(me) * (rx - xmesh(me)) + wedeltauy(me) * (ry - ymesh(me))
            !if (reu(me, k) * uum(me) < 0.0d0) reu(me, k) = 0.0d0 !check
            rev(me, k) = vvm(me) + wedeltavx(me) * (rx - xmesh(me)) + wedeltavy(me) * (ry - ymesh(me))
            !if (rev(me, k) * vvm(me) < 0.0d0) rev(me, k) = 0.0d0 !check
        endif
    enddo         
enddo   
endsubroutine reconstruction

!============================================================
!                          Linear Reconstruction (h)
!============================================================ 
  subroutine SRM
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer li, me, tme, k1t, k2t, k, k2, k3, kk, meid
  real(8) xx, yy, bl, br, tb
  real(8) :: hl(mesh) 

  do li = 1, link
!--- during calculation, water depth less than th equals to zero ---
      if (h(limesh(li, 1)) <= th) then
          hl(limesh(li, 1)) = 0.0d0
          hb(limesh(li, 1)) = baseo(limesh(li, 1)) + 0.0d0
      else
          hl(limesh(li, 1)) = h(limesh(li, 1))
      endif
      if (h(limesh(li, 2)) <= th) then
          hl(limesh(li, 2)) = 0.0d0
          hb(limesh(li, 2)) = baseo(limesh(li, 2)) + 0.0d0
      else
          hl(limesh(li, 2)) = h(limesh(li, 2))
      endif
      do k1t = 1, ko(limesh(li, 1))
          if (melink(limesh(li, 1), k1t) == li) then
              k = k1t
              exit
          endif
      enddo
      k2 = mod(k, ko(limesh(li, 1))) + 1
      call get_meid(limesh(li, 1), k, meid, tme, xx, yy)
      !if (buildmark(li) == 1) then
      !    do k2t = 1, ko(tme)
      !        if (melink(tme, k2t) == li) then
      !            kk = k2t
      !            exit
      !        endif
      !    enddo
      !    bf(li) = baseo(tme) 
      !    reh(limesh(li, 1), k) = max(0.0d0, hl(limesh(li, 1)))
      !    reh(tme, kk) = max(0.0d0, hl(tme)) 
      !else
          do k2t = 1, ko(meid)
            if (melink(meid, k2t) == li) then
                kk = k2t
                exit
            endif
        enddo  
!--- H ---      
        db(li) = reb(meid, kk) - reb(limesh(li, 1), k)
        rehh(limesh(li, 1), k) = hb(limesh(li, 1)) + max(0.0d0, min((baseo(meid) - baseo(limesh(li, 1)) - db(li)), (hb(meid) - hb(limesh(li, 1)))))
        rehh(meid, kk) = hb(meid) + max(0.0d0, min((baseo(limesh(li, 1)) - baseo(meid) + db(li)), (hb(limesh(li, 1)) - hb(meid))))
!--- bL&bR ---
        bl = rehh(limesh(li, 1), k) - hl(limesh(li, 1))
        br = rehh(meid, kk) - hl(meid)
!--- bf ---
        bf(li) = max(bl, br)
!--- hL&hR ---
        reh(limesh(li, 1), k) = max(0.0d0, (rehh(limesh(li, 1), k) - bf(li)))
        reh(meid, kk) = max(0.0d0, (rehh(meid, kk) - bf(li)))
      !endif
 enddo
!--- bbar Check ---
    do me = 1, mesh
        do k = 1, ko(me)
            call get_meid(me, k, meid, tme, xx, yy)
            k2 = mod(k, ko(me)) + 1
            k3 = 0.5*ko(me)*(ko(me)+1) - k - k2
            if (bf(melink(me, k)) > rehh(me, k2) .or. bf(melink(me, k)) > rehh(me, k3)) then
                if (h(meid) < th) then
                    tb = max(0.0d0, (bf(melink(me, k)) - hb(me)))
                else
                    tb = max(0.0d0, min(abs(db(melink(me, k))), (bf(melink(me, k)) - hb(me))))
                endif
                bfk(me, k) = bf(melink(me, k)) - tb
            else
                bfk(me, k) = bf(melink(me, k))
            endif
        enddo
    enddo  
  endsubroutine SRM

!=== Obtain the MESHID and centroid coordinates of neighbor meshes ===
subroutine get_meid(me, index, meshid, trueid, x, y)
use globals_condition
use globals_meshdata
implicit none
integer :: me, index, meshid, trueid
real(8) x, y

if (limesh(melink(me, index), 1) == me) then
    if (buildmark(melink(me, index)) == 0 .and. limesh(melink(me, index), 2) /= 0) then
        meshid = limesh(melink(me, index), 2)
        trueid = limesh(melink(me, index), 2)
        x = xmesh(limesh(melink(me, index), 2))
        y = ymesh(limesh(melink(me, index), 2))
    elseif (buildmark(melink(me, index)) == 1 .or. limesh(melink(me, index), 2) == 0) then
        meshid = limesh(melink(me, index), 1)
        trueid = limesh(melink(me, index), 2)  !include 0 possible
        x = 0.5 * (dnox(linode(melink(me, index), 1)) + dnox(linode(melink(me, index), 2)))
        y = 0.5 * (dnoy(linode(melink(me, index), 1)) + dnoy(linode(melink(me, index), 2)))
    endif
elseif (limesh(melink(me, index), 2) == me) then
    if (buildmark(melink(me, index)) == 0 .and. limesh(melink(me, index), 2) /= 0) then
        meshid = limesh(melink(me, index), 1)
        trueid = limesh(melink(me, index), 1)
        x = xmesh(limesh(melink(me, index), 1))
        y = ymesh(limesh(melink(me, index), 1))
    elseif (buildmark(melink(me, index)) == 1 .or. limesh(melink(me, index), 2) == 0) then
        meshid = limesh(melink(me, index), 2)
        trueid = limesh(melink(me, index), 1)
        x = 0.5 * (dnox(linode(melink(me, index), 1)) + dnox(linode(melink(me, index), 2)))
        y = 0.5 * (dnoy(linode(melink(me, index), 1)) + dnoy(linode(melink(me, index), 2)))
    endif
endif
endsubroutine get_meid     

!=== Calculate wedelta of u, v ===
subroutine calculate_wedelta(me1, me2, me3, deltax, deltay, wedeltax, wedeltay)
    implicit none
    integer, intent(in) :: me1, me2, me3
    real(8), intent(in) :: deltax(:), deltay(:)
    real(8), intent(out) :: wedeltax, wedeltay
    real(8) :: g1, g2, g3, sum_g, we1, we2, we3
    g1 = deltax(me1)**2 + deltay(me1)**2
    g2 = deltax(me2)**2 + deltay(me2)**2
    g3 = deltax(me3)**2 + deltay(me3)**2
    sum_g = g1**2 + g2**2 + g3**2 + 3.0d-14
    we1 = (g2 * g3 + 1.0d-14) / sum_g
    we2 = (g1 * g3 + 1.0d-14) / sum_g
    we3 = (g1 * g2 + 1.0d-14) / sum_g
    wedeltax = we1 * deltax(me1) + we2 * deltax(me2) + we3 * deltax(me3)
    wedeltay = we1 * deltay(me1) + we2 * deltay(me2) + we3 * deltay(me3)  
end subroutine calculate_wedelta

!============================================================
!                          Flux + Riemann Solver
!============================================================
  subroutine FluxRiemann
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer li, me, tme, k, kk, k2, k1t, k2t, meid
  real(8) ustar, hstar, sl, sr, sm, ulver, urver, ulpara, urpara
  real(8) E1, E2, E3, E1L, E1R, E2L, E2R, E1bar, E2bar, E3bar
  real(8) dx, dy, dx1, dy1, nx, ny, xx, yy
  
  do li = 1, link
  if(limesh(li, 2) == 0) goto 300 !Boundary condition
  if(h(limesh(li, 1)) <= th .and. h(limesh(li, 2)) <= th) goto 303 !th  
  !if(buildmark == 1) goto 302
  do k1t = 1, ko(limesh(li, 1))
      if (melink(limesh(li, 1), k1t) == li) then
          k = k1t
          exit
      endif
  enddo
  k2 = mod(k, ko(limesh(li, 1))) + 1
  call get_meid(limesh(li, 1), k, meid, tme, xx, yy)
  do k2t = 1, ko(meid)
      if (melink(meid, k2t) == li) then
          kk = k2t
          exit
      endif
  enddo
!--- outward normal vector ---
      dx1 = dnox(menode(limesh(li, 1), k2)) - dnox(menode(limesh(li, 1), k))
      dy1 = dnoy(menode(limesh(li, 1), k2)) - dnoy(menode(limesh(li, 1), k))
      dx = dx1 / sqrt(dx1**2 + dy1**2)
      dy = dy1 / sqrt(dx1**2 + dy1**2)
      nx = dy 
      ny = -dx
!--- U⊥ & U// ---
      ulver = reu(limesh(li, 1), k) * nx + rev(limesh(li, 1), k) * ny
      ulpara = - reu(limesh(li, 1), k) * ny + rev(limesh(li, 1), k) * nx
      urver = reu(meid, kk) * nx + rev(meid, kk) * ny
      urpara = - reu(meid, kk) * ny + rev(meid, kk) * nx       
      !ulver = uum(limesh(li, 1)) * nx + vvm(limesh(li, 1)) * ny
      !ulpara = - uum(limesh(li, 1)) * ny + vvm(limesh(li, 1)) * nx
      !urver = uum(meid) * nx + vvm(meid) * ny
      !urpara = - uum(meid) * ny + vvm(meid) * nx  
!--- u* & h* ---          
      ustar = 0.5d0 * (ulver + urver) + sqrt(gg * reh(limesh(li, 1), k)) - sqrt(gg * reh(meid, kk))
      hstar = ((ulver - urver + 2.0d0 * (sqrt(gg * reh(limesh(li, 1), k)) + sqrt(gg * reh(meid, kk))))**2) / (16.0d0 * gg)
!--- SL & SR & SM---  
      if (reh(limesh(li, 1), k) > 0.0d0) then
          sl = min(ulver - sqrt(gg * reh(limesh(li, 1), k)), ustar - sqrt(gg * hstar))
      else
          sl = urver - 2.0d0 * sqrt(gg * reh(meid, kk))
      endif
      if (reh(meid, kk) > 0.0d0) then
          sr = max(urver + sqrt(gg * reh(meid, kk)), ustar + sqrt(gg * hstar))
      else
          sr = ulver + 2.0d0 * sqrt(gg * reh(limesh(li, 1), k))
      endif
      sm = (sl * reh(meid, kk) * (urver - sr) - sr * reh(limesh(li, 1), k) * (ulver - sl)) / (reh(meid, kk) * (urver - sr) - reh(limesh(li, 1), k) * (ulver - sl))
!--- Divided into hu^2 and 1/2gh^2 ---
!--- E1L E1R & E2L E2R ---
      E1L = reh(limesh(li, 1), k) * ulver
      E1R = reh(meid, kk) * urver
      E2L = reh(limesh(li, 1), k) * ulver**2
      E2R = reh(meid, kk) * urver**2
!--- E1&E2&E3 by SL&SR ---
      if (sl >= 0.0d0) then
          E1bar = E1L
          E2bar = E2L
          E3bar = E1L * ulpara
      elseif (sr <= 0.0d0) then
          E1bar = E1R
          E2bar = E2R
          E3bar = E1R * urpara
      elseif (sl <= 0.0d0 .and. sm >= 0.0d0) then
          E1bar = (sr * E1L - sl * E1R + sr * sl * (reh(meid, kk) - reh(limesh(li, 1), k))) / (sr - sl)
          E2bar = (sr * E2L - sl * E2R + sr * sl * (reh(meid, kk) * urver - reh(limesh(li, 1), k) * ulver)) / (sr - sl)
          E3bar = E1bar * ulpara
      elseif (sm <= 0.0d0 .and. sr >= 0.0d0) then
          E1bar = (sr * E1L - sl * E1R + sr * sl * (reh(meid, kk) - reh(limesh(li, 1), k))) / (sr - sl)
          E2bar = (sr * E2L - sl * E2R + sr * sl * (reh(meid, kk) * urver - reh(limesh(li, 1), k) * ulver)) / (sr - sl) 
          E3bar = E1bar * urpara          
      endif
      E1 = E1bar * sqrt(dx1**2 + dy1**2)
      E2 = (nx * E2bar - ny * E3bar) * sqrt(dx1**2 + dy1**2)
      E3 = (ny * E2bar + nx * E3bar) * sqrt(dx1**2 + dy1**2)
      h11(li) = rbeta(li) * E1
      u11(li) = rbeta(li) * E2
      v11(li) = rbeta(li) * E3  
!--- 1/2gh^2 ---
      E1L = rbeta(li) * 0.0d0
      E1R = rbeta(li) * 0.0d0
      E2L = rbeta(li) * 0.0d0 + 0.5d0 * gg * reh(limesh(li, 1), k)**2
      E2R = rbeta(li) * 0.0d0 + 0.5d0 * gg * reh(meid, kk)**2
       if (sl >= 0.0d0) then
          E1bar = E1L
          E2bar = E2L
          E3bar = E1L * ulpara
      elseif (sr <= 0.0d0) then
          E1bar = E1R
          E2bar = E2R
          E3bar = E1R * urpara
      elseif (sl <= 0.0d0 .and. sm >= 0.0d0) then
          E1bar = (sr * E1L - sl * E1R + sr * sl * (0.0d0)) / (sr - sl)
          E2bar = (sr * E2L - sl * E2R + sr * sl * (0.0d0)) / (sr - sl)
          E3bar = E1bar * ulpara
      elseif (sm <= 0.0d0 .and. sr >= 0.0d0) then
          E1bar = (sr * E1L - sl * E1R + sr * sl * (0.0d0)) / (sr - sl)
          E2bar = (sr * E2L - sl * E2R + sr * sl * (0.0d0)) / (sr - sl) 
          E3bar = E1bar * urpara          
      endif
      E1 = E1bar * sqrt(dx1**2 + dy1**2)
      E2 = (nx * E2bar - ny * E3bar) * sqrt(dx1**2 + dy1**2)
      E3 = (ny * E2bar + nx * E3bar) * sqrt(dx1**2 + dy1**2)
      u11g(li) = E2
      v11g(li) = E3
      goto 301
! --- slip boundary ---    
300    do k1t = 1, ko(limesh(li, 1))
                if (melink(limesh(li, 1), k1t) == li) then
                    k = k1t
                    exit
                endif
       enddo
      k2 = mod(k, ko(limesh(li, 1))) + 1
!--- outward normal vector ---
      dx1 = dnox(menode(limesh(li, 1), k2)) - dnox(menode(limesh(li, 1), k))
      dy1 = dnoy(menode(limesh(li, 1), k2)) - dnoy(menode(limesh(li, 1), k))
      dx = dx1 / sqrt(dx1**2 + dy1**2)
      dy = dy1 / sqrt(dx1**2 + dy1**2)
      nx = dy 
      ny = -dx
      h11(li) = 0.0d0
      u11(li) = 0.0d0
      v11(li) = 0.0d0 
      u11g(li) = 0.5d0 * gg * reh(limesh(li, 1), k)**2 * nx * sqrt(dx1**2 + dy1**2)
      v11g(li) = 0.5d0 * gg * reh(limesh(li, 1), k)**2 * ny * sqrt(dx1**2 + dy1**2)  
        goto 301   
303     h11(li) = 0.0d0
        u11(li) = 0.0d0
        v11(li) = 0.0d0
        u11g(li) = 0.0d0
        v11g(li) = 0.0d0
301 enddo   
  endsubroutine FluxRiemann

!====================================
!         H & HU
!====================================
  subroutine correctorstep
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer me, meid, li, k, k2, tme
  integer flgc, iflg
  real(8) sqx, rt, dx1, dy1, dx, dy, nx, ny
  real(8) :: ht(mesh)

  flgc = 0
  do li = 1, link
      corrh(li) = 1.0d0
      corru(li) = 1.0d0
      corrv(li) = 1.0d0
      corrc(li) = 1.0d0  
  enddo
  ht = h
403 iflg = 0
    flgc = flgc + 1
  do me = 1, mesh
    sumh(me) = 0.0d0
    sumu(me) = 0.0d0
    sumv(me) = 0.0d0
    u13h(me) = 0.0d0
    v13h(me) = 0.0d0
    do k = 1, ko(me)
        k2 = mod(k, ko(me)) + 1
        dx1 = dnox(menode(me, k2)) - dnox(menode(me, k))
        dy1 = dnoy(menode(me, k2)) - dnoy(menode(me, k))
        dx = dx1 / sqrt(dx1**2 + dy1**2)
        dy = dy1 / sqrt(dx1**2 + dy1**2)
        nx = dy
        ny = -dx
!--- bed slope term calculation ---
        !if (buildmark(melink(me, k)) == 1) then    
        !    u13h(me) = u13h(me) + 0.0d0
        !    v13h(me) = v13h(me) + 0.0d0
        !else
            u13h(me) = u13h(me) + (ht(me) + reh(me, k)) * (baseo(me) - bfk(me, k)) * nx * sqrt(dx1**2 + dy1**2) * lambda(me)
            v13h(me) = v13h(me) + (ht(me) + reh(me, k)) * (baseo(me) - bfk(me, k)) * ny * sqrt(dx1**2 + dy1**2) * lambda(me)
        !endif

      if (limesh(melink(me, k), 1) == me) then
          !if (buildmark(melink(me, k)) == 1 .and. inf(me) == 2) then
          !    sumh(me) = sumh(me) + h11(melink(me, k)) * corrh(melink(me, k))
          !    sumu(me) = sumu(me) - u22(melink(me, k)) * corru(melink(me, k))
          !    sumv(me) = sumv(me) - v22(melink(me, k)) * corrv(melink(me, k))
          !elseif (buildmark(melink(me, k)) == 1 .and. inf(me) /= 2) then
          !    sumh(me) = sumh(me) + h11(melink(me, k)) * corrh(melink(me, k))
          !    sumu(me) = sumu(me) + u11(melink(me, k)) * corru(melink(me, k))
          !    sumv(me) = sumv(me) + v11(melink(me, k)) * corrv(melink(me, k))
          !elseif(buildmark(melink(me, k)) /= 1) then !---
              sumh(me) = sumh(me) + h11(melink(me, k)) * corrh(melink(me, k))
              sumu(me) = sumu(me) + u11(melink(me, k)) * corru(melink(me, k)) + u11g(melink(me, k))* lambda(me)
              sumv(me) = sumv(me) + v11(melink(me, k)) * corrv(melink(me, k)) + v11g(melink(me, k))* lambda(me)
          !endif    
      else
          !if (buildmark(melink(me, k)) == 1 .and. inf(me) == 2) then
          !    sumh(me) = sumh(me) - h11(melink(me, k)) * corrh(melink(me, k))
          !    sumu(me) = sumu(me) - u22(melink(me, k)) * corru(melink(me, k))
          !    sumv(me) = sumv(me) - v22(melink(me, k)) * corrv(melink(me, k))
          !elseif (buildmark(melink(me, k)) == 1 .and. inf(me) /= 2) then
          !    sumh(me) = sumh(me) - h11(melink(me, k)) * corrh(melink(me, k))
          !    sumu(me) = sumu(me) + u11(melink(me, k)) * corru(melink(me, k))
          !    sumv(me) = sumv(me) + v11(melink(me, k)) * corrv(melink(me, k))
          !elseif(buildmark(melink(me, k)) /= 1) then
              sumh(me) = sumh(me) - h11(melink(me, k)) * corrh(melink(me, k))
              sumu(me) = sumu(me) - u11(melink(me, k)) * corru(melink(me, k)) - u11g(melink(me, k))* lambda(me)
              sumv(me) = sumv(me) - v11(melink(me, k)) * corrv(melink(me, k)) - v11g(melink(me, k))* lambda(me)
          !endif    
      endif
    enddo
  enddo 
    do me = 1, mesh
        h(me) = ho(me) - dt2*(sumh(me)/smesh(me)/lambda(me))
        um(me) = umo(me) - dt2*(sumu(me)/smesh(me)/lambda(me)) + dt2*0.5d0*gg*u13h(me)/smesh(me)/lambda(me)
        vn(me) = vno(me) - dt2*(sumv(me)/smesh(me)/lambda(me)) + dt2*0.5d0*gg*v13h(me)/smesh(me)/lambda(me)
!--- velocity ---
        if (h(me) < th) then
            uum(me) = 0.0d0
            vvm(me) = 0.0d0
        else
            uum(me) = um(me) / h(me)
            vvm(me) = vn(me) / h(me)
            sqx = lambda(me)*gg*(mn(me)**2)*sqrt(uum(me)**2+vvm(me)**2)*dt2/h(me)**(1.333333) + (1-lambda(me))*0.5d0*cd(me)*aj(me)*sqrt(uum(me)**2+vvm(me)**2)*dt2
!--- update ---
            uum(me) = lambda(me)*uum(me)/(lambda(me)+sqx)
            vvm(me) = lambda(me)*vvm(me)/(lambda(me)+sqx)
        endif
!--- minus water depth corrector ---        
        if (h(me) < 0.0d0) then
            iflg = 1
            rt = (ho(me) - th) / (ho(me) - h(me))
            do k = 1, ko(me)
                corrh(melink(me, k)) = min(corrh(melink(me, k)), rt)
                corru(melink(me, k)) = 0.0d0
                corrv(melink(me, k)) = 0.0d0
            enddo
            h(me) = 0.0d0
            um(me) = 0.0d0
            vn(me) = 0.0d0
            uum(me) = 0.0d0
            vvm(me) = 0.0d0
        endif
    enddo
    if (iflg ==1 .and. flgc <= 5) goto 403
  endsubroutine correctorstep

!============================================================
!                   Replace Old >> New
!============================================================
!=== GROUND VARIABLE ================
  subroutine replace_ground
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer me, i
  real(8) sum, rmse
  
  do me = 1, mesh
    umo(me) = h(me) * uum(me)
    vno(me) = h(me) * vvm(me)
    ho(me) = h(me)
    hbo(me) = hb(me)
    hmax(me) = max(hmax(me), h(me))
  enddo
!--- rmse output -----
    if (time < 97200.0d0 .and. abs(time-97200.0d0) < (1.05*dt)) then 
        sum = 0
        do i = 1, op
            sum = sum + (hmax(opmesh(i)) - fm(i))**2
        enddo
        rmse = sqrt(sum/op)
        write(12, *) rmse
    endif
  endsubroutine replace_ground
  
!============================================================
!                       From Outside
!============================================================
!====================================
!             Rainfall
!====================================
  subroutine rainfall
  use globals_condition
  use globals_meshdata
  use globals_rain
  implicit none
  integer me, it
  
  do me = 1, mesh
    it = int(time/dtrain) + 1
    rr(me) = rain(it) * rnof(me)
    if(inf(me) /= 0) then
      vrain = vrain + rr(me)/1000.d0/dtrain*dt2*smesh(me)*(1.0d0 - lambda(me))
    endif
  enddo
  endsubroutine rainfall
!====================================
!            Qin Kyokai
!====================================
  subroutine lkyokai
  use globals_condition
  use globals_meshdata
  use globals_ground
  use globals_kyokai
  implicit none
  integer ii
  real(8) lqin1, lqin2
  
  ii = int(time/dtq) + 1
  lqin1 = qin1(ii) + (time - dtq*dble(ii - 1))/dtq*(qin1(ii + 1) - qin1(ii)) 
  lqin2 = qin2(ii) + (time - dtq*dble(ii - 1))/dtq*(qin2(ii + 1) - qin2(ii)) 
  h11(inl1) = -lqin1
  h11(inl2) = -lqin2
  vkyokai = vkyokai + -h11(inl1)*dt2 + -h11(inl2)*dt2
  endsubroutine lkyokai
!====================================
  subroutine okyokai
  use globals_condition
  use globals_meshdata
  use globals_ground
  use globals_okyokai
  implicit none
  integer i, k, ok, ok2
  real(8) hstar, ustar, vstar, uver
  real(8) dx1, dy1, dx, dy, nx, ny
  
  do i = 1, oid
      do k = 1, ko(limesh(outl(i), 1))
          if (melink(limesh(outl(i), 1), k) == outl(i)) then
              ok = k
              exit
          endif
      enddo
      hstar = reh(limesh(outl(i), 1), ok)
      ustar = reu(limesh(outl(i), 1), ok)
      vstar = rev(limesh(outl(i), 1), ok)
      !hstar = h(limesh(outl(i), 1))
      !ustar = uum(limesh(outl(i), 1))
      !vstar = vvm(limesh(outl(i), 1))
      ok2 = mod(ok, ko(limesh(outl(i), 1))) + 1
!--- outward normal vector ---
      dx1 = dnox(menode(limesh(outl(i), 1), ok2)) - dnox(menode(limesh(outl(i), 1), ok))
      dy1 = dnoy(menode(limesh(outl(i), 1), ok2)) - dnoy(menode(limesh(outl(i), 1), ok))
      dx = dx1 / sqrt(dx1**2 + dy1**2)
      dy = dy1 / sqrt(dx1**2 + dy1**2)
      nx = dy 
      ny = -dx
!--- U⊥ ----
      uver = ustar * nx + vstar * ny
      !uver = 0.5d0 * nx + 0.0d0 * ny
      !h11(outl(i)) = 
      !u11(outl(i)) = 
      !v11(outl(i)) = 
  enddo
  endsubroutine okyokai

!============================================================
!                   Water Volume Calculate
!============================================================
  subroutine sumqa(sv,sa,saj,vin,vout,in_out)
  use globals_condition
  use globals_meshdata
  use globals_ground
  use globals_kyokai
  implicit none
  real(8) sv, sa, saj, sumr, suml
  real(8) vin, vout, in_out
  integer me
  
  sv = 0.0d0
  sa = 0.0d0
  saj = 0.0d0
  do me = 1, mesh
    if(inf(me) /= 0) then
      sv = sv + h(me)*smesh(me)*lambda(me)
      if(h(me) > 0.0d0) sa  = sa  + smesh(me)*lambda(me)
      if(h(me) > 0.0d0) saj = saj + smesh(me)*(1.0d0 - lambda(me)) 
    endif
  enddo
  sumr  = 0.0d0
  suml  = 0.0d0
  vin   = 0.0d0
  vout  = 0.0d0
  in_out = 0.0d0
  if(type_r  == 1) call v_rain(sumr)
  if(type_fl == 1) call v_kyokai(suml)  
  vin  = sumr + suml
  vout = sv
  in_out = vin - vout
  endsubroutine sumqa
!=== RAINFALL VOLUME ================
  subroutine v_rain(sumr)
  use globals_rain
  implicit none
  real(8) sumr
  sumr = vrain
  endsubroutine v_rain
!=== KYOKAI VOLUME ==================
  subroutine v_kyokai(suml)
  use globals_kyokai
  implicit none
  real(8) suml
  suml = vkyokai
  endsubroutine v_kyokai
!============================================================
!                    Results Writing
!============================================================
!=== DISPLAY ========================
  subroutine dispwrite
  use globals_condition
  implicit none
  real(8) sv, sa, saj, vin, vout, in_out
  call sumqa(sv,sa,saj,vin,vout,in_out)

  write(*,'("TIME = ",f7.1,"     SV = ",f20.3,"     SA = ",f15.3, &
      & "     Vin = ",f20.4,"     Vout = ",f20.4,"     in-out = ",f20.4)') time, sv, sa, vin, vout, in_out
  write(104,'("TIME = ",f7.1,"     SA = ",f15.3, &
      & "     Vin = ",f20.3,"     Vout = ",f20.3,"     in-out = ",f20.3)') time, sa, vin, vout, in_out
  write(101,'(f20.3,"   = SV","   TIME = ",f7.1)') sv, time
  end subroutine dispwrite
!=== OUTPUT H =======================
  subroutine diskwrite_h
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer me
  write(90, '("TIME = ", f10.1)') time
  write(90, '(10(f15.7,2x))') (h(me), me = 1, mesh)
  end subroutine diskwrite_h
!=== OUTPUT HMAX ====================
  subroutine diskwrite_hmax
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer me, i
  write(91, '("TIME = ", f10.1)') time
  write(91, '(10(f15.7,2x))') (hmax(me), me = 1, mesh)
  write(203, *) 'hmaxopmesh'
  write(203, '(*(f15.7))') (hmax(opmesh(i)), i = 1, op)
  end subroutine diskwrite_hmax
!=== OUTPUT TECPLOT h ===============
subroutine tecwrite_h
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer no, me, nvar
  nvar = 0
  if(type_tech   == 1) nvar = 0
  if(type_tecbs  == 1) nvar = nvar + 1
  if(type_tecinf == 1) nvar = nvar + 2
  if(time == 0) then
    if(nvar == 3) write(95, *) 'variables = "x", "y", "h(m)", "um(m/s)", "vn(m/s)", "baseo(m)", "inf(m)"'
    if(nvar == 1) write(95, *) 'variables = "x", "y", "h(m)", "baseo(m)"'
    if(nvar == 2) write(95, *) 'variables = "x", "y", "h(m)", "inf(m)"'
    if(nvar == 0) write(95, *) 'variables = "x", "y", "h(m)"'
    write(95, 991) time
991 format('zone t = "',f10.0 ,'"')
    write(95, '("solutiontime = ",f7.1)') time
    write(95, *) 'datapacking=block'
    if(nvar == 3) write(95, *) 'varlocation=([3-7]=cellcentered)'
    if(nvar == 1) write(95, *) 'varlocation=([3-4]=cellcentered)'
    if(nvar == 2) write(95, *) 'varlocation=([3-4]=cellcentered)'
    if(nvar == 0) write(95, *) 'varlocation=([3]=cellcentered)'
    write(95, '("nodes =", i10)')    node
    write(95, '("elements =", i10)') mesh
    write(95, *) 'zonetype=fetriangle'
    do no = 1, node
      write(95, *) dnox(no)
    enddo
    do no = 1, node
      write(95, *) dnoy(no)
    enddo
    do me = 1, mesh
      write(95, *) h(me)
    enddo
!=== Velocity output ===   
    do me = 1, mesh
      write(95, *) uum(me)
    enddo
    do me = 1, mesh
      write(95, *) vvm(me)
    enddo
    if(type_tecbs   == 1) then
      do me = 1, mesh
        write(95, *) baseo(me)
      enddo
    endif
    if(type_tecinf  == 1) then
      do me = 1, mesh
        write(95, *) inf(me)
      enddo
    endif
    do me = 1, mesh
      write(95, *) menode(me,1), menode(me,2), menode(me,3)
    enddo
  endif
  if(time /= 0)then
    if(nvar == 3) write(95, *) 'variables = "x", "y", "h(m)", "um(m/s)", "vn(m/s)", "baseo(m)", "inf(m)"'
    if(nvar == 2) write(95, *) 'variables = "x", "y", "h(m)", "baseo(m)"'
    if(nvar == 1) write(95, *) 'variables = "x", "y", "h(m)", "inf(m)"'
    if(nvar == 0) write(95, *) 'variables = "x", "y", "h(m)"'
    write(95, 992) time
992 format('zone t = "',f10.0 ,'"')
    write(95, '("solutiontime = ",f7.1)') time
    write(95, *) 'datapacking=block'
    if(nvar == 3) write(95, *) 'varlocation=([3-7]=cellcentered)'
    if(nvar == 1) write(95, *) 'varlocation=([3-4]=cellcentered)'
    if(nvar == 2) write(95, *) 'varlocation=([3-4]=cellcentered)'
    if(nvar == 0) write(95, *) 'varlocation=([3]=cellcentered)'
    write(95, '("nodes =", i10)')    node
    write(95, '("elements =", i10)') mesh
    write(95, *) 'zonetype=fetriangle'
    if(nvar == 3) write(95, *) 'Varsharelist=([1,2,6,7]=1)'
    if(nvar == 2) write(95, *) 'Varsharelist=([1,2,4]=1)'
    if(nvar == 1) write(95, *) 'Varsharelist=([1,2,4]=1)'
    if(nvar == 0) write(95, *) 'Varsharelist=([1,2]=1)'
    write(95, *) 'connectivitysharezone=1'
    do me = 1, mesh
      write(95, *) h(me)
    enddo 
    !=== Velocity output ===   
    do me = 1, mesh
      write(95, *) uum(me)
    enddo
    do me = 1, mesh
      write(95, *) vvm(me)
    enddo
  endif
  end subroutine tecwrite_h  
  
!=== OUTPUT TECPLOT hmax ============
  subroutine tecwrite_hmax
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer no, me, nvar
  nvar = 0
  if(type_techm  == 1) nvar = 0
  if(type_tecbs  == 1) nvar = nvar + 1
  if(type_tecinf == 1) nvar = nvar + 2
    if(nvar == 3) write(96, *) 'variables = "x", "y", "hmax(m)", "baseo(m)", "inf(m)"'
    if(nvar == 2) write(96, *) 'variables = "x", "y", "hmax(m)", "baseo(m)"'
    if(nvar == 1) write(96, *) 'variables = "x", "y", "hmax(m)", "inf(m)"'
    if(nvar == 0) write(96, *) 'variables = "x", "y", "hmax(m)"'
  write(96, *) 'zone t= "geo data"'
  write(96, *) 'datapacking=block'
    if(nvar == 3) write(96, *) 'varlocation=([3-5]=cellcentered)'
    if(nvar == 1) write(96, *) 'varlocation=([3-4]=cellcentered)'
    if(nvar == 2) write(96, *) 'varlocation=([3-4]=cellcentered)'
    if(nvar == 0) write(96, *) 'varlocation=([3]=cellcentered)'
  write(96, '("nodes =", i10)')    node
  write(96, '("elements =", i10)') mesh
  write(96, *) 'zonetype=fetriangle'
  do no = 1, node
    write(96, *) dnox(no)
  enddo
  do no = 1, node
    write(96, *) dnoy(no)
  enddo
  if(type_techm == 1) then
    do me = 1, mesh
      write(96, *) hmax(me)
    enddo
  endif
  if(type_tecbs == 1) then
    do me = 1, mesh
      write(96, *) baseo(me)
    enddo
  endif
  if(type_tecinf == 1) then
    do me = 1, mesh
      write(96, *) inf(me)
    enddo
  endif
  do me = 1, mesh
    write(96, *) menode(me,1), menode(me,2), menode(me,3)
  enddo
  end subroutine tecwrite_hmax
!=== OUTPUT TECPLOT bs ============
  subroutine tecwrite_bs
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer no, me
  write(97, *) 'variables = "x", "y", "baseo(m)"'
  write(97, *) 'zone t= "geo data"'
  write(97, *) 'datapacking=block'
  write(97, *) 'varlocation=([3]=cellcentered)'
  write(97, '("nodes =", i10)')    node
  write(97, '("elements =", i10)') mesh
  write(97, *) 'zonetype=fetriangle'
  do no = 1, node
    write(97, *) dnox(no)
  enddo
  do no = 1, node
    write(97, *) dnoy(no)
  enddo
  do me = 1, mesh
    write(97, *) baseo(me)
  enddo
  do me = 1, mesh
    write(97, *) menode(me,1), menode(me,2), menode(me,3)
  enddo
  end subroutine tecwrite_bs
!=== OUTPUT TECPLOT inf ============
  subroutine tecwrite_inf
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer no, me
  write(98, *) 'variables = "x", "y", "inf(m)"'
  write(98, *) 'zone t= "geo data"'
  write(98, *) 'datapacking=block'
  write(98, *) 'varlocation=([3]=cellcentered)'
  write(98, '("nodes =", i10)')    node
  write(98, '("elements =", i10)') mesh
  write(98, *) 'zonetype=fetriangle'
  do no = 1, node
    write(98, *) dnox(no)
  enddo
  do no = 1, node
    write(98, *) dnoy(no)
  enddo
  do me = 1, mesh
    write(98, *) inf(me)
  enddo
  do me = 1, mesh
    write(98, *) menode(me,1), menode(me,2), menode(me,3)
  enddo
  end subroutine tecwrite_inf
  
  !=== OUTPUT TECPLOT endh ===
  subroutine tecwrite_endh
  use globals_condition
  use globals_meshdata
  use globals_ground
  implicit none
  integer no, me
  write(102, *) 'variables = "x", "y", "h(m)"'
  write(102, *) 'zone t= "geo data"'
  write(102, *) 'datapacking=block'
  write(102, *) 'varlocation=([3]=cellcentered)'
  write(102, '("nodes =", i10)')    node
  write(102, '("elements =", i10)') mesh
  write(102, *) 'zonetype=fetriangle'
  do no = 1, node
    write(102, *) dnox(no)
  enddo
  do no = 1, node
    write(102, *) dnoy(no)
  enddo
  do me = 1, mesh
    write(102, *) h(me)
  enddo
  do me = 1, mesh
    write(102, *) menode(me,1), menode(me,2), menode(me,3)
  enddo
  end subroutine tecwrite_endh

    !=== OUTPUT endh ===
  subroutine diskwrite_endh
    use globals_condition
    use globals_meshdata
    use globals_ground
    implicit none
    integer me
    write(103, '("TIME = ", f10.1)') time
    do me = 1, mesh
      write(103, *) h(me)
    enddo
    end subroutine diskwrite_endh



  end module

!     #######################################################
!     ##                                                                                                      ##
!     ##                   main program                                                            ##
!     ##                                                                                                      ##
!     #######################################################
!
program main
  use globals_condition
  use globals_ground
  use subprogs
implicit none
  character(40) fnode, finf, fbase, fmesh, flink, fmesho
  character(40) frain, ffl
  character(40) fh, fhmax, fth, fthm, ftbs, ftinf
  character(40) fsv, ftendh, fendh, fvsearch
  character(len=256) :: input_filename, output_filename
  integer i, status, ierr
  integer lpout, lkout
  real(8) start, finish
  
  call GET_COMMAND_ARGUMENT(1, input_filename)
  call GET_COMMAND_ARGUMENT(2, output_filename)
  
  write(*, *) ' === PROGRAM START !! === '
  write(*, *)
  
!============================================================
!                          READ FILE
!============================================================
  open(10, file = 'urbanfloodmodel/readfile/condition.dat', action = 'read')
 !=== Check Files ======================== 
  open(11, file = trim(input_filename), action = 'read', iostat=ierr)
  if (ierr /= 0) then
    write(*,*) 'Error opening dynamic input file: ', trim(input_filename)
    stop 1
  endif
  open(12, file = trim(output_filename), action = 'write', iostat=ierr)
  if (ierr /= 0) then
    write(*,*) 'Error opening dynamic output file: ', trim(output_filename)
    stop 1
  endif
  open(201, file = 'urbanfloodmodel/readfile/data/2D/porosity.dat', action = 'read')
  open(202, file = 'urbanfloodmodel/readfile/data/2D/op.dat', action = 'read')
  open(203, file = 'urbanfloodmodel/out/hmaxopmesh.dat', action = 'write')
  
!====================================
!       Calcultion Condition
!====================================
!=== CALCULATION TYPE ===============
  read(10, *)
  read(10, *)
  read(10, *)
  read(10, '(i1)') type_r
  read(10, '(i1)') type_fl
  read(10, '(i1)') type_tech
  read(10, '(i1)') type_techm
  read(10, '(i1)') type_tecbs
  read(10, '(i1)') type_tecinf
!=== CAL TIME & OUTPUT TIME =========
  read(10, *)
  read(10, *)
  read(10, *)
  read(10, *) timmax
  read(10, *) dt
  read(10, *) dkout
  read(10, *) dpout
!=== MESH INFORMATION ===============
  read(10, *)
  read(10, *)
  read(10, *) num_inf
  allocate(n_inf(num_inf), inf_mn(num_inf), inf_cd(num_inf))
  read(10, *)
  do i = 1, num_inf
    read(10, *) n_inf(i)
  enddo
!====================================
!          Data File Name
!====================================
!=== INPUT FILE NAME ================
  do i = 1, 4
    read(10, *)
  enddo
  read(10, '(a40)') fnode
  read(10, '(a40)') finf
  read(10, '(a40)') fbase
  read(10, *)
  read(10, *) status
  read(10, *)
  read(10, '(a40)') fmesh
  read(10, '(a40)') flink
  read(10, *)
  if(status == 1) then
    read(10, '(a40)') fmesho
  else
    read(10, *)
  endif
  read(10, *)
  if(type_r == 1) then
    read(10, '(a40)') frain
  else
    read(10, *)
  endif
  if(type_fl == 1) then
    read(10, '(a40)') ffl
  else
    read(10, *)
  endif
  
  do i = 1, 4
    read(10, *)
  enddo
  read(10, '(a40)') fh
  read(10, '(a40)') fhmax
  read(10, '(a40)') fth
  read(10, '(a40)') fthm
  read(10, '(a40)') ftbs
  read(10, '(a40)') ftinf
  read(10, '(a40)') fsv
  read(10, '(a40)') ftendh
  read(10, '(a40)') fendh
  read(10, '(a40)') fvsearch
!!=== OUTPUT FILE NAME ===============
!  write(*, *) '   -=-=- FINISH READING CONDITION FILE -=-=-'
!  write(*, *)
  !=== CONSTANT =======================
  dt2  = 2.0d0 * dt
  gg   = 9.8d0
  fita = 0.5d0
  pi   = 4.0d0 * atan(1.0d0)
!=== WRITING TIME ===================
  lpout = int(dpout/dt)
  lkout = int(dkout/dt)
!============================================================
!                      Preparation
!============================================================
!====================================
!          Open Data File
!====================================
  open(50, file = fnode, action = 'read')
  open(51, file = finf,  action = 'read')
  open(52, file = fbase, action = 'read')
  open(53, file = fmesh)
  open(54, file = flink)
  if(status == 1) open(55, file = fmesho , action = 'read')
  if(type_r == 1) open(56, file = frain  , action = 'read')
  if(type_fl == 1) open(57, file = ffl   , action = 'read')
  open(90, file = fh, action = 'write')
  open(91, file = fhmax, action = 'write')
  if(type_tech   == 1) open(95, file = fth, action = 'write')
  if(type_techm  == 1) open(96, file = fthm, action = 'write')
  if(type_tecbs  == 1) open(97, file = ftbs, action = 'write')
  if(type_tecinf == 1) open(98, file = ftinf, action = 'write')
  open(101, file = fsv, action = 'write')
  open(102, file = ftendh, action = 'write')
  open(103, file = fendh, action = 'write')
  open(104, file = fvsearch, action = 'write')
  
!====================================
!         Various Settings
!====================================
!=== READ DATA ======================
  call rdat
  if(type_r  == 1) call rraindat
  if(type_fl == 1) call rkyokaidat
  !call rokyokaidat 
  call observationpoints
!=== ALLOCATE ARRAY SIZE ============
  call allocate_ground
  if(type_r  == 1) call allocate_rain
!=== PREPARE SEWERAGE DATA ==========
!=== INITIAL VALUE ==================
  call initial
!=== WRITE INITIAL VALUE ============
!  call diskwrite_h
!  !call dispwrite
!  if(type_tech == 1) then
!    call tecwrite_h
!  endif
  
  !write(*, *) '   -=-=- FINISH SETTING INITIAL CONDITION -=-=-'
  !write(*, *)
  
!============================================================
!                         Loop
!============================================================
!+++ LOOP START ! +++++++++++++++++++
  call CPU_TIME(start)
!====================================
!               Flux
!====================================
!1 call predictor
1  call SRM
!=== TIMESTEP =======================
  time = time + dt
  mstep = mstep + 1
!===================================  
  call reconstruction
  call FluxRiemann
  if(type_fl == 1) then
    call lkyokai
  endif
  !call okyokai 
!====================================
!              Suisin & Velocity
!====================================
  if(type_r == 1) call rainfall
  call correctorstep
!====================================
!             Replace
!====================================
  call replace_ground
!=== TIMESTEP =======================
  time = time + dt
  mstep = mstep + 1
!====================================
!             Output
!====================================
!  if(mod(mstep, lkout) == 0) then
!    call diskwrite_h
!    if(type_tech == 1) then
!      call tecwrite_h
!    endif
!  endif
!  if(mod(mstep, lpout) == 0) then
!    call dispwrite
!  endif

!!!!=== IF YOU WANT TEST-OUTPUT-DATA ! =
!    if(mod(mstep, lkout) == 0) call test_output
!====================================
!            Time Jugde
!====================================
  if(time + dt <= timmax) goto 1
!+++ LOOP END ! +++++++++++++++++++++
!============================================================
!                    All Results Output
!============================================================
!=== FILE OUTPUT ====================
  !call diskwrite_h
  !call diskwrite_hmax
!=== TECPLOT HMAX BS INF ============
  !if(type_techm == 1) then
  !  call tecwrite_hmax
  !endif
  !if(type_tecbs == 1) then
  !  call tecwrite_bs
  !endif
  !if(type_tecinf == 1) then
  !  call tecwrite_inf
  !endif
!!=== ENDH OUTPUT ===
!  call tecwrite_endh
!  call diskwrite_endh

call CPU_TIME(finish)
print *, 'CPU Time Elapsed = ', finish - start, ' seconds.'

write(*,*) ' === normal end === '

end program main