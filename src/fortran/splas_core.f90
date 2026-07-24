module splas_core
  use, intrinsic :: iso_fortran_env, only: real64, int64
  implicit none
  private

  integer, parameter, public :: nphy = 70
  integer, parameter, public :: nzoo = 70
  integer, parameter, public :: nnut = 3
  integer, parameter :: nparams = 26
  real(real64), parameter :: pi = acos(-1.0_real64)
  real(real64), parameter :: tiny_value = 1.0e-30_real64

  public :: run_model, model_traits

contains

  subroutine model_traits(phyto_sink_coeff, prey_tolerance, phyto_size, zoo_size, mu_max, &
                          vmax, qmin, half_sat, phyto_sink, carbon_cell, zoo_ingestion, &
                          optimal_prey, preference)
    double precision, intent(in) :: phyto_sink_coeff, prey_tolerance
    double precision, intent(out) :: phyto_size(nphy), zoo_size(nzoo)
    double precision, intent(out) :: mu_max(nphy), vmax(nphy), qmin(nphy), half_sat(nphy)
    double precision, intent(out) :: phyto_sink(nphy), carbon_cell(nphy)
    double precision, intent(out) :: zoo_ingestion(nzoo), optimal_prey(nzoo)
    double precision, intent(out) :: preference(nphy, nzoo)
    real(real64) :: volume(nphy)
    integer :: i, j

    call logspace(1.0_real64, 200.0_real64, phyto_size)
    call logspace(2.0_real64, 1000.0_real64, zoo_size)
    volume = (pi / 6.0_real64) * phyto_size**3
    mu_max = 2.6_real64 * phyto_size**(-0.45_real64)
    vmax = 9.10e-12_real64 * volume**0.67_real64
    qmin = 1.36e-12_real64 * volume**0.77_real64
    half_sat = 0.1_real64 * phyto_size
    phyto_sink = phyto_sink_coeff * phyto_size**1.17_real64
    carbon_cell = 0.0109_real64 * volume**0.991_real64
    zoo_ingestion = 26.0_real64 * zoo_size**(-0.4_real64)
    optimal_prey = 0.16_real64 * zoo_size * exp(-0.02_real64 * log(zoo_size)**2)
    do j = 1, nzoo
      do i = 1, nphy
        preference(i, j) = exp(-((log10(phyto_size(i)) - log10(optimal_prey(j))) / &
                                 prey_tolerance)**2)
      end do
    end do
  end subroutine model_traits


  subroutine run_model(n_steps, output_every, nout, m, dt, dz, system_mode, forcing_mode, &
                       profile_mode, seed, params, p_out, z_out, q_out, n_out, det_out, &
                       surface_out, irradiance_out, nutlim_out, grazlim_out, lightlim_out, &
                       growth_out, grazing_out, export_p_out, export_det_out, mass_error_out, &
                       correction_out, status)
    integer, intent(in) :: n_steps, output_every, nout, m
    integer, intent(in) :: system_mode, forcing_mode, profile_mode, seed
    double precision, intent(in) :: dt, dz, params(26)
    double precision, intent(out) :: p_out(nphy, nout, m), z_out(nzoo, nout, m)
    double precision, intent(out) :: q_out(nphy, nout, m), n_out(nnut, nout, m)
    double precision, intent(out) :: det_out(nout, m), surface_out(nout)
    double precision, intent(out) :: irradiance_out(nout, m)
    double precision, intent(out) :: nutlim_out(nphy, nout, m)
    double precision, intent(out) :: grazlim_out(nphy, nout, m)
    double precision, intent(out) :: lightlim_out(nphy, nout, m)
    double precision, intent(out) :: growth_out(nphy, nout, m)
    double precision, intent(out) :: grazing_out(nphy, nout, m)
    double precision, intent(out) :: export_p_out(nout), export_det_out(nout)
    double precision, intent(out) :: mass_error_out(nout), correction_out(nout)
    integer, intent(out) :: status

    real(real64) :: p(nphy, m), z(nzoo, m), q(nphy, m), nutrient(m), det(m)
    real(real64) :: dp1(nphy, m), dz1(nzoo, m), dq1(nphy, m), dn1(m), dd1(m)
    real(real64) :: dp2(nphy, m), dz2(nzoo, m), dq2(nphy, m), dn2(m), dd2(m)
    real(real64) :: pp(nphy, m), zp(nzoo, m), qp(nphy, m), npred(m), detp(m)
    real(real64) :: phyto_size(nphy), zoo_size(nzoo), mu_max(nphy), vmax(nphy)
    real(real64) :: qmin(nphy), half_sat(nphy), phyto_sink(nphy), carbon_cell(nphy)
    real(real64) :: zoo_ingestion(nzoo), optimal_prey(nzoo), preference(nphy, nzoo)
    real(real64) :: irr1(m), irr2(m), nutlim1(nphy, m), nutlim2(nphy, m)
    real(real64) :: grazlim1(nphy, m), grazlim2(nphy, m)
    real(real64) :: lightlim1(nphy, m), lightlim2(nphy, m)
    real(real64) :: growth1(nphy, m), growth2(nphy, m)
    real(real64) :: grazing1(nphy, m), grazing2(nphy, m)
    real(real64) :: surface, initial_mass, current_mass, correction
    integer(int64) :: rng_state
    integer :: step, output_index, k

    p_out = 0.0_real64
    z_out = 0.0_real64
    q_out = 0.0_real64
    n_out = 0.0_real64
    det_out = 0.0_real64
    surface_out = 0.0_real64
    irradiance_out = 0.0_real64
    nutlim_out = 0.0_real64
    grazlim_out = 0.0_real64
    lightlim_out = 0.0_real64
    growth_out = 0.0_real64
    grazing_out = 0.0_real64
    export_p_out = 0.0_real64
    export_det_out = 0.0_real64
    mass_error_out = 0.0_real64
    correction_out = 0.0_real64
    status = 0

    if (n_steps < 1 .or. output_every < 1 .or. m < 3) then
      status = 10
      return
    end if
    if (nout /= n_steps / output_every + 1 .or. mod(n_steps, output_every) /= 0) then
      status = 11
      return
    end if
    if (dt <= 0.0_real64 .or. dz <= 0.0_real64) then
      status = 12
      return
    end if

    call model_traits(params(7), params(20), phyto_size, zoo_size, mu_max, vmax, qmin, &
                      half_sat, phyto_sink, carbon_cell, zoo_ingestion, optimal_prey, preference)
    if (system_mode == 0) phyto_sink = 0.0_real64

    do k = 1, m
      p(:, k) = params(21)
      z(:, k) = params(22)
      q(:, k) = max(params(25), 1.0_real64) * qmin
      nutrient(k) = real(k, real64) * params(23)
      det(k) = params(24)
    end do
    correction = 0.0_real64
    rng_state = max(1_int64, int(abs(seed), int64))
    initial_mass = inventory(nutrient, p, z, det, dz)

    surface = forcing_value(0.0_real64, n_steps * dt, forcing_mode, params, rng_state)
    call rhs(m, dz, system_mode, profile_mode, surface, params, mu_max, vmax, qmin, half_sat, &
             phyto_sink, zoo_ingestion, preference, p, z, q, nutrient, det, dp1, dz1, dq1, &
             dn1, dd1, irr1, nutlim1, grazlim1, lightlim1, growth1, grazing1)
    output_index = 1
    call save_output(output_index, m, dz, surface, params, phyto_sink, carbon_cell, p, z, q, &
                     nutrient, det, irr1, nutlim1, grazlim1, lightlim1, growth1, grazing1, &
                     initial_mass, correction, p_out, z_out, q_out, n_out, det_out, surface_out, &
                     irradiance_out, nutlim_out, grazlim_out, lightlim_out, growth_out, grazing_out, &
                     export_p_out, export_det_out, mass_error_out, correction_out)

    do step = 1, n_steps
      surface = forcing_value(real(step, real64) * dt, n_steps * dt, forcing_mode, params, rng_state)
      call rhs(m, dz, system_mode, profile_mode, surface, params, mu_max, vmax, qmin, half_sat, &
               phyto_sink, zoo_ingestion, preference, p, z, q, nutrient, det, dp1, dz1, dq1, &
               dn1, dd1, irr1, nutlim1, grazlim1, lightlim1, growth1, grazing1)

      pp = p + dt * dp1
      zp = z + dt * dz1
      qp = q + dt * dq1
      npred = nutrient + dt * dn1
      detp = det + dt * dd1
      call enforce_state(m, qmin, pp, zp, qp, npred, detp, dz, correction, status)
      if (status /= 0) return

      call rhs(m, dz, system_mode, profile_mode, surface, params, mu_max, vmax, qmin, half_sat, &
               phyto_sink, zoo_ingestion, preference, pp, zp, qp, npred, detp, dp2, dz2, dq2, &
               dn2, dd2, irr2, nutlim2, grazlim2, lightlim2, growth2, grazing2)

      p = p + 0.5_real64 * dt * (dp1 + dp2)
      z = z + 0.5_real64 * dt * (dz1 + dz2)
      q = q + 0.5_real64 * dt * (dq1 + dq2)
      nutrient = nutrient + 0.5_real64 * dt * (dn1 + dn2)
      det = det + 0.5_real64 * dt * (dd1 + dd2)
      call enforce_state(m, qmin, p, z, q, nutrient, det, dz, correction, status)
      if (status /= 0) return

      if (mod(step, output_every) == 0) then
        output_index = output_index + 1
        call rhs(m, dz, system_mode, profile_mode, surface, params, mu_max, vmax, qmin, half_sat, &
                 phyto_sink, zoo_ingestion, preference, p, z, q, nutrient, det, dp1, dz1, dq1, &
                 dn1, dd1, irr1, nutlim1, grazlim1, lightlim1, growth1, grazing1)
        call save_output(output_index, m, dz, surface, params, phyto_sink, carbon_cell, p, z, q, &
                         nutrient, det, irr1, nutlim1, grazlim1, lightlim1, growth1, grazing1, &
                         initial_mass, correction, p_out, z_out, q_out, n_out, det_out, surface_out, &
                         irradiance_out, nutlim_out, grazlim_out, lightlim_out, growth_out, grazing_out, &
                         export_p_out, export_det_out, mass_error_out, correction_out)
      end if
    end do

    current_mass = inventory(nutrient, p, z, det, dz)
    if (system_mode == 0 .and. profile_mode == 0 .and. abs(current_mass - initial_mass) > 1.0e-6_real64 * &
        max(abs(initial_mass), 1.0_real64)) status = 20
  end subroutine run_model


  subroutine rhs(m, dzgrid, system_mode, profile_mode, surface, params, mu_max, vmax, qmin, &
                 half_sat, phyto_sink, zoo_ingestion, preference, p, z, q, nutrient, det, &
                 dpdt, dzdt, dqdt, dndt, dddt, irradiance, nutlim, grazlim, lightlim, growth, &
                 grazing_loss)
    integer, intent(in) :: m, system_mode, profile_mode
    real(real64), intent(in) :: dzgrid, surface, params(nparams)
    real(real64), intent(in) :: mu_max(nphy), vmax(nphy), qmin(nphy), half_sat(nphy)
    real(real64), intent(in) :: phyto_sink(nphy), zoo_ingestion(nzoo)
    real(real64), intent(in) :: preference(nphy, nzoo)
    real(real64), intent(in) :: p(nphy, m), z(nzoo, m), q(nphy, m)
    real(real64), intent(in) :: nutrient(m), det(m)
    real(real64), intent(out) :: dpdt(nphy, m), dzdt(nzoo, m), dqdt(nphy, m)
    real(real64), intent(out) :: dndt(m), dddt(m), irradiance(m)
    real(real64), intent(out) :: nutlim(nphy, m), grazlim(nphy, m)
    real(real64), intent(out) :: lightlim(nphy, m), growth(nphy, m)
    real(real64), intent(out) :: grazing_loss(nphy, m)

    real(real64) :: monod(nphy), production(nphy), uptake_flux(nphy)
    real(real64) :: grazing_by_z(nzoo), total_prey, graze, ztotal, zloss_total
    real(real64) :: cumulative_p, depth, diffusion_d, nutrient_source
    real(real64) :: pleft, pright, zleft, zright, dleft, dright, nleft, nright
    real(real64) :: wdet, supply
    integer :: i, j, k

    dpdt = 0.0_real64
    dzdt = 0.0_real64
    dqdt = 0.0_real64
    dndt = 0.0_real64
    dddt = 0.0_real64
    grazing_loss = 0.0_real64
    cumulative_p = 0.0_real64

    supply = params(5)
    diffusion_d = params(6) * 86400.0_real64
    wdet = params(8)
    if (system_mode == 0) then
      supply = 0.0_real64
      diffusion_d = 0.0_real64
      wdet = 0.0_real64
    end if

    do k = 1, m
      depth = (real(k, real64) - 0.5_real64) * dzgrid
      cumulative_p = cumulative_p + sum(max(p(:, k), 0.0_real64)) * dzgrid
      irradiance(k) = surface * exp(-params(11) * depth - params(12) * cumulative_p)
      lightlim(:, k) = 1.0_real64 - exp(-irradiance(k) * params(13) * params(14) / mu_max)
      monod = max(nutrient(k), 0.0_real64) / (max(nutrient(k), 0.0_real64) + half_sat)
      nutlim(:, k) = 1.0_real64 - monod
      growth(:, k) = mu_max * (1.0_real64 - qmin / max(q(:, k), qmin))
      production = growth(:, k) * lightlim(:, k) * p(:, k)
      uptake_flux = vmax * monod * p(:, k) / max(q(:, k), qmin)

      grazing_by_z = 0.0_real64
      do j = 1, nzoo
        total_prey = sum(preference(:, j) * p(:, k))
        do i = 1, nphy
          graze = zoo_ingestion(j) * z(j, k) * preference(i, j) * p(i, k) / &
                  (params(19) + total_prey)
          grazing_loss(i, k) = grazing_loss(i, k) + graze
          grazing_by_z(j) = grazing_by_z(j) + graze
        end do
      end do
      grazlim(:, k) = grazing_loss(:, k) / max(mu_max * max(p(:, k), tiny_value), tiny_value)
      ztotal = sum(z(:, k))
      zloss_total = params(18) * ztotal * ztotal

      dpdt(:, k) = production - grazing_loss(:, k) - params(15) * p(:, k)
      dzdt(:, k) = params(16) * grazing_by_z - params(18) * z(:, k) * ztotal
      dqdt(:, k) = vmax * monod - growth(:, k) * q(:, k)
      if (profile_mode == 0) then
        dndt(k) = -sum(production)
      else
        dndt(k) = -sum(uptake_flux)
      end if
      dndt(k) = dndt(k) + (1.0_real64 - params(16) - params(17)) * sum(grazing_loss(:, k)) &
                 + params(9) * det(k)
      dddt(k) = -params(9) * det(k) + params(15) * sum(p(:, k)) + params(10) * zloss_total &
                + params(17) * sum(grazing_loss(:, k))

      if (system_mode == 1) then
        if (profile_mode == 0) then
          nutrient_source = params(26) * supply
        else
          nutrient_source = real(k, real64) * params(26) * supply
        end if
        if (k > 1 .and. k < m) dndt(k) = dndt(k) + nutrient_source
        if (k == m) dndt(k) = dndt(k) + supply
      end if

      if (k == 1) then
        nleft = nutrient(k)
        dleft = det(k)
      else
        nleft = nutrient(k - 1)
        dleft = det(k - 1)
      end if
      if (k == m) then
        nright = nutrient(k)
        dright = det(k)
      else
        nright = nutrient(k + 1)
        dright = det(k + 1)
      end if
      dndt(k) = dndt(k) + diffusion_d * (nright - 2.0_real64 * nutrient(k) + nleft) / dzgrid**2
      dddt(k) = dddt(k) + diffusion_d * (dright - 2.0_real64 * det(k) + dleft) / dzgrid**2
      if (k == 1) then
        dddt(k) = dddt(k) - wdet * det(k) / dzgrid
      else
        dddt(k) = dddt(k) - wdet * (det(k) - det(k - 1)) / dzgrid
      end if

      do i = 1, nphy
        if (k == 1) then
          pleft = p(i, k)
        else
          pleft = p(i, k - 1)
        end if
        if (k == m) then
          pright = p(i, k)
        else
          pright = p(i, k + 1)
        end if
        dpdt(i, k) = dpdt(i, k) + diffusion_d * (pright - 2.0_real64 * p(i, k) + pleft) / dzgrid**2
        if (k == 1) then
          dpdt(i, k) = dpdt(i, k) - phyto_sink(i) * p(i, k) / dzgrid
        else
          dpdt(i, k) = dpdt(i, k) - phyto_sink(i) * (p(i, k) - p(i, k - 1)) / dzgrid
        end if
      end do
      do j = 1, nzoo
        if (k == 1) then
          zleft = z(j, k)
        else
          zleft = z(j, k - 1)
        end if
        if (k == m) then
          zright = z(j, k)
        else
          zright = z(j, k + 1)
        end if
        dzdt(j, k) = dzdt(j, k) + diffusion_d * (zright - 2.0_real64 * z(j, k) + zleft) / dzgrid**2
      end do
    end do
  end subroutine rhs


  subroutine enforce_state(m, qmin, p, z, q, nutrient, det, dzgrid, correction, status)
    integer, intent(in) :: m
    real(real64), intent(in) :: qmin(nphy), dzgrid
    real(real64), intent(inout) :: p(nphy, m), z(nzoo, m), q(nphy, m)
    real(real64), intent(inout) :: nutrient(m), det(m), correction
    integer, intent(inout) :: status
    real(real64) :: negative_mass
    integer :: k

    if (any(.not. ieee_finite_array(p)) .or. any(.not. ieee_finite_array(z)) .or. &
        any(.not. ieee_finite_array(q)) .or. any(.not. ieee_finite_array(nutrient)) .or. &
        any(.not. ieee_finite_array(det))) then
      status = 30
      return
    end if
    if (minval(p) < -1.0e-9_real64 .or. minval(z) < -1.0e-9_real64 .or. &
        minval(nutrient) < -1.0e-9_real64 .or. minval(det) < -1.0e-9_real64) then
      status = 31
      return
    end if
    negative_mass = -sum(min(p, 0.0_real64)) - sum(min(z, 0.0_real64)) &
                    - sum(min(nutrient, 0.0_real64)) - sum(min(det, 0.0_real64))
    correction = correction + negative_mass * dzgrid
    p = max(p, 0.0_real64)
    z = max(z, 0.0_real64)
    nutrient = max(nutrient, 0.0_real64)
    det = max(det, 0.0_real64)
    do k = 1, m
      q(:, k) = max(q(:, k), qmin)
    end do
  end subroutine enforce_state


  elemental logical function ieee_finite_array(value) result(ok)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(real64), intent(in) :: value
    ok = ieee_is_finite(value)
  end function ieee_finite_array


  subroutine save_output(index, m, dzgrid, surface, params, phyto_sink, carbon_cell, p, z, q, &
                         nutrient, det, irradiance, nutlim, grazlim, lightlim, growth, grazing, &
                         initial_mass, correction, p_out, z_out, q_out, n_out, det_out, surface_out, &
                         irradiance_out, nutlim_out, grazlim_out, lightlim_out, growth_out, grazing_out, &
                         export_p_out, export_det_out, mass_error_out, correction_out)
    integer, intent(in) :: index, m
    real(real64), intent(in) :: dzgrid, surface, params(nparams), phyto_sink(nphy), carbon_cell(nphy)
    real(real64), intent(in) :: p(nphy, m), z(nzoo, m), q(nphy, m), nutrient(m), det(m)
    real(real64), intent(in) :: irradiance(m), nutlim(nphy, m), grazlim(nphy, m)
    real(real64), intent(in) :: lightlim(nphy, m), growth(nphy, m), grazing(nphy, m)
    real(real64), intent(in) :: initial_mass, correction
    real(real64), intent(inout) :: p_out(:, :, :), z_out(:, :, :), q_out(:, :, :)
    real(real64), intent(inout) :: n_out(:, :, :), det_out(:, :), surface_out(:)
    real(real64), intent(inout) :: irradiance_out(:, :), nutlim_out(:, :, :)
    real(real64), intent(inout) :: grazlim_out(:, :, :), lightlim_out(:, :, :)
    real(real64), intent(inout) :: growth_out(:, :, :), grazing_out(:, :, :)
    real(real64), intent(inout) :: export_p_out(:), export_det_out(:)
    real(real64), intent(inout) :: mass_error_out(:), correction_out(:)
    integer :: k

    p_out(:, index, :) = p
    z_out(:, index, :) = z
    q_out(:, index, :) = q
    do k = 1, m
      n_out(1, index, k) = nutrient(k)
      n_out(2, index, k) = nutrient(k) / 16.0_real64
      n_out(3, index, k) = 106.0_real64 * nutrient(k) / 16.0_real64
    end do
    det_out(index, :) = det
    surface_out(index) = surface
    irradiance_out(index, :) = irradiance
    nutlim_out(:, index, :) = nutlim
    grazlim_out(:, index, :) = grazlim
    lightlim_out(:, index, :) = lightlim
    growth_out(:, index, :) = growth
    grazing_out(:, index, :) = grazing
    export_p_out(index) = sum(phyto_sink * p(:, m) / max(q(:, m), tiny_value) * carbon_cell) * 1.0e-9_real64
    export_det_out(index) = params(8) * det(m) * 6.625_real64 * 12.0_real64
    mass_error_out(index) = inventory(nutrient, p, z, det, dzgrid) - initial_mass
    correction_out(index) = correction
  end subroutine save_output


  real(real64) function inventory(nutrient, p, z, det, dzgrid) result(total)
    real(real64), intent(in) :: nutrient(:), p(:, :), z(:, :), det(:), dzgrid
    total = dzgrid * (sum(nutrient) + sum(p) + sum(z) + sum(det))
  end function inventory


  real(real64) function forcing_value(time, duration, mode, params, state) result(value)
    real(real64), intent(in) :: time, duration, params(nparams)
    integer, intent(in) :: mode
    integer(int64), intent(inout) :: state
    real(real64) :: baseline

    baseline = params(1) * params(2)
    select case (mode)
    case (0)
      value = baseline
    case (1)
      value = baseline * (0.5_real64 + 0.5_real64 * sin(2.0_real64 * pi * time / 365.0_real64)**2)
    case (2)
      value = baseline * sin(pi * time)**2
    case (3)
      value = baseline + params(3) * time / max(duration, tiny_value)
    case (4)
      value = max(0.0_real64, baseline + params(4) * normal_draw(state))
    case default
      value = -1.0_real64
    end select
  end function forcing_value


  real(real64) function normal_draw(state) result(value)
    integer(int64), intent(inout) :: state
    real(real64) :: u1, u2
    state = modulo(16807_int64 * state, 2147483647_int64)
    u1 = real(state, real64) / 2147483647.0_real64
    state = modulo(16807_int64 * state, 2147483647_int64)
    u2 = real(state, real64) / 2147483647.0_real64
    value = sqrt(-2.0_real64 * log(max(u1, tiny_value))) * cos(2.0_real64 * pi * u2)
  end function normal_draw


  subroutine logspace(xmin, xmax, values)
    real(real64), intent(in) :: xmin, xmax
    real(real64), intent(out) :: values(:)
    integer :: i, n
    n = size(values)
    do i = 1, n
      values(i) = 10.0_real64**(log10(xmin) + real(i - 1, real64) * &
                  (log10(xmax) - log10(xmin)) / real(n - 1, real64))
    end do
  end subroutine logspace

end module splas_core
