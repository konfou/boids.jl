## VISIBLE SECTION
#

function flock_modded!(B; nfun=neighbors_n, wind_accel=zeros(2), tau=0.1)
  B0 = deepcopy(B)
  @simd for bi in B
    deli = nfun(bi, B0)

    si = seperate(bi, deli)
    ai = align(bi, deli)
    ci = cohere(bi, deli)

    wr = wander(bi)
    bp = bound_position(bi)

    bi.vel += bi.ws*si + bi.wa*ai + bi.wc*ci + bi.ww*wr + bi.wb*bp
    limit_velocity!(bi)
    bi.vel += wind(wind_accel)
    bi.pos += tau*bi.vel
  end
end

function neighbors_n(bi, B)
  deli = partialsort(collect(B), 2:(bi.nc + 1), by = bj -> dis(bi, bj))
  return deli
end

function limit_velocity!(bi)
  vmag = norm(bi.vel)
  if vmag > bi.vmax
    bi.vel *= bi.vmax/vmag
  end
  if vmag < bi.vmin
    bi.vel *= bi.vmin/vmag
  end
end

function wind(accel)
  return accel
end

function wander(bi)
  if bi.ww == 0.
    return zeros(dim(bi))
  else
    return bi.vel .* rand(dim(bi))
  end
end

function bound_position(bi)
  dv = zeros(dim(bi))
  for (i, pi) in enumerate(bi.pos)
    if pi < bi.bmin[i]
      dv[i] = bi.vmax
    elseif pi > bi.bmax[i]
      dv[i] = -bi.vmax
    end
  end
  return dv
end

## HIDDEN SECTION
#

function run_modded(;
  dim=2,
  nb=100, dv=5., av=3.05, nc=4, ws=.1, wa=.05, wc=.01, ww=0., nfun=neighbors_n,
  wind_accel=zeros(dim),
  nmax=1000
)
  B, W = initialize(nb=nb, dv=dv, av=av, nc=nc, ws=ws, wa=wa, wc=wc, ww=ww, dim=dim)
  t = time()
  display(W.sc)

  try
    @fastmath for _ = 1:nmax
      t = draw!(W, B, t=t, dim=dim)
      #stat!(W, B)
      flock_modded!(B, nfun=nfun, wind_accel=wind_accel)
    end
  catch InterruptException
    println("\nSimulation interrupted.")
  end
end
