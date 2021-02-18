## VISIBLE SECTION
#

function flock_basic!(B; tau=0.1)
  B0 = deepcopy(B)
  @simd for bi in B
    deli = neighbors(bi, B0)

    si = seperate(bi, deli)
    ai = align(bi, deli)
    ci = cohere(bi, deli)

    bi.vel += bi.ws*si + bi.wa*ai + bi.wc*ci
    bi.pos += tau*bi.vel
  end
end

function neighbors(bi, B)
  deli = [bj for bj in B if
    bi != bj && dis(bi, bj) < bi.dv && ang(bi, bj) < bi.av]
  return deli
end

function seperate(bi, deli)
  si = zeros(dim(bi))
  for bj in deli
    if dis(bi, bj) < bi.ds
      si -= bj.pos - bi.pos
    end
  end
  return si
end

function align(bi, deli)
  ai = zeros(dim(bi))
  if isempty(deli)
    return ai
  end
  for bj in deli
    ai += bj.vel
  end
  return ai/length(deli) - bi.vel
end

function cohere(bi, deli)
  ci = zeros(dim(bi))
  if isempty(deli)
    return ci
  end
  for bj in deli
    ci += bj.pos
  end
  return ci/length(deli) - bi.pos
end

## HIDDEN SECTION
#

function initialize(;
  dim=2,
  nb=100, dv=5., av=3.05, ds=dv/2., nc=4, de=3*dv, dp=3*dv,
  w0=1., ws=.1, wa=.05, wc=.01, ww=0., wb=.1, wp=1., we=1.,
  vmax=10., vmin=.4*vmax,
  bmin=-(200/dim^2)*ones(dim), bmax=(200/dim^2)*ones(dim)
)
  Random.seed!(0)
  pfun = x -> bmin + x .* (bmax .- bmin)
  vfun = x -> (x ./ norm(x)) .* vmax / 2
  B = [Boid(i, pfun(rand(dim)), vfun(rand(dim)*2 - ones(dim)),
    dv, av, ds, w0, ws, wa, wc, vmax, vmin, bmin, bmax, ww, wb, wp, we, nc, dp, de) for i = 1:nb]
  W = World()
  if dim == 2
    W.xy = Node(collect(Point2f0(b.pos) for b in B))
  else
    W.xy = Node(collect(Point3f0(b.pos) for b in B))
  end
  W.sc = scatter(W.xy)
  return (B, W)
end

function draw!(W, B; t, dim=2, fps=60)
  if dim == 2
    W.xy[] = collect(Point2f0(b.pos) for b in B)
  else
    W.xy[] = collect(Point3f0(b.pos) for b in B)
  end
  dt = time() - t
  if dt > 1/fps
    sleep(1/fps)
    return time()
  else
    sleep(1/fps - dt)
    return time()
  end
  #sleep(1/fps)
  #return time()
end

function stat!(W, B)
  # TODO: make relevant draw function
  n = dim(first(B))
  pgv = zeros(n)
  cg = zeros(n)
  for bi in B
    pgv += bi.vel/norm(bi.vel)
    cg += bi.pos
  end
  mgv = zeros(n, n)
  for bi in B
    mgv += (bi.pos - cg)*bi.vel'
  end
  push!(W.sc, norm(pgv)/n, norm(mgv)/n)
  for bi in B
  end
end

function run_basic(; nb=100, nmax=1000, dim=2)
  B, W = initialize(nb=nb, dim=dim)
  t = time()
  display(W.sc)

  try
    @fastmath for _ = 1:nmax
      t = draw!(W, B, t=t, dim=dim)
      flock_basic!(B)
    end
  catch InterruptException
    println("\nSimulation interrupted.")
  end
end
