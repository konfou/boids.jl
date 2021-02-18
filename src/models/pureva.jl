## VISIBLE SECTION
#

function flock!(Ba, Bb; W, nfun=neighbors_n, tau=0.1)
  Ba0 = deepcopy(Ba)
  Bb0 = deepcopy(Bb)
  @simd for bi in Ba
    deli = neighbors(bi, Ba0)

    si = seperate(bi, deli)
    ai = align(bi, deli)
    ci = cohere(bi, deli)
    ev = evade(bi, Bb0, c=tau)
    bp = bound_position(bi)

    bi.vel += bi.ws*si + bi.wa*ai + bi.wc*ci + bi.we*ev + bi.wb*bp
    limit_velocity!(bi)
    bi.pos += tau*bi.vel
  end
  @simd for bi in Bb
    si = seperate(bi, Bb0)
    pu = pursuit(bi, Ba0, c=tau)
    bp = bound_position(bi)

    bi.vel += bi.ws*si + bi.wp*pu + bi.wb*bp
    limit_velocity!(bi)
    bi.pos += tau*bi.vel
  end
end

function pursuit(bi, O; c)
  j, cg = 0, zeros(dim(bi))
  for bj in O
    if dis(bi, bj) < bi.dp
      cg += bj.pos + c*dis(bi, bj)*bj.vel
      j += 1
    end
  end
  if j != 0
    pit = cg/j - bi.pos
    dv = (pit/norm(pit))*bi.vmax - bi.vel
    return dv
  else
    dv = cg #zeros(dim(bi))
    return dv
  end
end

function evade(bi, O; c)
  bj = partialsort(collect(O), 1, by = bj -> dis(bi, bj))
  if dis(bi, bj) < bi.de
    pit = (bj.pos + c*dis(bi, bj)*bj.vel) - bi.pos
    dv = -(pit/norm(pit))*bi.vmax - bi.vel
    return dv
  else
    dv = zeros(dim(bi))
    return dv
  end
end

## HIDDEN SECTION
#

function initialize_pureva(;
  dim=2, nbb=2,
  nba=100, dv=5., av=3.05, ds=dv/2., nc=4, de=3*dv, dp=3*dv,
  w0=1., ws=.1, wa=.05, wc=.01, ww=0., wb=.1, wp=.01, we=.1,
  vmax=5., vmin=.4*vmax,
  bmin=-(200/dim^2)*ones(dim), bmax=(200/dim^2)*ones(dim)
)
  Random.seed!(0)
  pfun = x -> bmin + x .* (bmax .- bmin)
  vfun = x -> (x ./ norm(x)) .* vmax / 2
  Bfun = nb -> [Boid(i, pfun(rand(dim)), vfun(rand(dim)*2 - ones(dim)),
    dv, av, ds, w0, ws, wa, wc, vmax, vmin, bmin, bmax, ww, wb, wp, we, nc, dp, de) for i = 1:nb]
  if dim == 2
    Wxy = B -> Node(collect(Point2f0(b.pos) for b in B))
  else
    Wxy = B -> Node(collect(Point3f0(b.pos) for b in B))
  end
  B = map(Bfun, [nba, nbb])
  W = World()
  W.xy = map(Wxy, B)
  W.sc = scatter(W.xy[1])
  scatter!(W.xy[2], color = "red")
  return (B[1], B[2], W)
end

function draw!(W, Ba, Bb; t, dim=2, fps=60)
  if dim == 2
    Wxy = B -> collect(Point2f0(b.pos) for b in B)
  else
    Wxy = B -> collect(Point3f0(b.pos) for b in B)
  end
  W.xy[1][] = Wxy(Ba)
  W.xy[2][] = Wxy(Bb)
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

function run_pureva(; nba=100, nbb=2, nfun=neighbors_n, nmax=1000, dim=2)
  Ba, Bb, W = initialize_pureva(nba=nba, nbb=nbb, dim=dim)
  t = time()
  display(W.sc)

  try
    @fastmath for _ = 1:nmax
      t = draw!(W, Ba, Bb, t=t, dim=dim)
      #stat!(W, Ba, Bb)
      flock!(Ba, Bb, W=W, nfun=nfun)
    end
  catch InterruptException
    println("\nSimulation interrupted.")
  end
end
