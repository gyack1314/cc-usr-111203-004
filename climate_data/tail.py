"""超阈值尾部分析：runs去簇、GPD拟合与Poisson簇重现水平。

约定与docs/data-contract.md一致：
- 段内缺行或null即整段拒绝，不补干日；段重叠、日期重复、负雨量由domain.summarize拒绝。
- 严格超过u的日按runs法成簇：相邻超阈值日之间至少有run_length个连续已观测
  非超阈值日才分成两簇；簇峰取最大值，同峰取最早日。
- 段边缘不足run_length个连续已观测非超阈值日界定的簇按边界簇排除；暴露不变。
- 暴露=通过最短完整段长度门槛的完整段日数/365.2425，含边界观察日。
- GPD参数(sigma>0, xi)，xi趋零用稳定极限；支持域1+xi*z/sigma>0为硬约束，
  越界记不可行（NLL=inf），不是可被容差忽略的小误差。
- 重现水平：年最大值不超概率P(annual max<=u+z)=exp(-lambda*S(z))，目标1-1/T；
  所需生存概率>1时水平落在阈值以下，超出尾部适用域，不外推；有限支持上端点不可越过。
- 求解只给有限矩形内的多起点局部候选、网格扫描与KKT/边界诊断，不宣称全局最优。
"""
import math
from datetime import timedelta

from .domain import summarize
from .io import date, integer, number

DAYS_PER_YEAR = 365.2425

_XI_TINY = 1e-6    # |xi|小于此值时NLL/梯度改用xi->0稳定级数
_XI_ZERO = 1e-12   # |xi|小于此值时生存/逆函数取指数极限
_GRID = 41         # 矩形粗网格每维点数（全局扫描证据）
_START_FRAC_S = (0.2, 0.5, 0.8)
_START_FRAC_X = (0.0, 0.25, 0.5, 0.75, 1.0)
_NM_MAXITER = 3000
_NM_XTOL = 1e-11
_NM_FTOL = 1e-12
_KKT_GTOL = 1e-5
_BOUND_RTOL = 1e-8

CONVENTIONS = [
    '边界簇排除仅为估计约定：段边缘不足run_length个连续已观测非超阈值日的簇不入选，'
    '不声称因此消除簇间全部相关性',
    '暴露=通过最短完整段长度门槛的完整段日数/365.2425，包含段边界观察日；'
    '缺日或含null的段整段拒绝，不补干日',
    'GPD以(sigma,xi)参数化，xi=0取指数极限；支持域1+xi*z/sigma>0为硬约束，'
    '越界记为不可行而非容差忽略',
    '重现水平基于Poisson簇到达+GPD尾部：P(年最大值<=u+z)=exp(-lambda*S(z))；'
    '所需生存概率>1时水平在阈值以下，超出尾部适用域，不外推到分布主体',
    '合成资料，不涉及真实预警发布',
]

LIMITATIONS = [
    '全局最优未证明：报告矩形内全部起点结果、粗网格扫描与KKT/边界诊断，'
    '候选仅为最佳局部结果',
    '敏感性各组(u,run_length)结果全部保留，未择优剔除',
]

_CONFIG_KEYS = {'u', 'run_length', 'min_segment_days', 'sigma_bounds', 'xi_bounds',
                'return_periods', 'sensitivity'}


def check_config(c):
    """校验分析前给定的配置，返回规范化副本。"""
    if not isinstance(c, dict):
        raise ValueError('配置应为对象')
    unknown = set(c) - _CONFIG_KEYS
    if unknown:
        raise ValueError('配置含未知键: ' + ','.join(sorted(unknown)))
    u = number(c['u'])
    if u < 0:
        raise ValueError('阈值u不能负')
    r = integer(c['run_length'], 0)
    min_days = integer(c['min_segment_days'], 1)
    sb = c['sigma_bounds']
    if not (isinstance(sb, list) and len(sb) == 2):
        raise ValueError('sigma_bounds需为[下界,上界]')
    s_lo, s_hi = number(sb[0]), number(sb[1])
    if not 0 < s_lo < s_hi:
        raise ValueError('sigma边界需为正有限且递增')
    xb = c['xi_bounds']
    if not (isinstance(xb, list) and len(xb) == 2):
        raise ValueError('xi_bounds需为[下界,上界]')
    x_lo, x_hi = number(xb[0]), number(xb[1])
    if not x_lo < x_hi:
        raise ValueError('xi边界需为有限且递增')
    periods = c['return_periods']
    if not isinstance(periods, list) or not periods:
        raise ValueError('return_periods应为非空列表')
    ts = []
    for t in periods:
        t = number(t)
        if t <= 1:
            raise ValueError('重现期T需>1')
        ts.append(t)
    sens = c.get('sensitivity', [])
    if not isinstance(sens, list):
        raise ValueError('sensitivity应为列表')
    pairs = []
    for e in sens:
        if not isinstance(e, dict) or set(e) != {'u', 'run_length'}:
            raise ValueError('敏感性条目需且仅需含u与run_length')
        eu = number(e['u'])
        if eu < 0:
            raise ValueError('阈值u不能负')
        pairs.append((eu, integer(e['run_length'], 0)))
    return {'u': u, 'run_length': r, 'min_segment_days': min_days,
            'sigma_bounds': (s_lo, s_hi), 'xi_bounds': (x_lo, x_hi),
            'return_periods': ts, 'sensitivity': pairs}


def _audit_segments(v, min_days):
    """逐段核对完整性：缺行或null即拒绝该段；短于门槛的段拒绝。"""
    day_mm = {}
    for row in v['rain']:
        day_mm[date(row['day'])] = row['mm']
    accepted, rejected = [], []
    for s in sorted(v['segments'], key=lambda x: x['start']):
        lo, hi = date(s['start']), date(s['end'])
        n = (hi - lo).days + 1
        missing, null_days = [], []
        cur = lo
        while cur <= hi:
            if cur not in day_mm:
                missing.append(cur)
            elif day_mm[cur] is None:
                null_days.append(cur)
            cur += timedelta(days=1)
        base = {'id': s['id'], 'start': s['start'], 'end': s['end'], 'calendar_days': n}
        if missing or null_days:
            rejected.append({**base, 'reason': 'incomplete',
                             'missing_count': len(missing),
                             'missing_first': [d.isoformat() for d in missing[:5]],
                             'null_count': len(null_days),
                             'null_first': [d.isoformat() for d in null_days[:5]]})
            continue
        if n < min_days:
            rejected.append({**base, 'reason': 'below-min-segment-days',
                             'min_segment_days': min_days})
            continue
        series = [day_mm[lo + timedelta(days=i)] for i in range(n)]
        accepted.append({**base, 'series': series})
    return accepted, rejected


def _decluster(seg, u, r):
    """runs法成簇并按段边缘run_length规则区分入选簇与边界簇。"""
    xs = seg['series']
    total = len(xs)
    groups = []
    for i, x in enumerate(xs):
        if x > u:
            if groups and i - groups[-1][-1] - 1 < r:
                groups[-1].append(i)
            else:
                groups.append([i])
    accepted, excluded = [], []
    day0 = date(seg['start'])
    for k, g in enumerate(groups):
        f, l = g[0], g[-1]
        before = f if k == 0 else f - groups[k - 1][-1] - 1
        after = (total - 1 - l) if k == len(groups) - 1 else groups[k + 1][0] - l - 1
        pk = g[0]
        for i in g[1:]:
            if xs[i] > xs[pk]:
                pk = i
        rec = {'segment': seg['id'],
               'first_day': (day0 + timedelta(days=f)).isoformat(),
               'last_day': (day0 + timedelta(days=l)).isoformat(),
               'exceedance_days': len(g),
               'peak_day': (day0 + timedelta(days=pk)).isoformat(),
               'peak_mm': xs[pk],
               'excess_mm': xs[pk] - u,
               'non_exceedance_run_before': before,
               'non_exceedance_run_after': after}
        if before >= r and after >= r:
            accepted.append(rec)
        else:
            sides = []
            if before < r:
                sides.append('段前缘')
            if after < r:
                sides.append('段后缘')
            rec['excluded_reason'] = ('边界簇：' + '与'.join(sides)
                                      + '连续已观测非超阈值日不足run_length=' + str(r))
            excluded.append(rec)
    return accepted, excluded


def gpd_nll(sigma, xi, z):
    """GPD负对数似然；支持域外返回inf（硬约束，不容差）。xi趋零用稳定级数。"""
    if not (math.isfinite(sigma) and math.isfinite(xi)) or sigma <= 0:
        return math.inf
    n = len(z)
    if abs(xi) < _XI_TINY:
        total = 0.0
        for zi in z:
            w = zi / sigma
            if 1 + xi * w <= 0:
                return math.inf
            total += w + xi * (w - 0.5 * w * w) + xi * xi * (w * w * w / 3.0 - 0.5 * w * w)
        return n * math.log(sigma) + total
    s = 0.0
    for zi in z:
        t = 1 + xi * zi / sigma
        if t <= 0:
            return math.inf
        s += math.log1p(xi * zi / sigma)
    return n * math.log(sigma) + (1 + 1 / xi) * s


def gpd_grad(sigma, xi, z):
    """NLL对(sigma,xi)的解析梯度；不可行点返回(nan,nan)。"""
    n = len(z)
    if abs(xi) < _XI_TINY:
        gs, gx = n / sigma, 0.0
        for zi in z:
            w = zi / sigma
            if 1 + xi * w <= 0:
                return (math.nan, math.nan)
            gs += (-w + xi * (w * w - w)) / sigma
            gx += (w - 0.5 * w * w) + 2 * xi * (w * w * w / 3.0 - 0.5 * w * w)
        return (gs, gx)
    s1, swt = 0.0, 0.0
    for zi in z:
        w = zi / sigma
        t = 1 + xi * w
        if t <= 0:
            return (math.nan, math.nan)
        s1 += math.log1p(xi * w)
        swt += w / t
    return (n / sigma - (xi + 1) / sigma * swt,
            -s1 / (xi * xi) + (1 + 1 / xi) * swt)


def gpd_survival(z, sigma, xi):
    """超量z的GPD生存概率；xi趋零取指数极限；越过有限端点为0。"""
    if z < 0:
        raise ValueError('超量不能负')
    if abs(xi) < _XI_ZERO:
        return math.exp(-z / sigma)
    t = 1 + xi * z / sigma
    if t <= 0:
        return 0.0
    return math.exp(-math.log1p(xi * z / sigma) / xi)


def gpd_inv_survival(s, sigma, xi):
    """由生存概率s求超量z；xi趋零取指数极限。"""
    if not 0 < s <= 1:
        raise ValueError('生存概率需在(0,1]')
    if abs(xi) < _XI_ZERO:
        return -sigma * math.log(s)
    return sigma / xi * (math.exp(-xi * math.log(s)) - 1)


def _nelder_mead(f, x0, box):
    """矩形约束Nelder-Mead；越出矩形裁剪回边界，支持域外f=inf。"""
    def clip(p):
        return (min(max(p[0], box[0][0]), box[0][1]),
                min(max(p[1], box[1][0]), box[1][1]))
    sx = 0.05 * (box[0][1] - box[0][0])
    sy = 0.05 * (box[1][1] - box[1][0])
    simplex = [clip(x0), clip((x0[0] + sx, x0[1])), clip((x0[0], x0[1] + sy))]
    if simplex[1] == simplex[0]:
        simplex[1] = clip((x0[0] - sx, x0[1]))
    if simplex[2] == simplex[0]:
        simplex[2] = clip((x0[0], x0[1] - sy))
    vals = [f(p) for p in simplex]
    nfev = 3
    converged = False
    for _ in range(_NM_MAXITER):
        order = sorted(range(3), key=lambda i: vals[i])
        simplex = [simplex[i] for i in order]
        vals = [vals[i] for i in order]
        if not any(math.isfinite(v) for v in vals):
            break
        fspread = abs(vals[2] - vals[0]) if all(math.isfinite(v) for v in vals) else math.inf
        diam = max(abs(simplex[i][0] - simplex[0][0]) / (1 + abs(simplex[0][0]))
                   + abs(simplex[i][1] - simplex[0][1]) / (1 + abs(simplex[0][1]))
                   for i in (1, 2))
        if fspread <= _NM_FTOL * (1 + abs(vals[0])) and diam <= _NM_XTOL:
            converged = True
            break
        cx = (simplex[0][0] + simplex[1][0]) / 2
        cy = (simplex[0][1] + simplex[1][1]) / 2
        xr = clip((2 * cx - simplex[2][0], 2 * cy - simplex[2][1]))
        fr = f(xr)
        nfev += 1
        if fr < vals[0]:
            xe = clip((2 * xr[0] - cx, 2 * xr[1] - cy))
            fe = f(xe)
            nfev += 1
            if fe < fr:
                simplex[2], vals[2] = xe, fe
            else:
                simplex[2], vals[2] = xr, fr
        elif fr < vals[1]:
            simplex[2], vals[2] = xr, fr
        else:
            if fr < vals[2]:
                xc = clip(((cx + xr[0]) / 2, (cy + xr[1]) / 2))
                fc = f(xc)
                nfev += 1
                if fc <= fr:
                    simplex[2], vals[2] = xc, fc
                    continue
            else:
                xc = clip(((cx + simplex[2][0]) / 2, (cy + simplex[2][1]) / 2))
                fc = f(xc)
                nfev += 1
                if fc < vals[2]:
                    simplex[2], vals[2] = xc, fc
                    continue
            simplex[1] = clip(((simplex[0][0] + simplex[1][0]) / 2,
                               (simplex[0][1] + simplex[1][1]) / 2))
            simplex[2] = clip(((simplex[0][0] + simplex[2][0]) / 2,
                               (simplex[0][1] + simplex[2][1]) / 2))
            vals[1] = f(simplex[1])
            vals[2] = f(simplex[2])
            nfev += 2
    order = sorted(range(3), key=lambda i: vals[i])
    return simplex[order[0]], vals[order[0]], nfev, converged


def _pwm_start(z):
    """概率加权矩起点（仅作多起点之一，不作估计量）。"""
    n = len(z)
    zs = sorted(z)
    m0 = sum(zs) / n
    m1 = sum((1 - (i + 1 - 0.35) / n) * v for i, v in enumerate(zs)) / n
    denom = m0 - 2 * m1
    if denom == 0:
        return None
    sigma = 2 * m0 * m1 / denom
    xi = 2 - m0 / denom
    if not (math.isfinite(sigma) and math.isfinite(xi)) or sigma <= 0:
        return None
    return (sigma, xi)


def _kkt(sigma, xi, z, box):
    """候选点的KKT/边界诊断：解析梯度、参数边界状态与支持余量。"""
    gs, gx = gpd_grad(sigma, xi, z)

    def state(v, lo, hi):
        tol = _BOUND_RTOL * (hi - lo)
        if v - lo <= tol:
            return 'lower-bound'
        if hi - v <= tol:
            return 'upper-bound'
        return 'interior'

    def ok(g, st):
        if not math.isfinite(g):
            return False
        if st == 'interior':
            return abs(g) <= _KKT_GTOL
        if st == 'lower-bound':
            return g >= -_KKT_GTOL
        return g <= _KKT_GTOL

    ss = state(sigma, box[0][0], box[0][1])
    sx = state(xi, box[1][0], box[1][1])
    margin = min(1 + xi * zi / sigma for zi in z)
    return {'gradient': {'d_sigma': gs if math.isfinite(gs) else None,
                         'd_xi': gx if math.isfinite(gx) else None},
            'sigma_bound_state': ss,
            'xi_bound_state': sx,
            'kkt_sign_ok': {'sigma': ok(gs, ss), 'xi': ok(gx, sx)},
            'support_min_margin': margin,
            'gradient_tolerance': _KKT_GTOL,
            'bound_relative_tolerance': _BOUND_RTOL,
            'satisfied': ok(gs, ss) and ok(gx, sx) and margin > 0}


def fit_gpd(z, s_lo, s_hi, x_lo, x_hi):
    """有限矩形内多起点拟合GPD，报告全部起点结果、网格扫描与KKT诊断。"""
    zs = [float(v) for v in z]
    zmax = max(zs)
    box = ((s_lo, s_hi), (x_lo, x_hi))

    def objective(p):
        return gpd_nll(p[0], p[1], zs)

    def feasible(p):
        s, x = p
        s = min(max(s, s_lo), s_hi)
        x = min(max(x, x_lo), x_hi)
        if x < 0:
            need = -x * zmax
            if s <= need:
                s = need * (1 + 1e-9) + 1e-12
                if s > s_hi:
                    return None
        return (s, x)

    starts = [(s_lo + fs * (s_hi - s_lo), x_lo + fx * (x_hi - x_lo))
              for fs in _START_FRAC_S for fx in _START_FRAC_X]
    pwm = _pwm_start(zs)
    if pwm is not None:
        starts.append(pwm)

    runs = []
    for s0 in starts:
        p0 = feasible(s0)
        if p0 is None:
            runs.append({'start': [s0[0], s0[1]], 'status': 'skipped',
                         'reason': '起点在矩形内无可行化（支持域与矩形不相交）'})
            continue
        p, fv, nfev, conv = _nelder_mead(objective, p0, box)
        finite = math.isfinite(fv)
        runs.append({'start': [p0[0], p0[1]], 'end': [p[0], p[1]],
                     'nll': fv if finite else None,
                     'converged': conv, 'evaluations': nfev,
                     'status': 'ok' if conv and finite
                     else ('finite-not-converged' if finite else 'failed')})

    feasible_points = 0
    grid_best = None
    for i in range(_GRID):
        for j in range(_GRID):
            s = s_lo + (i + 0.5) / _GRID * (s_hi - s_lo)
            x = x_lo + (j + 0.5) / _GRID * (x_hi - x_lo)
            v = gpd_nll(s, x, zs)
            if math.isfinite(v):
                feasible_points += 1
                if grid_best is None or v < grid_best['nll']:
                    grid_best = {'sigma': s, 'xi': x, 'nll': v}

    finite_runs = [rn for rn in runs if rn.get('nll') is not None]
    candidate = None
    if finite_runs:
        best = min(finite_runs, key=lambda rn: rn['nll'])
        candidate = {'sigma': best['end'][0], 'xi': best['end'][1], 'nll': best['nll']}

    groups = {}
    for rn in finite_runs:
        key = (round(rn['end'][0], 6), round(rn['end'][1], 6))
        g = groups.setdefault(key, {'sigma': rn['end'][0], 'xi': rn['end'][1],
                                    'nll': rn['nll'], 'hit_count': 0})
        g['hit_count'] += 1
        if rn['nll'] < g['nll']:
            g.update(sigma=rn['end'][0], xi=rn['end'][1], nll=rn['nll'])
    optima = sorted(groups.values(), key=lambda g: g['nll'])

    return {'n_excesses': len(zs),
            'excesses_mm': zs,
            'rectangle': {'sigma': [s_lo, s_hi], 'xi': [x_lo, x_hi]},
            'starts': runs,
            'grid_scan': {'points_per_dimension': _GRID,
                          'feasible_points': feasible_points,
                          'best': grid_best},
            'distinct_local_optima': optima,
            'candidate': candidate,
            'kkt': _kkt(candidate['sigma'], candidate['xi'], zs, box) if candidate else None,
            'solver_status': 'local-candidates-only',
            'global_optimality': '未证明全局最优：候选为有限矩形内多起点局部搜索'
                                 '与粗网格扫描的最佳结果'}


def return_level(lam, sigma, xi, u, t):
    """T年重现水平；目标在阈值下或越过有限端点时不给水平值。"""
    p = 1 - 1 / t
    s = -math.log(p) / lam
    rec = {'T': t, 'annual_non_exceedance_target': p, 'required_survival': s}
    if s > 1:
        rec['status'] = 'below-threshold'
        rec['note'] = '目标水平落在阈值u以下，超出尾部模型适用域，不外推'
        return rec
    z = gpd_inv_survival(s, sigma, xi)
    rec['excess_mm'] = z
    rec['level_mm'] = u + z
    if xi < 0:
        endpoint = -sigma / xi
        rec['endpoint_excess_mm'] = endpoint
        rec['endpoint_margin_mm'] = endpoint - z
        if endpoint - z <= 0:
            rec['status'] = 'endpoint-violation'
            rec['note'] = '有限支持上端点不可越过'
            return rec
    s2 = gpd_survival(z, sigma, xi)
    z_back = gpd_inv_survival(s2, sigma, xi) if s2 > 0 else None
    g = math.exp(-lam * s2)
    rec['roundtrip'] = {'survival_at_level': s2,
                        'survival_residual': s2 - s,
                        'level_roundtrip_residual_mm':
                            (z_back - z) if z_back is not None else None,
                        'annual_non_exceedance': g,
                        'annual_target_residual': g - p}
    rec['status'] = 'ok'
    return rec


def _run_once(v, cfg, u, r):
    accepted_segs, rejected_segs = _audit_segments(v, cfg['min_segment_days'])
    exposure_days = sum(s['calendar_days'] for s in accepted_segs)
    exposure_years = exposure_days / DAYS_PER_YEAR
    acc, exc = [], []
    for s in accepted_segs:
        a, e = _decluster(s, u, r)
        acc += a
        exc += e
    out = {'u': u, 'run_length': r,
           'segments': {'accepted': [{k: s[k] for k in
                                      ('id', 'start', 'end', 'calendar_days')}
                                     for s in accepted_segs],
                        'rejected': rejected_segs},
           'exposure': {'days': exposure_days, 'years': exposure_years,
                        'days_per_year': DAYS_PER_YEAR},
           'clusters': {'accepted': acc, 'excluded_boundary': exc},
           'cluster_count_accepted': len(acc)}
    if not acc:
        out['gpd'] = None
        out['return_levels'] = None
        out['note'] = '无入选簇，不能拟合GPD，不给重现水平'
        return out
    lam = len(acc) / exposure_years
    out['lambda_per_year'] = lam
    fit = fit_gpd([c['excess_mm'] for c in acc],
                  *cfg['sigma_bounds'], *cfg['xi_bounds'])
    out['gpd'] = fit
    if fit['candidate'] is None:
        out['return_levels'] = None
        out['note'] = '矩形内未得可行候选，不给重现水平'
        return out
    out['return_levels'] = [return_level(lam, fit['candidate']['sigma'],
                                         fit['candidate']['xi'], u, t)
                            for t in cfg['return_periods']]
    return out


def analyze(v, digest, config, config_digest=None):
    """完整分析：段审核、去簇、GPD拟合、重现水平与敏感性重算。"""
    cfg = check_config(config)
    raw = summarize(v)
    pairs = [(cfg['u'], cfg['run_length'])]
    for p in cfg['sensitivity']:
        if p not in pairs:
            pairs.append(p)
    runs = [_run_once(v, cfg, u, r) for u, r in pairs]
    report = {'analysis': 'peaks-over-threshold-tail',
              'source_sha256': digest,
              'units': {'rainfall': 'mm', 'exposure': 'year(365.2425d)'},
              'raw_summary': raw,
              'parameters': {'u': cfg['u'], 'run_length': cfg['run_length'],
                             'min_segment_days': cfg['min_segment_days'],
                             'sigma_bounds': list(cfg['sigma_bounds']),
                             'xi_bounds': list(cfg['xi_bounds']),
                             'return_periods': cfg['return_periods'],
                             'sensitivity': [{'u': u, 'run_length': r}
                                             for u, r in cfg['sensitivity']]},
              'conventions': CONVENTIONS,
              'limitations': LIMITATIONS,
              'primary': runs[0],
              'sensitivity': runs[1:]}
    if config_digest is not None:
        report['config_sha256'] = config_digest
    return report
