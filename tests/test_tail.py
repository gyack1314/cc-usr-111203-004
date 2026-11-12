"""超阈值尾部分析的行为与数值核对。"""
import json
import math
import unittest
from datetime import date as dt
from datetime import timedelta

from climate_data.io import load
from climate_data.tail import (analyze, check_config, fit_gpd, gpd_grad,
                               gpd_inv_survival, gpd_nll, gpd_survival,
                               return_level)

Z = [15, 16, 17, 18, 19, 20]


def cfg(**kw):
    base = {'u': 10, 'run_length': 3, 'min_segment_days': 30,
            'sigma_bounds': [0.5, 200], 'xi_bounds': [-0.6, 0.8],
            'return_periods': [10, 20, 50, 100]}
    base.update(kw)
    return base


def make_v(start, values, seg_id='s'):
    d0 = dt.fromisoformat(start)
    rain = [{'day': (d0 + timedelta(days=i)).isoformat(), 'mm': x}
            for i, x in enumerate(values)]
    end = (d0 + timedelta(days=len(values) - 1)).isoformat()
    return {'segments': [{'id': seg_id, 'start': start, 'end': end}], 'rain': rain}


def analyze_file(name, config):
    v, digest = load('data/' + name)
    return analyze(v, digest, config)


class SegmentAuditTests(unittest.TestCase):
    def test_missing_day_rejects_segment(self):
        v = make_v('2020-01-01', [0] * 10)
        del v['rain'][5]
        out = analyze(v, 'x', cfg(min_segment_days=1))
        seg = out['primary']['segments']
        self.assertEqual(seg['accepted'], [])
        self.assertEqual(seg['rejected'][0]['reason'], 'incomplete')
        self.assertEqual(seg['rejected'][0]['missing_count'], 1)
        self.assertEqual(out['primary']['exposure']['days'], 0)

    def test_null_is_not_dry_day(self):
        v = make_v('2020-01-01', [0] * 9 + [None])
        out = analyze(v, 'x', cfg(min_segment_days=1))
        seg = out['primary']['segments']
        self.assertEqual(seg['rejected'][0]['reason'], 'incomplete')
        self.assertEqual(seg['rejected'][0]['null_count'], 1)

    def test_min_segment_days_gate(self):
        v = make_v('2020-01-01', [0] * 10)
        out = analyze(v, 'x', cfg(min_segment_days=30))
        seg = out['primary']['segments']
        self.assertEqual(seg['rejected'][0]['reason'], 'below-min-segment-days')
        self.assertEqual(out['primary']['exposure']['days'], 0)
        self.assertIsNone(out['primary']['gpd'])

    def test_overlap_rejected(self):
        v = {'segments': [{'id': 'a', 'start': '2020-01-01', 'end': '2020-01-10'},
                          {'id': 'b', 'start': '2020-01-05', 'end': '2020-01-15'}],
             'rain': []}
        with self.assertRaises(ValueError):
            analyze(v, 'x', cfg())

    def test_exposure_includes_boundary_days(self):
        out = analyze_file('boundary-clusters.json', cfg())
        self.assertEqual(out['primary']['exposure']['days'], 90)
        self.assertAlmostEqual(out['primary']['exposure']['years'],
                               90 / 365.2425)


class DeclusterTests(unittest.TestCase):
    def test_normal_clusters(self):
        out = analyze_file('normal.json', cfg())
        p = out['primary']
        self.assertEqual(p['cluster_count_accepted'], 6)
        self.assertEqual(p['clusters']['excluded_boundary'], [])
        peaks = [(c['peak_day'], c['peak_mm'], c['excess_mm'])
                 for c in p['clusters']['accepted']]
        self.assertEqual(peaks, [('2019-01-08', 25, 15), ('2019-01-23', 26, 16),
                                 ('2019-02-07', 27, 17), ('2019-02-22', 28, 18),
                                 ('2019-03-09', 29, 19), ('2019-03-24', 30, 20)])
        first, last = p['clusters']['accepted'][0], p['clusters']['accepted'][-1]
        self.assertEqual((first['non_exceedance_run_before'],
                          first['non_exceedance_run_after']), (6, 12))
        self.assertEqual((last['non_exceedance_run_before'],
                          last['non_exceedance_run_after']), (12, 6))
        self.assertAlmostEqual(p['lambda_per_year'], 6 / (90 / 365.2425))

    def test_three_day_rain_is_one_event(self):
        out = analyze_file('normal.json', cfg())
        for c in out['primary']['clusters']['accepted']:
            self.assertEqual(c['exceedance_days'], 3)

    def test_boundary_clusters_excluded(self):
        out = analyze_file('boundary-clusters.json', cfg())
        p = out['primary']
        self.assertEqual(p['cluster_count_accepted'], 6)
        exc = p['clusters']['excluded_boundary']
        self.assertEqual(len(exc), 2)
        self.assertEqual(exc[0]['peak_day'], '2019-01-01')
        self.assertEqual(exc[0]['non_exceedance_run_before'], 0)
        self.assertIn('段前缘', exc[0]['excluded_reason'])
        self.assertEqual(exc[1]['peak_day'], '2019-03-31')
        self.assertEqual(exc[1]['non_exceedance_run_after'], 0)
        self.assertIn('段后缘', exc[1]['excluded_reason'])

    def test_long_run_length_turns_edge_clusters_into_boundary(self):
        out = analyze_file('normal.json', cfg(run_length=7))
        p = out['primary']
        self.assertEqual(p['cluster_count_accepted'], 4)
        self.assertEqual(len(p['clusters']['excluded_boundary']), 2)

    def test_peak_tie_takes_earliest_day(self):
        v = make_v('2020-01-01', [0, 0, 0, 0, 0, 20, 1, 20] + [0] * 12)
        out = analyze(v, 'x', cfg(min_segment_days=1))
        clusters = out['primary']['clusters']['accepted']
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]['peak_day'], '2020-01-06')

    def test_no_exceedance_no_fit(self):
        out = analyze_file('no-exceedance.json', cfg())
        p = out['primary']
        self.assertEqual(p['cluster_count_accepted'], 0)
        self.assertIsNone(p['gpd'])
        self.assertIsNone(p['return_levels'])

    def test_incomplete_segment_file_no_fit(self):
        out = analyze_file('incomplete-segment.json', cfg())
        p = out['primary']
        self.assertEqual(p['segments']['rejected'][0]['reason'], 'incomplete')
        self.assertIsNone(p['gpd'])

    def test_invalid_domain_raises(self):
        with self.assertRaises(ValueError):
            analyze_file('invalid-domain.json', cfg())


class ConfigTests(unittest.TestCase):
    def test_rejections(self):
        for bad in (cfg(u=-1), cfg(run_length=1.5), cfg(min_segment_days=0),
                    cfg(sigma_bounds=[0, 10]), cfg(sigma_bounds=[5, 5]),
                    cfg(xi_bounds=[1, -1]), cfg(return_periods=[1]),
                    cfg(return_periods=[]), cfg(unknown_key=1),
                    cfg(sensitivity=[{'u': 8}]),
                    cfg(sensitivity=[{'u': 8, 'run_length': 3, 'x': 1}])):
            with self.subTest(bad=bad):
                with self.assertRaises((ValueError, KeyError)):
                    check_config(bad)

    def test_sensitivity_deduped_against_primary(self):
        out = analyze_file('normal.json', cfg(sensitivity=[
            {'u': 10, 'run_length': 3}, {'u': 8, 'run_length': 3}]))
        self.assertEqual(len(out['sensitivity']), 1)
        self.assertEqual(out['sensitivity'][0]['u'], 8)


class GpdMathTests(unittest.TestCase):
    def test_nll_exponential_limit_known_value(self):
        self.assertAlmostEqual(gpd_nll(10, 0.0, Z),
                               6 * math.log(10) + sum(Z) / 10, places=12)

    def test_nll_xi_tends_to_zero_continuously(self):
        base = gpd_nll(10, 0.0, Z)
        for xi in (1e-8, -1e-8, 1e-12):
            self.assertAlmostEqual(gpd_nll(10, xi, Z), base, places=5)

    def test_nll_general_branch_known_value(self):
        sigma, xi = 20.0, 0.5
        expect = 6 * math.log(sigma) + (1 + 1 / xi) * sum(
            math.log1p(xi * z / sigma) for z in Z)
        self.assertAlmostEqual(gpd_nll(sigma, xi, Z), expect, places=12)

    def test_support_violation_is_infinite_not_tolerated(self):
        self.assertEqual(gpd_nll(5, -0.5, Z), math.inf)
        self.assertEqual(gpd_nll(10.0001, -0.5, Z), gpd_nll(10.0001, -0.5, Z))
        self.assertLess(gpd_nll(10.0001, -0.5, Z), math.inf)
        self.assertEqual(gpd_nll(0, 0.5, Z), math.inf)

    def test_gradient_matches_finite_difference(self):
        for sigma, xi in ((17.3, 0.35), (25.0, -0.2), (12.0, 5e-7), (30.0, 0.0)):
            with self.subTest(sigma=sigma, xi=xi):
                gs, gx = gpd_grad(sigma, xi, Z)
                hs = 1e-6 * sigma
                fds = (gpd_nll(sigma + hs, xi, Z) - gpd_nll(sigma - hs, xi, Z)) / (2 * hs)
                hx = 1e-6
                fdx = (gpd_nll(sigma, xi + hx, Z) - gpd_nll(sigma, xi - hx, Z)) / (2 * hx)
                self.assertAlmostEqual(gs, fds, places=4)
                self.assertAlmostEqual(gx, fdx, places=4)

    def test_survival_inverse_roundtrip(self):
        for sigma, xi in ((10, 0.0), (17.5, -0.3), (12, 0.4)):
            with self.subTest(sigma=sigma, xi=xi):
                s = gpd_survival(12.3, sigma, xi)
                self.assertAlmostEqual(gpd_inv_survival(s, sigma, xi), 12.3, places=9)
        self.assertAlmostEqual(gpd_survival(5, 10, 0.0), math.exp(-0.5), places=15)
        self.assertEqual(gpd_survival(40, 10, -0.5), 0.0)
        with self.assertRaises(ValueError):
            gpd_inv_survival(1.5, 10, 0.0)


class FitTests(unittest.TestCase):
    def setUp(self):
        self.fit = fit_gpd(Z, 0.5, 200, -0.6, 0.8)

    def test_candidate_is_best_of_all_starts(self):
        cand = self.fit['candidate']
        self.assertIsNotNone(cand)
        for run in self.fit['starts']:
            if run.get('nll') is not None:
                self.assertLessEqual(cand['nll'], run['nll'] + 1e-9)
        self.assertLessEqual(cand['nll'], self.fit['grid_scan']['best']['nll'] + 1e-6)

    def test_multistart_reports_distinct_optima(self):
        self.assertEqual(len(self.fit['starts']), 16)
        self.assertTrue(all(r['status'] == 'ok' for r in self.fit['starts']))
        optima = self.fit['distinct_local_optima']
        self.assertGreaterEqual(len(optima), 2)
        self.assertAlmostEqual(optima[0]['nll'], self.fit['candidate']['nll'], places=9)
        self.assertEqual(sum(o['hit_count'] for o in optima), 16)

    def test_kkt_boundary_diagnostics(self):
        kkt = self.fit['kkt']
        self.assertEqual(kkt['xi_bound_state'], 'lower-bound')
        self.assertEqual(kkt['sigma_bound_state'], 'interior')
        self.assertTrue(kkt['satisfied'])
        self.assertGreater(kkt['support_min_margin'], 0)
        self.assertEqual(self.fit['solver_status'], 'local-candidates-only')

    def test_candidate_values_deterministic(self):
        self.assertAlmostEqual(self.fit['candidate']['sigma'], 17.7451707597822, places=9)
        self.assertAlmostEqual(self.fit['candidate']['xi'], -0.6, places=12)
        self.assertAlmostEqual(self.fit['candidate']['nll'], 20.880510905986934, places=9)


class ReturnLevelTests(unittest.TestCase):
    def test_exponential_known_value(self):
        lam, sigma, u, t = 2.0, 10.0, 10.0, 20.0
        s = -math.log(1 - 1 / t) / lam
        rec = return_level(lam, sigma, 0.0, u, t)
        self.assertEqual(rec['status'], 'ok')
        self.assertAlmostEqual(rec['level_mm'], u - sigma * math.log(s), places=9)
        self.assertAlmostEqual(rec['roundtrip']['annual_target_residual'], 0.0, places=12)

    def test_below_threshold_not_extrapolated(self):
        rec = return_level(0.5, 10.0, 0.0, 10.0, 2.0)
        self.assertEqual(rec['status'], 'below-threshold')
        self.assertNotIn('level_mm', rec)
        self.assertGreater(rec['required_survival'], 1)

    def test_endpoint_margin_reported(self):
        rec = return_level(24.3495, 17.7451707597822, -0.6, 10.0, 100.0)
        self.assertEqual(rec['status'], 'ok')
        self.assertGreater(rec['endpoint_margin_mm'], 0)
        self.assertAlmostEqual(rec['endpoint_excess_mm'],
                               17.7451707597822 / 0.6, places=9)

    def test_levels_monotone_and_roundtrip_clean(self):
        out = analyze_file('normal.json', cfg())
        levels = out['primary']['return_levels']
        seq = [r['level_mm'] for r in levels]
        self.assertEqual(seq, sorted(seq))
        for r in levels:
            self.assertEqual(r['status'], 'ok')
            self.assertGreater(r['level_mm'], 10)
            self.assertLess(abs(r['roundtrip']['survival_residual']), 1e-9)
            self.assertLess(abs(r['roundtrip']['level_roundtrip_residual_mm']), 1e-6)
            self.assertLess(abs(r['roundtrip']['annual_target_residual']), 1e-9)


class ReportTests(unittest.TestCase):
    def test_sensitivity_all_kept(self):
        out = analyze_file('normal.json', load('configs/tail-normal.json')[0])
        self.assertEqual(len(out['sensitivity']), 6)
        seen = [(s['u'], s['run_length']) for s in out['sensitivity']]
        self.assertEqual(seen, [(8, 3), (12, 3), (20, 3), (10, 1), (10, 5), (10, 7)])
        by_pair = {(s['u'], s['run_length']): s for s in out['sensitivity']}
        self.assertEqual(by_pair[(10, 7)]['cluster_count_accepted'], 4)
        self.assertEqual(by_pair[(20, 3)]['cluster_count_accepted'], 6)

    def test_report_is_strict_json(self):
        for name in ('normal.json', 'boundary-clusters.json',
                     'incomplete-segment.json', 'no-exceedance.json'):
            with self.subTest(name=name):
                out = analyze_file(name, cfg())
                json.dumps(out, allow_nan=False)

    def test_report_identity_and_conventions(self):
        out = analyze_file('normal.json', cfg())
        self.assertEqual(len(out['source_sha256']), 64)
        self.assertEqual(out['units']['rainfall'], 'mm')
        self.assertTrue(any('边界簇' in s for s in out['conventions']))
        self.assertTrue(any('全局最优' in s for s in out['limitations']))


if __name__ == '__main__':
    unittest.main()
