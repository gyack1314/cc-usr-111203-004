"""核对原始资料读取和已有摘要的边界。"""
import json
import tempfile
import unittest
from pathlib import Path
from climate_data.io import date,number,integer,fraction,load,unique
from climate_data.domain import summarize

class MaterialTests(unittest.TestCase):
    def test_normal(self):
        original,digest=load('data/normal.json');value=summarize(original)
        self.assertEqual(value,summarize(json.loads(json.dumps(original))))
        self.assertEqual(len(digest),64);self.assertEqual(digest,load('data/normal.json')[1])
        self.assertEqual(value['raw_observed_days'],90);self.assertEqual(value['declared_segments'][0]['calendar_days'],90)
    def test_domain_rejections(self):
        for p in Path('data').glob('invalid-*.json'):
            with self.subTest(file=p.name):
                with self.assertRaises((ValueError,KeyError,TypeError)):summarize(load(p)[0])
    def test_additional_material(self):
        for p in Path('data').glob('*.json'):
            if p.name.startswith('invalid-'):continue
            with self.subTest(file=p.name):self.assertIsNotNone(summarize(load(p)[0]))
    def test_numeric_rejection(self):
        for x in (True,float('inf'),float('nan'),'2'):
            with self.assertRaises(ValueError):number(x)
        for x in (True,1.5,-1):
            with self.assertRaises(ValueError):integer(x)
        self.assertEqual(str(fraction('2/4')),'1/2')
        for x in ('1/0',True,0.1):
            with self.assertRaises(ValueError):fraction(x)
    def test_calendar_contract(self):
        self.assertEqual(date('2020-02-29').day,29)
        for x in ('2019-02-29','20200101','2020-2-01',None):
            with self.assertRaises(ValueError):date(x)
    def test_identity(self):
        with self.assertRaises(ValueError):unique([{'id':'same'},{'id':'same'}],'id')
    def test_json_contract(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'资料.json';p.write_text('{"海温":2}',encoding='utf8')
            self.assertEqual(load(p)[0],{'海温':2})
            for raw in ('[]','{"a":NaN}','{"a":Infinity}'):
                p.write_text(raw,encoding='utf8')
                with self.assertRaises(ValueError):load(p)
