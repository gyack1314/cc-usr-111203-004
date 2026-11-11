"""读取气候合成资料，固定原始字节摘要和基本数值口径。"""
import datetime
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

def number(x):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x):raise ValueError('需要有限数值')
    return x

def positive(x):
    if number(x)<=0:raise ValueError('需要正数')
    return x

def integer(x,minimum=0):
    if type(x) is not int or x<minimum:raise ValueError('无效整数')
    return x

def date(x):
    if not isinstance(x,str):raise ValueError('日期使用YYYY-MM-DD')
    try:d=datetime.date.fromisoformat(x)
    except ValueError:raise ValueError('公历日期无效') from None
    if d.isoformat()!=x:raise ValueError('日期格式必须固定')
    return d

def fraction(x):
    if isinstance(x,bool) or not isinstance(x,(str,int)):raise ValueError('精确量使用整数或分数字符串')
    try:return Fraction(x)
    except (ValueError,ZeroDivisionError):raise ValueError('有理数无效') from None

def unique(rows,key):
    seen=set()
    for row in rows:
        x=row[key]
        if not isinstance(x,str) or not x or x in seen:raise ValueError('身份应非空且唯一')
        seen.add(x)
    return seen

def load(path):
    raw=Path(path).read_bytes()
    v=json.loads(raw.decode('utf8'),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('JSON不得含NaN或Infinity')))
    if not isinstance(v,dict):raise ValueError('顶层资料应为对象')
    return v,hashlib.sha256(raw).hexdigest()
