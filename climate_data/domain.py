"""检查雨量原始值并列出声明观测段的跨度。"""
from .io import date,number,unique

def summarize(v):
    unique(v['segments'],'id');spans=[];prior=None
    for s in sorted(v['segments'],key=lambda x:x['start']):
        lo=date(s['start']);hi=date(s['end'])
        if hi<lo or prior is not None and lo<=prior:raise ValueError('观测段重叠或逆序')
        prior=hi;spans.append({'id':s['id'],'calendar_days':(hi-lo).days+1})
    seen=set();observed=0;total=0
    for r in v['rain']:
        d=date(r['day'])
        if d in seen:raise ValueError('雨量日期重复')
        seen.add(d)
        if r['mm'] is not None:
            x=number(r['mm'])
            if x<0:raise ValueError('雨量不能负')
            observed+=1;total+=x
    return {'declared_segments':spans,'raw_observed_days':observed,'raw_total_mm':total}
