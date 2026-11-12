"""显示已有原始摘要或运行超阈值尾部分析，来源文件保持不动。"""
import argparse
import json
from .io import load
from .domain import summarize
from .tail import analyze

def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['inspect','tail']);parser.add_argument('file');parser.add_argument('config',nargs='?');a=parser.parse_args()
    v,digest=load(a.file)
    if a.command=='inspect':
        print(json.dumps({'source_sha256':digest,'summary':summarize(v)},ensure_ascii=False,indent=2,allow_nan=False))
        return
    if a.config is None:parser.error('tail需要分析配置文件')
    cfg,cfg_digest=load(a.config)
    print(json.dumps(analyze(v,digest,cfg,cfg_digest),ensure_ascii=False,indent=2,allow_nan=False))
if __name__=='__main__':main()
