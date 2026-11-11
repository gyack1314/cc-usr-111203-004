"""显示已有原始摘要，来源文件保持不动。"""
import argparse
import json
from .io import load
from .domain import summarize

def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['inspect']);parser.add_argument('file');a=parser.parse_args()
    v,digest=load(a.file)
    print(json.dumps({'source_sha256':digest,'summary':summarize(v)},ensure_ascii=False,indent=2,allow_nan=False))
if __name__=='__main__':main()
