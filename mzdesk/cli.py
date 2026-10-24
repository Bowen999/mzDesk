"""mzdesk 启动命令。

用法：
  mzdesk                 启动，默认 http://localhost:8000
  mzdesk --port 8599     指定端口
  mzdesk --host 0.0.0.0  监听所有地址（Docker）
  mzdesk --reload        开发模式，改后端代码自动重启

前端在 web/，npm run build 后输出到 mzdesk/static/，由后端托管。
"""
import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mzdesk", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--reload", action="store_true")
    a = ap.parse_args(argv)
    import uvicorn
    uvicorn.run("mzdesk.api.app:app", host=a.host, port=a.port, reload=a.reload)
    return 0


if __name__ == "__main__":
    sys.exit(main())
