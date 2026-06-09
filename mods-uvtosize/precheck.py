"""
提交前自检：扫描 algo.py 的安全红线 + 检查 sample/ 是否放了测试样例。
用法：python precheck.py    （全绿才能提交）
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
FORBIDDEN = [
    (r"\bos\.system\b", "禁止 os.system"),
    (r"\bsubprocess\b", "禁止 subprocess"),
    (r"\beval\s*\(", "禁止 eval("),
    (r"\bexec\s*\(", "禁止 exec("),
    (r"[A-Za-z]:\\\\?Users", "禁止写死本地绝对路径 (C:\\Users...)"),
    (r"/home/[a-z]", "禁止写死本地绝对路径 (/home/...)"),
    (r"\bshutil\.rmtree\b|\bos\.remove\b|\bos\.unlink\b", "禁止删除文件操作"),
]

def scan_file(path):
    issues = []
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            for pat, msg in FORBIDDEN:
                if re.search(pat, line):
                    issues.append(f"  algo.py:{i}  {msg}  →  {line.strip()[:60]}")
    return issues

def main():
    problems = []

    algo = os.path.join(ROOT, "algo.py")
    if not os.path.exists(algo):
        problems.append("  缺少 algo.py")
    else:
        problems += scan_file(algo)

    sample_dir = os.path.join(ROOT, "sample")
    samples = [f for f in os.listdir(sample_dir)] if os.path.isdir(sample_dir) else []
    samples = [f for f in samples if not f.startswith(".") and f.lower() != "readme.txt"]
    if not samples:
        problems.append("  sample/ 里没有放测试样例文件（必须放一个真实输入）")

    print("=" * 50)
    if problems:
        print("[FAIL] precheck 未通过：\n" + "\n".join(problems))
        print("=" * 50)
        sys.exit(1)
    print("[OK] precheck 全绿，可以打包提交了！")
    print("=" * 50)

if __name__ == "__main__":
    main()
