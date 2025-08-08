import subprocess
import datetime
import re

def log(msg):
    """打印带时间戳的日志"""
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}")


def parse_output(output):
    """
    从 RNAinverse 输出中提取 d 值和序列。
    如果没有 d= 行，说明 difference = 0（完全匹配）。
    """
    match = re.search(r"d=\s*(\d+)", output)
    sequence = output.strip().split()[0] if output.strip() else None
    if match:
        d_value = int(match.group(1))
        return d_value, sequence
    return None, sequence  # 视为 d=0


def run_rna_inverse(input_file="inv.in", max_attempts=1000):
    found = False

    for attempt in range(1, max_attempts + 1):
        process = subprocess.Popen(
            ['RNAinverse'],
            stdin=open(input_file, 'r'),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            log(f"❌ Error (attempt {attempt}): {stderr.strip()}")
            continue

        d_value, sequence = parse_output(stdout)

        if d_value is None or d_value == 0:
            # 当没有 d= 行时，视为 d=0（完全匹配）
            log(f"✅ Found PERFECT MATCH (d=0) at attempt {attempt}:\nSequence: {sequence}\nFull Output:\n{stdout.strip()}")
            found = True
            break
        elif d_value < 5:
            log(f"🔍 Found d={d_value} at attempt {attempt}: {sequence}")

    if not found:
        log(f"❌ No d=0 result found in {max_attempts} attempts.")

if __name__ == "__main__":
    run_rna_inverse()
